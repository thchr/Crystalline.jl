const MULTIPLICITY_ATOL = 2e-2 # allow up to 2% error in irrep multiplicity by default
const MAXRESNORM_TOL = 1e-2    # allow up to 1% error in residual norm of symeigs vs. irreps

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Given a vector of symmetry eigenvalue data `symeigsv` and an associated set of band
representations, `brs`, return a set of bands, as a `Vector{SymmetryVector{D}}`, each
element of which is a minimum-occupation and compatibility-respecting band grouping. The
lowest-lying bands are returned first.

## Input arguments

- `symeigsv`: a triply-nested vector of symmetry eigenvalues, with the following indexing
  convention that `symeigsv[kidx][bandidx][opidx]` gives the symmetry eigenvalue of the
  `kidx`th **k**-point and the `bandidx`th band, under the action of the `opidx`th symmetry
   operation in the little group of the `kidx`th **k**-point. The sorting of little group
   operations must correspond to those in `group(irreps(brs)[kidx])`.
- `brs :: Collection{NewBandRep{D}}`: a collection of band representations, iterating a set
   of `NewBandRep{D}` objects, obtained from [`calc_bandreps`](@ref), and is expected to be
   provided in `primitivized` form (see [`primitivize(::Collection{<:NewBandRep})`](@ref)).
   The little group irreps are implicitly specified via `brs` as well (via `irreps(brs)`),
   as are the corresponding little groups and their associated operator sorting (via
   `group.(irreps(brs))`). It assumed that the sorting of symmetry eigenvalues in `symeigsv`
   is such that `group(irreps(brs)[kidx])[opidx]` matches the elements of 
   `symeigsv[kidx][:][opidx]`.
- `F::Smith{<:Integer}`: optional argument, corresponding to the Smith normal form of the
   band representation matrix `stack(brs)`. Can be supplied explicitly to avoid repeated
   recalculations for repeated calls to this function (see `smith`).

## Keyword arguments

Keyword arguments `kws` are forwarded to [`find_multiplicities`](@ref), which determines the
irrep multiplicities from the symmetry eigenvalues `symeigsv`.

For low-resolution calculations (i.e., if `symeigsv` has appreciable numerical error),
it can be worthwhile to increase the absolute tolerance `atol` (default, $MULTIPLICITY_ATOL)
used by [`find_multiplicities`](@ref).
"""
function symeigs_analysis(
    symeigsv::AbstractVector{<:AbstractVector{<:AbstractVector{<:Number}}},
    brs::Collection{NewBandRep{D}},
    F::Smith{<:Integer} = smith(stack(brs));
    kws...
) where D
    lgirsv = irreps(brs)
    bandirsv = map(zip(lgirsv, symeigsv)) do (lgirs, symeigs)
        find_multiplicities(symeigs, lgirs; kws...)
    end
    candidate_ns = build_candidate_symmetryvectors(bandirsv, lgirsv)

    ns = SymmetryVector{D}[]
    idx = 1
    while idx ≤ length(candidate_ns)
        n_and_idx = _find_next_separable_band_grouping(candidate_ns, F, idx)
        if !isnothing(n_and_idx)
            # found a separable band grouping with sym vec `n′ = sum(ns[idx:idx′])`
            n, idx′ = n_and_idx
            push!(ns, n)
            idx = idx′ + 1 # set to next "starting" index
        else
            break # could not find any more separable bands in `ns`; stop iteration
        end
    end
    return ns
end

function _find_next_separable_band_grouping(
    ns::AbstractVector{<:AbstractVector{<:Integer}}, F::Smith, idx::Integer=1
)
    idx > length(ns) && return nothing
    n′ = ns[idx]
    iscompatible(n′, F; allow_negative=true) && return n′, idx
    while (idx += 1) ≤ length(ns)
        n′ += ns[idx]
        iscompatible(n′, F; allow_negative=true) && return n′, idx
    end
    return nothing
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return `ns` with `ns` denoting symmetry vectors associated with a vector of
symmetry data `bandirsv` (see [`build_candidate_symmetryvectors`](@ref)), in the 
rrep-sorting of `lgirsv`. The returned vectors are sorted according to iteration through
`symeigsv[…][bandidx][…]` (i.e., low-lying bands first).

The returned symmetry vectors are _potentially_ separable, in the sense that they have
integer symmetry data and a consistent occupation number across **k**-points. Aside from
this, no account of compatibility relations is made. To test whether a returned symmetry
vector in fact fulfills compatibility relations, see [`iscompatible`](@ref).

If no potentially separable groupings can be found, an empty pair of vectors are returned.

## Input
The required irrep data `bandirsv` can be obtained from
[`build_candidate_symmetryvectors`](@ref).


## Keyword arguments
- `latestarts` (default, `nothing`): a dictionary of **k**-point labels and a late-start
  band index. Bands below this index at the corresponding **k**-label are not included In
  the tallying of symmetry vectors. This is useful for skipping bands whose symmetry data
  may be known to be incorrect (such as the zero-frequency data of photonic crystals).
"""
function build_candidate_symmetryvectors(
    bandirsv::Vector{Vector{Pair{UnitRange{Int}, Vector{Int}}}},
    lgirsv::AbstractVector{Collection{LGIrrep{D}}};
    latestarts::Union{Dict{String, Int}, Nothing} = nothing
) where D
    if length(bandirsv) ≠ length(lgirsv)
        throw(DimensionMismatch("incompatible lengths of `bandirsv` and `lgirsv`"))
    end
    
    # Check for empty vectors in bandirsv and return gracefully. Otherwise the mapreduce
    # line will throw an error.
    if any(isempty, bandirsv)
        # return empty `collectibles_bands, collectibles_symvecs`
        return Vector{SymmetryVector{D}}()
    end

    Nbands = mapreduce(last∘first∘last, min, bandirsv) # smallest "last" band-index
    include = [Int[] for _ in eachindex(bandirsv)]
    occupations = Vector{Int}()
    includesv = Vector{typeof(include)}()
    start = stop = 1
    while stop ≤ Nbands
        stable = true
        for (kidx, bandirs) in enumerate(bandirsv)
            idxs = include[kidx]
            latestart = isnothing(latestarts) ? nothing : get(latestarts, klabel(lgirsv[kidx]), nothing)
            for (i, (bands, _)) in enumerate(bandirs)
                minband, maxband = extrema(bands)
                minband < start && continue

                if latestart !== nothing && stop < latestart && minband ≤ latestart
                    break # allow `bands` to "not count" if `latestarts` indicates skips
                end
                i ∉ idxs && push!(idxs, i)
                
                if maxband > stop
                    stop = maxband
                    stable = false
                    break
                elseif maxband == stop
                    break
                end
            end
        end

        if stable
            push!(occupations, stop-start+1)
            push!(includesv, include)
            start = (stop += 1)
            include = [Int[] for _ in eachindex(bandirsv)]
        end
    end
    
    # create projected symmetry vectors for each, using knowledge of number of irreps in
    # each little group
    multsv = Vector{Vector{Vector{Int}}}()
    for include in includesv
        mults = [zeros(Int, length(lgirs)) for lgirs in lgirsv]
        for (kidx, idxs) in enumerate(include)
            bandirs = bandirsv[kidx]
            for i in idxs
                mults[kidx] += last(bandirs[i])
            end
        end
        push!(multsv, mults)
    end

    # construct symmetry vectors from multiplicity-data
    ns = map(zip(occupations, multsv)) do (μ, mults)
        SymmetryVector(lgirsv, mults, μ)
    end

    return ns
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return a `Vector{Pair{UnitRange{Int}, Vector{Int}}}` of band-position indices (`UnitRange`)
and associated irrep multiplicities (`Vector{Int}`) at a single **k**-point associated with
the provided symmetry eigenvalues `symeigs` and little group irreps `lgirs`. Each element
of the vector contains as few bands as needed to obtain integer multiplicities (and iterates
from lowest to highest bands in `symeigs`).

Indexing into `symeigs` is assumed to be such that `symeigs[bandidx][opidx]` returns the
the symmetry eigenvalue of `bandidx`th energy band and the `opidx`th little group operation,
corresponding to `group(lgirs)[opidx]`.

## Keyword arguments
- `atol` (default, `$MULTIPLICITY_ATOL`): absolute tolerance used in computing the irrep
  multiplicities. Passed to [`find_representation`](@ref) and also used as a maximum
  allowable deviation from computed floating point multiplicities to the associated nearest
  integer values.
- `αβγ` (default, `$TEST_αβγ`): fractional parameters provided to any little group irreps 
  with nonspecial **k**-points.
- `latestarts` (default, `nothing`): see description in
  [`build_candidate_symmetryvectors`](@ref).
- `maxresnorm` (default, `$MAXRESNORM_TOL`): passsed to [`find_representation`](@ref).
   Denotes the maximum allowable residual norm-difference between provided symmetry
   eigenvalues and the symmetry eigenvalues associated with the computed floating point
   multiplicities.
"""
function find_multiplicities(
    symeigs::AbstractVector{<:AbstractVector{<:Number}},
    lgirs::Collection{LGIrrep{D}};
    atol::Real = MULTIPLICITY_ATOL,
    αβγ::AbstractVector{<:Real} = TEST_αβγ,
    latestarts::Union{Dict{String,Int}, Nothing} = nothing,
    maxresnorm::Real = MAXRESNORM_TOL
) where D
    
    Nbands = length(symeigs)
    bandirs = Pair{UnitRange{Int}, Vector{Int}}[]
    start = stop = isnothing(latestarts) ? 1 : get(latestarts, klabel(lgirs), 1)
    while stop ≤ Nbands
        bands = start:stop
        n = find_representation(symeigs, lgirs, bands, Float64; atol, αβγ, maxresnorm)
        if !isnothing(n)
            # check if at least one multiplicity is a nonzero near-integer (up to
            # `atol`) and all other multiplcities are near-zero (also up to `atol`);
            # equivalently, whether `n` is nearly an integer vector, and has at least
            # one nonzero element
            idxs = findall(nᵢ -> abs(nᵢ) > atol, n) # nonzero elements (allow negative)
            if (!isempty(idxs) && # at least one nonzero multiplicity
                all(nᵢ -> abs(round(nᵢ) - nᵢ) < atol, n)) # all mults. are near-integer
                # `bands` makes up a valid whole-multiple irrep combination at
                # `klabel(lgirs)`; note that this can be a _combination_ i.e. multiple whole
                # irreps; this can arise in cases of near-degeneracies, where numerical
                # "mode-mixing" effectively spreads two irreps over two bands bands
                # (even if they _should_ be each in one band); in principle, this can
                # always be circumvented by increasing the resolution, but this is not
                # a practical solution in general
                push!(bandirs, bands => round.(Int, n)::Vector{Int})
                start = stop + 1 # prepare for next band grouping
            end
        end
        stop += 1
    end

    return bandirs
end