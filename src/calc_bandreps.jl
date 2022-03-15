using LinearAlgebra
using Crystalline
using Crystalline: irdim, constant, free, AbstractIrrep, iscorep,
                   _mulliken, DEFAULT_ATOL, formatirreplabel
import Crystalline: mulliken, realify, group
using StaticArrays

# The implementation here follows Elcoro et al., Phys. Rev. B 97, 035139 (2018)
# (https://doi.org/10.1103/PhysRevB.97.035139), specifically, Sections II.C-D

# ---------------------------------------------------------------------------------------- #
# Site symmetry irreps

struct SiteIrrep{D} <: AbstractIrrep{D}
    cdml     :: String
    g        :: SiteGroup{D}
    matrices :: Vector{Matrix{ComplexF64}}
    reality  :: Reality
    iscorep  :: Bool
    pglab    :: String
end
position(siteir::SiteIrrep) = position(group(siteir))

"""
    siteirreps(sitegroup::SiteGroup) --> ::Vector{SiteIrrep}

Return the site symmetry irreps associated with the provided `SiteGroup`, derived from a
search over isomorphic point groups.
"""
function siteirreps(siteg::SiteGroup{D}) where D
    parent_pg, Iᵖ²ᵍ, _ = find_isomorphic_parent_pointgroup(siteg)
    pglab = label(parent_pg)
    pgirs = pgirreps(pglab, Val(D))
    
    # note that we _have to_ make a copy when re-indexing `pgir.matrices` here, since
    # .jld files apparently cache accessed content; so if we modify it, we mess with the
    # underlying data (see https://github.com/JuliaIO/JLD2.jl/issues/277)
    siteirs = map(pgirs) do pgir
        SiteIrrep{D}(label(pgir), siteg, pgir.matrices[Iᵖ²ᵍ], reality(pgir), pgir.iscorep,
                     pglab)
    end
    return siteirs
end
pglabel(siteir::SiteIrrep)  = siteir.pglab # associated point group label
mulliken(siteir::SiteIrrep) = _mulliken(pglabel(siteir), label(siteir), iscorep(siteir))

# ---------------------------------------------------------------------------------------- #
# Misc utility functions 

function realify!(lgirsd::Dict{String, <:Vector{<:LGIrrep}})
    for (klab, lgirs) in lgirsd
        lgirsd[klab] = realify(lgirs)
    end
    return lgirs
end

"""
    reduce_dict_of_vectors([f=identity,] d::Dict{_, <:Vector})  -->  Vector{T}

Return the concatenated vector of all of vectors in `d` under the element-wise application
of `f`. Effectively flattens a `Dict` of `Vector`s to a single `Vector`.

Note that application of `f` to vector-elements of `d` must return a stable type `T`.
"""
function reduce_dict_of_vectors(f::F, d::Dict{<:Any, <:Vector}) where F
    # get element type of application of `f` to elements of vectors in `d`; assumed fixed
    eltype = typeof(f(first(first(values(d)))))
    # get total number of elements across vectors in `d`
    N = sum(((_, v),) -> length(v), d)
    # preallocate output vector
    w      = Vector{eltype}(undef, N)
    start  = 1
    for v in values(d)
        stop = start + length(v) - 1
        @inbounds for (i, j) in enumerate(start:stop)
            w[j] = f(v[i])
        end
        start = stop + 1
    end
    return w
end
reduce_dict_of_vectors(d::Dict{<:Any, <:Vector}) = reduce_dict_of_vectors(identity, d)

"""
    reduce_orbits!(orbits::Vector{WyckoffPosition}, cntr::Char, conv_or_prim::Bool=true])

Update `orbits` in-place to contain only those Wyckoff positions that are not equivalent
in a primitive basis (as determined by the centering type `cntr`). Returns the updated
`orbits`. 

If `conv_or_prim = true` (default), the Wyckoff positions are returned in the original,
conventional basis; if `conv_or_prim = false`, they are returned in a primitive basis.
"""
function reduce_orbits!(orbits::Vector{WyckoffPosition{D}}, cntr::Char,
            conv_or_prim::Bool=true) where D

    orbits′ = primitivize.(orbits, cntr)
    i = 2 # start from second element
    while i ≤ length(orbits)
        wp′ = parent(orbits′[i])
        was_removed = false
        for j in 1:i-1
            δᵢⱼ = wp′ - parent(orbits′[j])
            if (iszero(free(δᵢⱼ)) && 
                all(x->abs(rem(x, 1.0, RoundNearest)) < DEFAULT_ATOL, constant(δᵢⱼ)))
                # -> means it's equivalent to a previously "greenlit" orbit
                deleteat!(orbits, i)
                deleteat!(orbits′, i)
                was_removed = true
                break
            end
        end
        was_removed || (i += 1)
    end

    return conv_or_prim ? orbits : orbits′
end

# TODO: This is probably pretty fragile and depends on an implicit shared iteration order of
#       cosets and orbits; it also assumes that wp is the first position in the `orbits`
#       vector. Seems to be sufficient for now though...
function reduce_cosets!(ops::Vector{SymOperation{D}}, wp::WyckoffPosition{D}, 
            orbits::Vector{WyckoffPosition{D}}) where D
    wp ≈ first(orbits) || error("implementation depends on a shared iteration order; sorry...")
    i = 1
    while i ≤ length(ops) && i ≤ length(orbits)
        wpᵢ = orbits[i]
        opᵢ = ops[i]
        if ops[i]*parent(wp) ≈ wpᵢ
            i += 1 # then ops[i] is indeed a "generator" of wpᵢ
        else
            deleteat!(ops, i)
        end
    end
    i ≤ length(ops) && deleteat!(ops, i:length(ops))
    return ops
end

# ---------------------------------------------------------------------------------------- #
# Bandrep related functions: induction/subduction

"""
    induce_bandrep(siteir::SiteIrrep, h::SymOperation, kv::KVec)

Return the band representation induced by the provided `SiteIrrep` evaluated at `kv` and for
a `SymOperation` `h`.
"""
function induce_bandrep(siteir::SiteIrrep{D}, h::SymOperation{D}, kv::KVec{D}) where D

    siteg  = group(siteir)
    wp     = position(siteg)   
    kv′    = constant(h*kv) # <-- TODO: Why only constant part?
    χs     = characters(siteir)

    # only loop over the cosets/orbits that are not mere centering "copies"
    orbits = orbit(siteg)
    gαs    = cosets(siteg)
    cntr   = centering(num(siteg), D)
    reduce_orbits!(orbits, cntr, true) # remove centering "copies" from orbits
    reduce_cosets!(gαs, wp, orbits)    # remove the associated cosets

    # sum over all the (non-centering-equivalent) wyckoff positions/cosets in the orbit 
    χᴳₖ = zero(ComplexF64)
    for (wpα′, gα′) in zip(orbits, gαs)
        tα′α′ = constant(parent(h*wpα′) - parent(wpα′)) # TODO: <-- explain why we only need constant part here?
        opᵗ   = SymOperation(-tα′α′)

        gα′⁻¹     = inv(gα′)
        gα′⁻¹ggα′ = compose(gα′⁻¹, compose(opᵗ, compose(h, gα′, false), false ), false)

        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), siteg)
        if site_symmetry_index !== nothing
            χᴳₖ += cis(2π*dot(kv′, tα′α′)) * χs[site_symmetry_index]
            # TODO: The sign in this `cis(...)` above is different from in Elcoro's. Why?
        end
    end
    return χᴳₖ
end

function subduce_onto_lgirreps(siteir::SiteIrrep{D}, lgirs::Vector{LGIrrep{D}}) where D
    lg = group(first(lgirs))
    kv = position(lg)

    # characters of induced site representation and little group irreps
    site_χs  = induce_bandrep.(Ref(siteir), lg, Ref(kv))
    lgirs_χm = matrix(characters(lgirs))

    # little group irrep multiplicities, after subduction
    m  = lgirs_χm\site_χs
    m′ = round.(Int, real.(m)) # truncate to integers
    isapprox(m, m′, atol=DEFAULT_ATOL) || error(DomainError(m, "failed to convert to integers"))
    return m′
end

function calc_bandrep(siteir::SiteIrrep{D}, lgirsd::Dict{String, Vector{LGIrrep{D}}};
            irlabs::Vector{String} = reduce_dict_of_vectors(formatirreplabel∘label, lgirsd),
            irdims::Vector{Int}    = reduce_dict_of_vectors(irdim, lgirsd)) where D

    irvec = Int[]
    for (_, lgirs) in lgirsd
        kv_array = subduce_onto_lgirreps(siteir, lgirs)
        append!(irvec, Int.(kv_array))
    end    
    
    irdims_Γmask = [irlab[1]=='Γ' for irlab in irlabs] .* irdims # masked w/ Γ-irreps only
    brdim   = dot(irvec, irdims_Γmask)
    wycklab = label(position(group(siteir)))
    spinful      = false # NB: default; Crystalline currently doesn't have spinful irreps
    decomposable = false # NB: placeholder, because we don't know the true answer presently

    return BandRep(wycklab, pglabel(siteir), mulliken(siteir)*"↑G", brdim, decomposable,
                   spinful, irvec, irlabs)
end

# ---------------------------------------------------------------------------------------- #
"
    calc_bandreps(sgnum::Integer, Dᵛ::Val{D}=Val(3);
                  timereversal::Bool=true,
                  allpaths::Bool=false)

Compute the band representations of space group `sgnum` in dimension `D`, returning a
`BandRepSet`.

## Keyword arguments
- `timereversal` (default, `true`): whether the irreps used to induce the band
  representations are assumed to be time-reversal invariant (i.e., are coreps, see 
  [`realify`](@ref)).
- `allpaths` (default, `false`): whether the band representations are projected to all
  distinct **k**-points returned by `lgirreps` (`allpaths = false`), including high-symmetry
  **k**-lines and -plane, or only to the maximal **k**-points (`allpaths = true`), i.e.,
  just to high-symmetry points.

## implementation
The implementation is based on Elcoro et al., [Phys. Rev. B 97, 035139
(2018)](https://doi.org/10.1103/PhysRevB.97.035139), Sections II.C-D.
"
function calc_bandreps(sgnum::Integer, Dᵛ::Val{D}=Val(3);
                       timereversal::Bool=true,
                       allpaths::Bool=false) where D

    # get all the little group irreps that we want to subduce onto
    lgirsd = lgirreps(sgnum, Val(D))
    allpaths     || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && realify!(lgirsd)

    irlabs = reduce_dict_of_vectors(formatirreplabel∘label, lgirsd)
    irdims = reduce_dict_of_vectors(irdim, lgirsd)

    # get the bandreps induced by every maximal site symmetry irrep
    wps    = wyckoffs(sgnum, Dᵛ)
    sg     = spacegroup(sgnum, Dᵛ)
    sitegs = findmaximal(SiteGroup.(Ref(sg), wps))
    brs    = BandRep[]
    for siteg in sitegs
        siteirs = siteirreps(siteg)
        timereversal && (siteirs = realify(siteirs))
        append!(brs, calc_bandrep.(siteirs, Ref(lgirsd); irlabs=irlabs, irdims=irdims))
    end

    # TODO: There's potentially some kind of bad implicit dependence/overreliance on fixed
    #       on iteration order of `Dict`s here and related methods; probably OK, but unwise
    klabs  = collect(keys(lgirsd))
    kvs    = position.(first.(values(lgirsd)))
    spinful      = false # NB: default; Crystalline currently doesn't have spinful irreps
    decomposable = false # NB: placeholder, because we don't know the true answer presently

    return BandRepSet(sgnum, brs, kvs, klabs, irlabs, allpaths, spinful, timereversal)
end