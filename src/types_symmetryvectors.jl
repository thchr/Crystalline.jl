# ---------------------------------------------------------------------------------------- #
# AbstractSymmetryVector definition
abstract type AbstractSymmetryVector{D} <: AbstractVector{Int} end

# ---------------------------------------------------------------------------------------- #
# SymmetryVector

"""
    SymmetryVector{D} <: AbstractSymmetryVector{D}

A symmetry vector in dimension `D`, containing the featured irreps and their multiplicities
and overall band occupation number.

## Fields
- `lgirsv :: Vector{Collection{LGIrrep{D}}}` (`const`): a vector of `LGIrrep{D}` collections
  associated with each high-symmetry **k**-point of a space group.
- `multsv :: Vector{Vector{Int}}` (`const`): a vector of associated irrep multiplicities.
  The irrep `lgirsv[i][j]` occurs with multiplicity `multsv[i][j]` in the symmetry vector.
- `occupation :: Int`: the occupation (number of bands) associated with the symmetry vector.

## Interface
- The irreps, multiplicities, and band occupation of a `SymmetryVector` can be obtained via
  `irreps`, `multiplicities`, and `occupation`, respectively.
- The irrep-labels, **k**-labels, space group number of a `SymmetryVector` are available via
  `irreplabels`, `klabels`, and `num`, respectively.
- The `SymmetryVector` can be converted to a "raw" concatenated vector representation via
  `Vector`.

## Construction
- From strings: see `parse(::Type{SymmetryVector{D}}, ::AbstractString,
  ::Vector{Collection{LGIrrep{D}}})`.
- From "raw" concatenated vectors: see [`SymmetryVector`](@ref)`(::AbstractVector{<:Integer},
  ...)`.
"""
mutable struct SymmetryVector{D} <: AbstractSymmetryVector{D}
    const lgirsv :: Vector{Collection{LGIrrep{D}}}
    const multsv :: Vector{Vector{Int}}
    occupation   :: Int
end

# ::: AbstractSymmetryVector interface :::
irreps(n::SymmetryVector) = n.lgirsv
multiplicities(n::SymmetryVector) = n.multsv
occupation(n::SymmetryVector) = n.occupation
SymmetryVector(n::SymmetryVector) = n
SymmetryVector{D}(n::SymmetryVector{D}) where D = n
SymmetryVector{D′}(::SymmetryVector{D}) where {D′, D} = error("incompatible dimensions")

# ::: AbstractArray interface beyond AbstractSymmetryVector :::
function Base.similar(n::SymmetryVector{D}) where D
    SymmetryVector{D}(irreps(n), similar(multiplicities(n)), 0)
end
@propagate_inbounds function Base.setindex!(n::SymmetryVector, v::Int, i::Int)
    Nⁱʳ = length(n)
    @boundscheck i > Nⁱʳ && throw(BoundsError(n, i))
    i == Nⁱʳ && (n.occupation = v; return v)
    i₀ = i₁ = 0
    for mults in multiplicities(n)
        i₁ += length(mults)
        if i ≤ i₁
            mults[i-i₀] = v
            return v
        end
        i₀ = i₁
    end
end

# copy: want to copy just the multiplicities, but not the underlying irrep data
function Base.copy(n::SymmetryVector{D}) where D
    SymmetryVector{D}(n.lgirsv, [copy(mults) for mults in multiplicities(n)], n.occupation)
end

# ::: Optimizations and utilities :::
function Base.Vector(n::SymmetryVector)
    nv = Vector{Int}(undef, length(n))
    i = 1
    @inbounds for mults in multiplicities(n)
        N = length(mults)
        nv[i:i+N-1] .= mults
        i += N
    end
    nv[end] = occupation(n)
    return nv
end
irreplabels(n::SymmetryVector) = [label(ir) for ir in Iterators.flatten(irreps(n))]
klabels(n::SymmetryVector) = [klabel(first(irs)) for irs in irreps(n)]
num(n::SymmetryVector) = num(first(first(irreps(n))))

# ::: Parsing from string :::


""" 
    parse(::Type{SymmetryVector{D}}, 
          s::AbstractString,
          lgirsv::Vector{Collection{LGIrrep{D}}})  ->  SymmetryVector{D}

Parse a string `s` to a `SymmetryVector` over the irreps provided in `lgirsv`. 
The irrep labels of `lgirsv` and `s` must use the same convention.

## Example
```jldoctest
julia> brs = calc_bandreps(220);

julia> lgirsv = irreps(brs); # irreps at P, H, Γ, & PA

julia> s = "[2P₃, 4N₁, H₁H₂+H₄H₅, Γ₁+Γ₂+Γ₄+Γ₅, 2PA₃]";

julia> parse(SymmetryVector, s, lgirsv)
15-irrep SymmetryVector{3}:
 [2P₃, 4N₁, H₁H₂+H₄H₅, Γ₁+Γ₂+Γ₄+Γ₅, 2PA₃] (8 bands)
```
"""
function Base.parse(
            T::Type{<:SymmetryVector},
            s::AbstractString, 
            lgirsv::Vector{Collection{LGIrrep{D}}}) where D
    if !isnothing(dim(T)) && dim(T) != D
        # small dance to allow using both `T=SymmetryVector` & `T=SymmetryVector{D}`
        error("incompatible dimensions of requested SymmetryVector and provided `lgirsv`")
    end
    s′ = replace(s, " "=>"", "*"=>"")
    multsv = [zeros(Int, length(lgirs)) for lgirs in lgirsv]
    μ = typemax(Int) # sentinel for uninitialized
    for (i, lgirs) in enumerate(lgirsv)
        found_irrep = false
        for (j, lgir) in enumerate(lgirs)
            irlab = label(lgir)
            idxs = findfirst(irlab, s′)
            isnothing(idxs) && continue
            found_irrep = true
            pos₂ = first(idxs)
            m = searchpriornumerals(s′, pos₂, Int)
            multsv[i][j] = m
        end
        if !found_irrep
            error("could not identify any irreps associated with the \
                   $(klabel(first(lgirs)))-point in the input string $s")
        end
        μ′ = sum(multsv[i][j]*irdim(lgirs[j]) for j in eachindex(lgirs))
        if μ == typemax(Int) # uninitialized
            μ = μ′
        else
            # validate that occupation number is invariant across all k-points
            μ == μ′ || error("inconsistent occupation numbers across k-points")
        end
    end
    return SymmetryVector(lgirsv, multsv, μ)
end

# ::: Construction from multiplicities, irreps, and Collection{<:LGIrrep} :::
# ---------------------------------------------------------------------------------------- #

"""
    SymmetryVector(
        nv :: AbstractVector{<:Integer},
        irlabs_nv :: AbstractVector{<:AbstractString},
        lgirsd :: AbstractDict{String, <:AbstractVector{LGIrrep{D}}}) --> SymmetryVector{D}

Build a structured `SymmetryVector` representation of a "raw" vector `nv` of irrep
multiplicities, whose `i`th element gives the irrep multiplicity of the irrep whose label 
is `irlabs_nv[i]`. The corresponding full irrep information is inferred by comparison of
labels in `irlabs_nv` and a provided set of covering irreps in `lgirsd`.

The raw vector `nv` must contain its band occupation number as its last element.
The `irlabs_nv` vector must consequently have length `length(nv)-1`.

The sorting of irrep labels in `lgirsd` and (`nv`, `irlabs_nv`) is allowed to differ: this
is the main utility of the function: to map between differently sorted raw vectors and a
structured irrep storage in `lgirsd`.
"""
function SymmetryVector(
            nv::AbstractVector{<:Integer},
            irlabs_nv::AbstractVector{String},
            lgirsd::AbstractDict{String, <:AbstractVector{LGIrrep{D}}}) where D

    klabs = klabel.(unique(klabel, irlabs_nv))
    Nk = length(klabs)
    multsv = [Int[] for _ in 1:Nk]
    lgirsv = [LGIrrep{D}[] for _ in 1:Nk]
    j = 1
    for (nᵢ, irlabᵢ) in zip(nv, irlabs_nv)
        klabᵢ = klabel(irlabᵢ)
        if klabᵢ != klabs[j]
            j += 1
        end
        push!(multsv[j], nᵢ)

        # find associated irrep in `lgirsd[klabᵢ]`
        lgirsⱼ = lgirsd[klabᵢ]
        iridxᵢ = findfirst(lgir -> label(lgir) == irlabᵢ, lgirsⱼ)
        if isnothing(iridxᵢ)
            _throw_failed_to_find_irrep(irlabᵢ, lgirsⱼ)
        else
            push!(lgirsv[j], lgirsⱼ[something(iridxᵢ)])
        end
    end
    lgirsv = [Collection(lgirs) for lgirs in lgirsv]
    if length(nv) ≠ length(irlabs_nv)+1
        error("n must contain its band connectivity")
    end
    μ = nv[end]
    return SymmetryVector{D}(lgirsv, multsv, μ)
end

"""
    SymmetryVectors(
        nvs :: AbstractVector{<:Integer},
        irlabs_nv :: AbstractVector{<:AbstractString},
        lgirsd :: AbstractDict{String, <:AbstractVector{LGIrrep{D}}}) 
                                                            --> Vector{SymmetryVector{D}}

Similar to
[`SymmetryVector(::AbstractVector{<:Integer}, ::AbstractVector{<:AbstractString}, ::AbstractDict)](@ref),
but for a vector of distinct raw multiplicy vectors `nvs`, rather than a single vector,
returning a `Vector{SymmetryVector{D}}`.

The returned `SymmetryVector`s, `ns`, will share the same underlying irrep information such
that `irreps(n) === irreps(n′)` for all `n` and `n′` in `ns`.
"""
function SymmetryVectors(
            nvs::AbstractVector{<:AbstractVector{<:Integer}},
            irlabs_nv::AbstractVector{String},
            lgirsd::AbstractDict{String, <:AbstractVector{LGIrrep{D}}}) where D

    isempty(nvs) && return SymmetryVector{D}[]

    klabs = klabel.(unique(klabel, irlabs_nv))
    Nk, Nir = length(klabs), length(irlabs_nv)
    j = 1
    sortidxs = Vector{Tuple{Int, Int}}(undef, Nir)
    lgirsv = [Collection(LGIrrep{D}[]) for _ in 1:Nk] # to be shared across all `ns`
    max_qs = zeros(Int, Nk)
    for (i, irlabᵢ) in enumerate(irlabs_nv)
        klabᵢ = klabel(irlabᵢ)
        j = @something(findfirst(==(klabᵢ), klabs),
                       error(lazy"failed to associate irrep $irlabᵢ to k-labels $klabs"))
        q = (max_qs[j] += 1)
        sortidxs[i] = (j, q)

        # find associated irrep in `lgirsd[klabᵢ]`
        lgirsⱼ = lgirsd[klabᵢ]
        iridxᵢ = findfirst(lgir -> label(lgir) == irlabᵢ, lgirsⱼ)
        if isnothing(iridxᵢ)
            _throw_failed_to_find_irrep(irlabᵢ, lgirsⱼ)
        else
            push!(parent(lgirsv[j]), lgirsⱼ[something(iridxᵢ)])
        end
    end

    ns = Vector{SymmetryVector{D}}(undef, length(nvs))
    for (r, nv) in enumerate(nvs)
        if length(nv) ≠ length(irlabs_nv)+1
            error("`nv` must contain its band connectivity")
        end
        μ = nv[end]
        multsv = [Vector{Int}(undef, length(lgirs)) for lgirs in lgirsv]
        for i in eachindex(irlabs_nv)
            j, q = sortidxs[i]
            multsv[j][q] = nv[i]
        end
        ns[r] = SymmetryVector(lgirsv, multsv, μ)
    end

    return ns
end

@noinline function _throw_failed_to_find_irrep(irlabᵢ, lgirsⱼ)
    error("failed to find an irrep label \"$irlabᵢ\" in `lgirsd`; available irrep labels \
           at the considered k-point were $(label.(lgirsⱼ))")
end

# ---------------------------------------------------------------------------------------- #
# AbstractSymmetryVector interface & shared implementation

# ::: API :::
"""
    SymmetryVector(n::AbstractSymmetryVector) -> SymmetryVector

Return a `SymmetryVector` realization of `n`. If `n` is already a `SymmetryVector`,
return `n` directly; usually, the returned value will directly reference information in `n`
(i.e., will not be a copy).
"""
function SymmetryVector(::AbstractSymmetryVector) end

# ::: Optional API (if not `SymmetryVector`; otherwise fall-back to `SymmetryVector`) :::
"""
    irreps(n::AbstractSymmetryVector{D}) -> AbstractVector{<:Collection{<:AbstractIrrep{D}}}

Return the irreps referenced by `n`. 

The returned value is an `AbstractVector` of `Collection{<:AbstractIrrep}`s, with irreps for
distinct groups, usually associated with specific **k**-manifolds, belonging to the same
`Collection`.

See also [`multiplicities(::AbstractSymmetryVector)`](@ref).
"""
irreps(n::AbstractSymmetryVector) = irreps(SymmetryVector(n))

"""
    multiplicities(n::AbstractSymmetryVector) -> AbstractVector{<:AbstractVector{Int}}

Return the multiplicities of the irreps referenced by `n`.

See also [`irreps(::AbstractSymmetryVector)`](@ref).
"""
multiplicities(n::AbstractSymmetryVector) = multiplicities(SymmetryVector(n))

"""
    occupation(n::AbstractSymmetryVector) -> Int

Return the occupation of (i.e., number of bands contained within) `n`.
"""
occupation(n::AbstractSymmetryVector) = occupation(SymmetryVector(n))

# misc convenience accessors
irreplabels(n::AbstractSymmetryVector) = irreplabels(SymmetryVector(n))
klabels(n::AbstractSymmetryVector) = klabels(SymmetryVector(n))
num(n::AbstractSymmetryVector) = num(SymmetryVector(n))

# ::: AbstractArray interface :::
Base.size(n::AbstractSymmetryVector) = (mapreduce(length, +, multiplicities(n)) + 1,)
@propagate_inbounds function Base.getindex(n::AbstractSymmetryVector, i::Int)
    Nⁱʳ = length(n)
    @boundscheck i > Nⁱʳ && throw(BoundsError(n, i))
    i == Nⁱʳ && return occupation(n)
    i₀ = i₁ = 0
    for mults in multiplicities(n)
        i₁ += length(mults)
        if i ≤ i₁
            return mults[i-i₀]
        end
        i₀ = i₁
    end
    error("unreachable reached")
end
Base.iterate(n::AbstractSymmetryVector) = multiplicities(n)[1][1], (1, 1)
function Base.iterate(n::AbstractSymmetryVector, state::Tuple{Int, Int})
    i, j = state
    j += 1
    if i > length(multiplicities(n))
        return nothing
    end
    if j > length(multiplicities(n)[i])
        if i == length(multiplicities(n))
            return occupation(n), (i+1, j)
        end
        i += 1
        j = 1
    end
    return multiplicities(n)[i][j], (i, j)
end

# ::: Algebraic operations :::
function Base.:+(n::AbstractSymmetryVector{D}, m::AbstractSymmetryVector{D}) where D
    _n = SymmetryVector(n)
    _m = SymmetryVector(m)
    irreps(_n) === irreps(_m)
    return SymmetryVector(irreps(_n), 
                          multiplicities(_n) .+ multiplicities(_m),
                          occupation(_n) + occupation(_m))
end
function Base.:-(n::AbstractSymmetryVector{D}) where D
    SymmetryVector{D}(irreps(n), -multiplicities(n), -occupation(n))
end
Base.:-(n::AbstractSymmetryVector{D}, m::AbstractSymmetryVector{D}) where D = n + (-m)
function Base.:*(n::AbstractSymmetryVector{D}, k::Integer) where D
    SymmetryVector{D}(irreps(n), [ms .* k for ms in multiplicities(n)], occupation(n) * k)
end
Base.:*(k::Integer, n::AbstractSymmetryVector) = n * k
function Base.zero(n::AbstractSymmetryVector{D}) where D
    SymmetryVector{D}(irreps(n), zero.(multiplicities(n)), 0)
end
# make sure `sum(::AbstractSymmetryVector)` is type-stable (necessary since the + operation
# now may change the type of an `AbstractSymmetryVector` - so `n[1]` and `n[1]+n[2]` may
# be of different types) and always returns a `SymmetryVector`
Base.reduce_first(::typeof(+), n::AbstractSymmetryVector) = SymmetryVector(n)

# ::: Utilities & misc :::
dim(::AbstractSymmetryVector{D}) where D = D
dim(::Type{<:AbstractSymmetryVector{D}}) where D = D
dim(::Type{<:AbstractSymmetryVector}) = nothing

# ---------------------------------------------------------------------------------------- #
# NewBandRep

struct NewBandRep{D} <: AbstractSymmetryVector{D}
    siteir       :: SiteIrrep{D}
    n            :: SymmetryVector{D}
    timereversal :: Bool
    spinful      :: Bool
end

# ::: AbstractSymmetryVector interface :::
SymmetryVector(br::NewBandRep) = br.n

# ::: AbstractArray interface beyond AbstractSymmetryVector :::
Base.setindex!(br::NewBandRep{D}, v::Int, i::Int) where D = (br.n[i] = v)
function Base.similar(br::NewBandRep{D}) where D
    NewBandRep{D}(br.siteir, similar(br.n), br.timereversal, br.spinful)
end
Base.Vector(br::NewBandRep) = Vector(br.n)

# ::: Utilities :::
group(br::NewBandRep) = group(br.siteir)
Base.position(br::NewBandRep) = position(group(br))

# ::: Conversion to BandRep :::
function Base.convert(::Type{BandRep}, br::NewBandRep{D}) where D
    wyckpos     = label(position(br.siteir))
    sitesym     = br.siteir.pglabel
    siteirlabel = label(br.siteir)*"↑G"
    dim         = occupation(br)
    spinful     = br.spinful
    irvec       = collect(br)[1:end-1]
    irlabs      = irreplabels(br)
    return BandRep(wyckpos, sitesym, siteirlabel, dim, spinful, irvec, irlabs)
end

# ---------------------------------------------------------------------------------------- #
# Collection{<:NewBandRep}

# ::: Utilities :::
irreps(brs::Collection{<:NewBandRep}) = irreps(SymmetryVector(first(brs)))
irreplabels(brs::Collection{<:NewBandRep}) = irreplabels(SymmetryVector(first(brs)))
klabels(brs::Collection{<:NewBandRep}) = klabels(SymmetryVector(first(brs)))
littlegroups(::Collection{<:NewBandRep}) = group.(irreps(brs))

# ::: Conversion to BandRepSet :::
function Base.convert(::Type{BandRepSet}, brs::Collection{<:NewBandRep})
    sgnum = num(brs)
    bandreps = convert.(Ref(BandRep), brs)
    kvs = [position(lgirs) for lgirs in irreps(brs)]
    klabs = klabels(brs)
    irlabs = irreplabels(brs)
    spinful = first(brs).spinful
    timereversal = first(brs).timereversal
    return BandRepSet(sgnum, bandreps, kvs, klabs, irlabs, spinful, timereversal)
end


# ---------------------------------------------------------------------------------------- #
# CompositeBandRep

"""
    CompositeBandRep{D} <: AbstractSymmetryVector{D}

A type representing a linear rational-coefficient combination of `NewBandRep{D}`s. 

Although the coefficients may be rational numbers in general, their superposition must
correspond to integer-valued irrep multiplicities and band occupation numbers; in
particular, if they do not, conversion to a `SymmetryVector` will error.

See also [`CompositeBandRep_from_indices`](@ref) for construction from a vector included
indices into `brs`.

## Fields
- `coefs::Vector{Rational{Int}}`: a coefficient vector associated with each band
  representation in `brs`; the coefficient of the `i`th band representation `brs[i]` is
  `coefs[i]`.
- `brs::Collection{NewBandRep{D}}`: the band representations referenced by `coefs`.

## Example
### Fragile symmetry vector
As a first example, we build a `CompositeBandRep` representing a fragilely topological
configuration (i.e., featuring negative integer coefficients):
```julia
julia> brs = calc_bandreps(2);

julia> coefs = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1];

julia> cbr = CompositeBandRep(coefs, brs)
16-irrep CompositeBandRep{3}:
 (1g|Ag) + (1f|Aᵤ) + (1e|Ag) - (1a|Aᵤ) (2 bands)
```
We can build the associated [`SymmetryVector`](@ref) to inspect the associated irrep 
content:
```julia
julia> SymmetryVector(cbr)
 [2Z₁⁺, 2Y₁⁻, 2U₁⁻, 2X₁⁺, 2T₁⁺, 2Γ₁⁺, 2V₁⁺, 2R₁⁺] (2 bands)
```
Similarly, we can confirm that `cbr` is indeed a simple linear combination of the
associated band representations:
```julia
julia> sum(c * brs[i] for (i, c) in enumerate(coefs)) == cbr
true
```
And we can even do simple arithmetic with `cbr` (a `CompositeBandRep` is an element of a
monoid:
```julia
julia> cbr + 2cbr
3(1g|Ag) + 3(1f|Aᵤ) + 3(1e|Ag) - 3(1a|Aᵤ) (6 bands)
```

### Nontrivial symmetry vector
The coefficients of a `CompositeBandRep` can be non-integer, provided the associated irrep
multiplicities of the overall summation are integer. Accordingly, a `CompositeBandRep` may
also be used to represent a topologically nontrivial symmetry vector (and, by extension, any
physical symmetry vector):
```julia
julia> coefs = [-1//4, 0, -1//4, 0, -1//4, 0, 1//4, 0, 1//4, 0, 1//4, 0, -1//4, 0, 1//4, 1];

julia> cbr = CompositeBandRep{3}(coefs, brs)
16-irrep CompositeBandRep{3}:
 -(1/4)×(1h|Ag) - (1/4)×(1g|Ag) - (1/4)×(1f|Ag) + (1/4)×(1e|Ag) + (1/4)×(1d|Ag) + (1/4)×(1c|Ag) - (1/4)×(1b|Ag) + (1/4)×(1a|Ag) + (1a|Aᵤ) (1 band)

julia> SymmetryVector(cbr)
16-irrep SymmetryVector{3}:
 [Z₁⁺, Y₁⁻, U₁⁻, X₁⁻, T₁⁻, Γ₁⁻, V₁⁻, R₁⁻] (1 band)
```
"""
struct CompositeBandRep{D} <: AbstractSymmetryVector{D}
    coefs :: Vector{Rational{Int}}
    brs   :: Collection{NewBandRep{D}}
    function CompositeBandRep{D}(coefs, brs) where D
        if length(coefs) ≠ length(brs)
            error("length of provided coefficients do not match length of provided band \
                   representations")
        end
        new{D}(coefs, brs)
    end
end
function CompositeBandRep(coefs, brs::Collection{NewBandRep{D}}) where D
    return CompositeBandRep{D}(coefs, brs)
end

# ::: Convenience constructor :::
"""
    CompositeBandRep_from_indices(idxs::Vector{Int}, brs::Collection{<:NewBandRep})

Return a [`CompositeBandRep`](@ref) whose symmetry content is equal to the sum of the band
representations in `brs` over `idxs`.
In terms of irrep multiplicity, this is equivalent to `sum(brs[idxs])` in the sense that
`CompositeBandRep(idxs, brs)` is equal to `sum(brs[idxs])` for each irrep multiplicity. 

The difference, and primary motivation for using `CompositeBandRep`, is that  
`CompositeBandRep` retains information about which band representations are included and
with what multiplicity (the multiplicity of the `i`th `brs`-element being equaling to
`count(==(i), idxs)`).

## Example
```julia
julia> brs = calc_bandreps(2);

julia> cbr = CompositeBandRep_from_indices([1, 1, 2, 6], brs)
16-irrep CompositeBandRep{3}:
 (1f|Aᵤ) + (1h|Aᵤ) + 2(1h|Ag) (4 bands)

julia> cbr == brs[1] + brs[1] + brs[2] + brs[6]
true

julia> SymmetryVector(cbr)
16-irrep SymmetryVector{3}:
 [2Z₁⁺+2Z₁⁻, Y₁⁺+3Y₁⁻, 2U₁⁺+2U₁⁻, 2X₁⁺+2X₁⁻, 3T₁⁺+T₁⁻, 2Γ₁⁺+2Γ₁⁻, 3V₁⁺+V₁⁻, R₁⁺+3R₁⁻] (4 bands)
```
"""
function CompositeBandRep_from_indices(idxs::Vector{Int}, brs::Collection{<:NewBandRep})
    coefs = zeros(Rational{Int}, length(brs))
    for i in idxs
        @boundscheck checkbounds(brs, i)
        coefs[i] += 1
    end
    CompositeBandRep(coefs, brs)
end

# ::: AbstractSymmetryVector interface :::
function SymmetryVector(cbr::CompositeBandRep{D}) where {D}
    brs = cbr.brs
    lgirsv = irreps(brs)
    multsv_r = zeros.(eltype(cbr.coefs), size.(multiplicities(first(brs))))
    μ = 0
    for (j, c) in enumerate(cbr.coefs)
        iszero(c) && continue
        multsv_r .+= c * multiplicities(brs[j])
        μ += c * occupation(brs[j])
    end
    multsv = similar.(multsv_r, Int)
    for (i, (mults_r, mults)) in enumerate(zip(multsv_r, multsv))
        for (j, m_r) in enumerate(mults_r)
            m_int = round(Int, m_r)
            if m_int ≠ m_r
                error(LazyString("CompositeBandRep has non-integer multiplicity ",
                                 m_r, " for irrep ", label(lgirsv[i][j]), 
                                 ": cannot convert to SymmetryVector"))
            end
            mults[j] = m_int
        end
    end
    multsv = [Int.(mults) for mults in multsv_r]
    μ = Int(μ)
    return SymmetryVector{D}(lgirsv, multsv, μ)
end
function occupation(cbr::CompositeBandRep)
    μ_r = sum(c * occupation(cbr.brs[j]) for (j, c) in enumerate(cbr.coefs) if !iszero(c);
              init=zero(eltype(cbr.coefs)))
    isinteger(μ_r) || error(lazy"CompositeBandRep has non-integer occupation (= $μ_r)")
    return Int(μ_r)
end

# ::: overloads to avoid materializing a SymmetryVector unnecessarily :::
irreps(cbr::CompositeBandRep) = irreps(first(cbr.brs))
klabels(cbr::CompositeBandRep) = klabels(first(cbr.brs))
irreplabels(cbr::CompositeBandRep) = irreplabels(first(cbr.brs))
num(cbr::CompositeBandRep) = num(first(cbr.brs))

# ::: AbstractArray interface :::
Base.size(cbr::CompositeBandRep) = size(first(cbr.brs))
function Base.getindex(cbr::CompositeBandRep{D}, i::Int) where {D}
    m_r = sum(c * cbr.brs[j][i] for (j, c) in enumerate(cbr.coefs) if !iszero(c);
              init=zero(eltype(cbr.coefs)))
    isinteger(m_r) || error(lazy"CompositeBandRep has non-integer multiplicity (= $m_r)")
    return Int(m_r)
end

# ::: Arithmetic operations :::
function Base.:+(cbr1::CompositeBandRep{D}, cbr2::CompositeBandRep{D}) where D
    if cbr1.brs !== cbr2.brs
        error("provided CompositeBandReps must reference egal `brs`")
    end
    return CompositeBandRep{D}(cbr1.coefs + cbr2.coefs, cbr1.brs)
end
Base.:-(cbr::CompositeBandRep{D}) where D = CompositeBandRep{D}(-cbr.coefs, cbr.brs)
Base.:-(cbr1::CompositeBandRep{D}, cbr2::CompositeBandRep{D}) where D = cbr1 + (-cbr2)
function Base.:*(cbr::CompositeBandRep{D}, n::Integer) where D
    return CompositeBandRep{D}(cbr.coefs .* n, cbr.brs)
end
Base.zero(cbr::CompositeBandRep{D}) where D = CompositeBandRep{D}(zero(cbr.coefs), cbr.brs)
