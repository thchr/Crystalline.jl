# ---------------------------------------------------------------------------------------- #
# AbstractSymmetryVector definition
abstract type AbstractSymmetryVector{D} <: AbstractVector{Int} end

# ---------------------------------------------------------------------------------------- #
# SymmetryVector

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
    SymmetryVector{D}(n.lgirsv, similar(n.multsv), 0)
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

## Fields
- `coefs::Dict{Int, Rational{Int}}`: the coefficients of the band representations, stored as
  a dictionary mapping key-value pairs to index-coefficient pairs, corresponding to indices
  into `brs`.
- `brs::Collection{NewBandRep{D}}`: the band representations referenced by `coefs`.

## Example
```julia
julia> brs = calc_bandreps(2);

julia> cbr = CompositeBandRep([1, 1, 2, 6], brs)

# example: a fragilely topological combination of band representations
julia> cbr = CompositeBandRep(Dict(16=>-1//1, 7=>1//1, 6=>1//1, 3=>1//1), brs)
16-irrep CompositeBandRep{3}:
 -(1a|Aᵤ) + (1e|Ag) + (1f|Aᵤ) + (1g|Ag) (2 bands)

julia> SymmetryVector(cbr)
 [2Z₁⁺, 2Y₁⁻, 2U₁⁻, 2X₁⁺, 2T₁⁺, 2Γ₁⁺, 2V₁⁺, 2R₁⁺] (2 bands)

# example: nontrivial topology & fractional coefficients (but integer irrep multiplicities)
julia> cbr = CompositeBandRep{3}(Dict(16=>1, 15=>1//4, 13=>-1//4, 11=>1//4, 9=>1//4,
                                      7=>1//4, 5=>-1//4, 3=>-1//4, 1=>-1//4),
                                 brs)
16-irrep CompositeBandRep{3}:
 -(1/4)×(1f|Ag) - (1/4)×(1b|Ag) + (1a|Aᵤ) + (1/4)×(1a|Ag) + (1/4)×(1c|Ag) + (1/4)×(1e|Ag) + (1/4)×(1d|Ag) - (1/4)×(1g|Ag) - (1/4)×(1h|Ag) (1 band)

julia> SymmetryVector(cbr)
16-irrep SymmetryVector{3}:
 [Z₁⁺, Y₁⁻, U₁⁻, X₁⁻, T₁⁻, Γ₁⁻, V₁⁻, R₁⁻] (1 band)
```
"""
struct CompositeBandRep{D} <: AbstractSymmetryVector{D}
    coefs :: Dict{Int, Rational{Int}}
    brs   :: Collection{NewBandRep{D}}
end

# ::: Convenience constructor :::
"""
    CompositeBandRep(idxs::Vector{Int}, brs::Collection{<:NewBandRep})

Return a `CompositeBandRepSet` whose symmetry content is equal to the sum the band
representations of `brs` over `idxs`.
In terms of irrep multiplicity, this is equivalent to `sum(brs[idxs])`, in the sense that
`SymmetryVector(CompositeBandRep(idxs, brs))` is equal to `sum(brs[idxs])`. 

The difference, and primary motivation for using `CompositeBandRep`, is that  
`CompositeBandRep` retains information about which band representations are included and
with what multiplicity (the multiplicity of the `i`th `brs`-element being equaling to
`count(==(i), idxs)`).

## Example
```julia
julia> brs = calc_bandreps(2);

julia> cbr = CompositeBandRep([1, 1, 2, 6], brs)
16-irrepCompositeBandRep{3}:
 (1f|Aᵤ) + (1h|Aᵤ) + 2(1h|Ag) (4 bands)

julia> cbr == brs[1] + brs[1] + brs[2] + brs[6]
true

julia> SymmetryVector(cbr)
16-irrep SymmetryVector{3}:
 [2Z₁⁺+2Z₁⁻, Y₁⁺+3Y₁⁻, 2U₁⁺+2U₁⁻, 2X₁⁺+2X₁⁻, 3T₁⁺+T₁⁻, 2Γ₁⁺+2Γ₁⁻, 3V₁⁺+V₁⁻, R₁⁺+3R₁⁻] (4 bands)
```
"""
function CompositeBandRep(idxs::Vector{Int}, brs::Collection{<:NewBandRep})
    coefs = Dict{Int, Rational{Int}}()
    for i in idxs
        @boundscheck checkbounds(brs, i)
        coefs[i] = get(coefs, i, zero(valtype(coefs))) + 1
    end
    CompositeBandRep(coefs, brs)
end

# ::: AbstractSymmetryVector interface :::
function SymmetryVector(cbr::CompositeBandRep{D}) where {D}
    brs = cbr.brs
    lgirsv = irreps(brs)
    multsv_r = zeros.(Rational{Int}, size.(multiplicities(first(brs))))
    μ = 0
    for (j, c) in cbr.coefs
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
    μ_r = sum(c * occupation(cbr.brs[j]) for (j, c) in cbr.coefs; 
              init=zero(valtype(cbr.coefs)))
    isinteger(μ_r) || error(lazy"CompositeBandRep has non-integer occupation (= $μ_r)")
    return round(Int, μ_r)
end

# ::: overloads to avoid materializing a SymmetryVector unnecessarily :::
irreps(cbr::CompositeBandRep) = irreps(first(cbr.brs))
klabels(cbr::CompositeBandRep) = klabels(first(cbr.brs))
irreplabels(cbr::CompositeBandRep) = irreplabels(first(cbr.brs))
num(cbr::CompositeBandRep) = num(first(cbr.brs))

# ::: AbstractArray interface :::
Base.size(cbr::CompositeBandRep) = size(first(cbr.brs))
Base.getindex(cbr::CompositeBandRep{D}, i::Int) where {D} = Int(sum(c * cbr.brs[j][i] for (j, c) in cbr.coefs; init=zero(valtype(cbr.coefs))))

# ::: Arithmetic operations :::
function Base.:+(cbr1::CompositeBandRep{D}, cbr2::CompositeBandRep{D}) where D
    if cbr1.brs !== cbr2.brs
        error("provided CompositeBandReps must reference egal `brs`")
    end
    coefs = mergewith(+, cbr1.coefs, cbr2.coefs)
    foreach(coefs) do (k, c)
        iszero(c) && delete!(coefs, k)
    end
    return CompositeBandRep{D}(coefs, cbr1.brs)
end
function Base.:-(cbr::CompositeBandRep{D}) where D
    coefs = Dict(k=>-c for (k,c) in cbr.coefs)
    CompositeBandRep{D}(coefs, cbr.brs)
end
Base.:-(cbr1::CompositeBandRep{D}, cbr2::CompositeBandRep{D}) where D = cbr1 + (-cbr2)
function Base.:*(cbr::CompositeBandRep{D}, n::Integer) where D
    coefs = Dict(k=>c*n for (k,c) in cbr.coefs)
    return CompositeBandRep{D}(coefs, cbr.brs)
end
function Base.zero(cbr::CompositeBandRep{D}) where D
    CompositeBandRep{D}(typeof(cbr.coefs)(), cbr.brs)
end