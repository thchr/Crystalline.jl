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

Parse a string `s` to a `SymmetryVector` over the irreps provided in `lgirsv`. The labels of
`lgirsv` and `s` must be provided in a shared convention.

## Example
```julia
julia> brs = calc_bandreps(220);
julia> lgirsv = irreps(brs); # irreps at P, H, Γ, & PA
julia> s = "[2P₃, 4N₁, H₁H₂+H₄H₅, Γ₁+Γ₂+Γ₄+Γ₅, 2PA₃]"
julia> parse(SymmetryVector, s, lgirsv)
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
function Base.zero(n::AbstractSymmetryVector{D}) where D
    SymmetryVector{D}(irreps(n), zero.(multiplicities(n)), 0)
end

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