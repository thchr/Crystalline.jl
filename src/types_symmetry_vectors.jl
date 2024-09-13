# ---------------------------------------------------------------------------------------- #
# AbstractSymmetryVector
abstract type AbstractSymmetryVector{D} <: AbstractVector{Int} end

# ::: API :::
"""
    irreps(n::AbstractSymmetryVector{D}) -> AbstractVector{<:Collection{<:AbstractIrrep{D}}}

Return the irreps referenced by `n`. 

The returned value is an `AbstractVector` of `Collection{<:AbstractIrrep}`s, with irreps for
distinct groups, usually associated with specific **k**-manifolds, belonging to the same
`Collection`.

See also [`multiplicities(::AbstractSymmetryVector)`](@ref).
"""
function irreps(::AbstractSymmetryVector) end

"""
    multiplicities(n::AbstractSymmetryVector) -> AbstractVector{<:AbstractVector{Int}}

Return the multiplicities of the irreps referenced by `n`.

See also [`irreps(::AbstractSymmetryVector)`](@ref).
"""
function multiplicities(::AbstractSymmetryVector) end

"""
    occupation(n::AbstractSymmetryVector) -> Int

Return the occupation of (i.e., number of bands contained within) `n`.
"""
function occupation(::AbstractSymmetryVector) end

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
    for mults in n.multsv
        N = length(mults)
        copyto!(nv, i, mults, 1, N)
        i += N
    end
    nv[end] = n.occupation
    return nv
end
irreplabels(n::SymmetryVector) = [label(ir) for ir in Iterators.flatten(n.lgirsv)]
klabels(n::SymmetryVector) = [klabel(first(irs)) for irs in n.lgirsv]

# ---------------------------------------------------------------------------------------- #
# NewBandRep

struct NewBandRep{D} <: AbstractSymmetryVector{D}
    siteir       :: SiteIrrep{D}
    n            :: SymmetryVector{D}
    timereversal :: Bool
    spinful      :: Bool
end

# ::: AbstractSymmetryVector interface :::
irreps(br::NewBandRep) = irreps(br.n)
multiplicities(br::NewBandRep) = multiplicities(br.n)
occupation(br::NewBandRep) = occupation(br.n)

# ::: AbstractArray interface beyond AbstractSymmetryVector :::
Base.setindex!(br::NewBandRep{D}, v::Int, i::Int) where D = (br.n[i] = v)
function Base.similar(br::NewBandRep{D}) where D
    NewBandRep{D}(br.siteir, similar(br.n), br.timereversal, br.spinful)
end
Base.Vector(br::NewBandRep) = Vector(br.n)

# ::: Utilities :::
num(br::NewBandRep) = num(br.siteir)
irreplabels(br::NewBandRep) = irreplabels(br.n)
klabels(br::NewBandRep) = klabels(br.n)

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
irreplabels(brs::Collection{<:NewBandRep}) = irreplabels(first(brs).n)
klabels(brs::Collection{<:NewBandRep}) = klabels(first(brs).n)

# ::: Conversion to BandRepSet :::
function Base.convert(::Type{BandRepSet}, brs::Collection{<:NewBandRep})
    sgnum = num(first(brs))
    bandreps = convert.(Ref(BandRep), brs)
    kvs = [position(first(lgirs)) for lgirs in irreps(first(brs))]
    klabs = klabels(first(brs))
    irlabs = irreplabels(first(brs))
    spinful = first(brs).spinful
    timereversal = first(brs).timereversal
    return BandRepSet(sgnum, bandreps, kvs, klabs, irlabs, spinful, timereversal)
end