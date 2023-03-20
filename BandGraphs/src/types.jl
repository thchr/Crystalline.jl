# ---------------------------------------------------------------------------------------- #

struct SymVector{D}
    mults  :: Vector{Vector{Int}}
    lgirsv :: Vector{Vector{LGIrrep{D}}}
    klabs  :: Vector{String}
    μ      :: Int # connectivity
end
Crystalline.irreplabels(n::SymVector) = [label.(lgirs) for lgirs in n.lgirsv]

"""
    SymVector(nv, irlabs_nv, lgirsd)

Build a structured representation of a symmetry vector `nv` whose entries correspond to the
irrep labels in `irlabs_nv`, mapping these labels to the full irrep information provided in
`lgirsd`.
"""
function SymVector(
            nv::AbstractVector{<:Integer},
            irlabs_nv::AbstractVector{String},
            lgirsd::Dict{String, Vector{LGIrrep{D}}}) where D
    klabs = unique(klabel.(irlabs_nv))
    Nk = length(klabs)
    mults = [Int[] for _ in 1:Nk]
    lgirsv = [LGIrrep{D}[] for _ in 1:Nk]
    j = 1
    for (i, (nᵢ, irlabᵢ)) in enumerate(zip(nv, irlabs_nv))
        klabᵢ = klabel(irlabᵢ)
        if klabᵢ != klabs[j]
            j += 1
        end
        push!(mults[j], nᵢ)

        # find associated irrep in `lgirsd[klabᵢ]`
        lgirsⱼ = lgirsd[klabᵢ]
        iridxᵢ = findfirst(lgir -> label(lgir) == irlabᵢ, lgirsⱼ)
        if isnothing(iridxᵢ)
            _throw_failed_to_find_irrep(irlabᵢ, lgirsⱼ)
        else
            push!(lgirsv[j], lgirsⱼ[something(iridxᵢ)])
        end
    end
    if length(nv) ≠ length(irlabs_nv)+1
        error("n must contain its band connectivity")
    end
    μ = nv[end]
    return SymVector(mults, lgirsv, string.(klabs), μ)
end
function _throw_failed_to_find_irrep(irlabᵢ, lgirsⱼ)
    error("failed to find an irrep label \"$irlabᵢ\" in `lgirsd`; available irrep labels \
           at the considered k-point were $(label.(lgirsⱼ))")
end

# ---------------------------------------------------------------------------------------- #

struct Partition{D}
    klab      :: String
    lgirs     :: Vector{LGIrrep{D}}
    multiples :: Vector{UnitRange{Int}} # nodes to potentially shuffle
    maximal   :: Bool
    kidx      :: Int # index in the block form (i.e., index of k-manifold)
    iridxs    :: UnitRange{Int} # _global_ irrep indices
end
Base.length(p::Partition) = length(irreplabels(p))
Crystalline.irreplabels(p::Partition) = label.(p.lgirs)
irdims(p::Partition) = Crystalline.irdim.(p.lgirs)

# ---------------------------------------------------------------------------------------- #

abstract type AbstractSubGraph{D} end

mutable struct SubGraph{D} <: AbstractSubGraph{D} # a 2-partite subgraph
    const p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    const p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    const A        :: Matrix{Int}  # a canonically ordered subgraph adjacency matrix
    pinned         :: Bool
end
function SubGraph(p_max::Partition, p_nonmax::Partition, A::AbstractMatrix)
    return SubGraph(p_max, p_nonmax, A, false)
end

# ---------------------------------------------------------------------------------------- #

mutable struct SubGraphPermutations{D} <: AbstractSubGraph{D}
    const p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    const p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    const As       :: Vector{Matrix{Int}} # valid/relevant permutations of the subgraph adj mat
    pinned         :: Bool
end
Base.length(subgraphp::SubGraphPermutations) = length(subgraphp.As)

# ---------------------------------------------------------------------------------------- #