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
    for (nᵢ, irlabᵢ) in zip(nv, irlabs_nv)
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

function Base.show(io :: IO, ::MIME"text/plain", n :: SymVector)
    summary(io, n)
    print(io, " over ", length(n.mults), " k-points:\n ")
    show(io, n)
end

function Base.show(io :: IO, n :: SymVector)
    print(io, "[")
    for (i, (mults_k, lgirs_k)) in enumerate(zip(n.mults, n.lgirsv))
        printstyled(io, 
                    Crystalline.symvec2string(mults_k, label.(lgirs_k); braces=false);
                    color=iseven(i) ? :light_blue : :normal)
        i ≠ length(n.mults) && print(io, ", ")
    end
    print(io, "]")
    printstyled(io, " (", n.μ, " band", n.μ ≠ 1 ? "s" : "", ")"; color=:light_black)
end

Base.collect(n :: SymVector) = reduce(vcat, n.mults) # flatten `n` to a simple vector `nv`

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

function Base.show(io :: IO, p :: Partition)
    Nir = length(p.lgirs)
    print(io, Nir, " irrep", Nir ≠ 1 ? "s" : "", " (", join(label.(p.lgirs), ", "), ")",
              " at ", p.klab, " = ", position(group(first(p.lgirs))))
    printstyled(io, " (", p.maximal ? "" : "non", "maximal)"; color=:light_black)
end


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

function Base.show(io :: IO, s :: SubGraph)
    print(io, "[")
    for i in 1:size(s.A, 1)
        print(io, label(s.p_max.lgirs[i]), " ↓ ")
        js = findall(!iszero, @view s.A[i,:])
        for (nⱼ,j) in enumerate(js)
            lgirⱼ = s.p_nonmax.lgirs[j]
            multiplicity = s.A[i,j] // Crystalline.irdim(lgirⱼ)
            isone(multiplicity) || print(io, multiplicity)
            print(io, label(lgirⱼ))
            if !s.pinned
                # print a superscript index to identify subduction to nonmaximal irreps for 
                # which there may be multiple copies; only print an identifier if there
                # are multiple "copies", and only if the subgraph is not pinned; the index
                # is a local counter within the set of identical/"copied" nonmaximal irreps
                for m in s.p_nonmax.multiples
                    length(m) == 1 && continue # don't print a local index if no copies
                    local_iridx = findfirst(==(j), m)
                    isnothing(local_iridx) && continue
                    printstyled(io, Crystalline.supscriptify(string(local_iridx)); 
                                    color=:light_black)
                end
            end
            nⱼ == length(js) || print(io, " + ")
        end
        i == size(s.A, 1) || print(io, ", ")
    end
    print(io, "]")
    s.pinned && printstyled(io, " (pinned)"; color=:light_black)
end
# ---------------------------------------------------------------------------------------- #

struct BandGraph{D}
    subgraphs  :: Vector{SubGraph{D}}
    partitions :: Vector{Partition{D}}
end
function BandGraph(
            subgraphs::AbstractVector{Partition{D}},
            partitions::AbstractVector{Partition{D}}
            ) where D
    BandGraph{D}(Vector(subgraphs), Vector(partitions))
end

function Base.show(io :: IO, bandg ::BandGraph)
    print(io, "{")
    Ns = length(bandg.subgraphs)
    for (i, s) in enumerate(bandg.subgraphs)
        print(io, s)
        i == Ns || print(io, ", ")
    end
    print(io, "}")
end
function Base.show(io :: IO, ::MIME"text/plain", bandg ::BandGraph)
    Base.summary(io, bandg)
    Np = length(bandg.partitions)
    Np_max = count(p->p.maximal, bandg.partitions)
    Ns = length(bandg.subgraphs)
    println(io, " with ", Ns, " subgraph", Ns ≠ 1 ? "s" : "",
                " over ", Np, " partition", Np ≠ 1 ? "s" : "",
                " (", Np_max, " maximal):")
    for (i, s) in enumerate(bandg.subgraphs)
        print(io, " ", s)
        i == Ns || println(io)
    end
end

# ---------------------------------------------------------------------------------------- #

mutable struct SubGraphPermutations{D} <: AbstractSubGraph{D}
    const p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    const p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    const As       :: Vector{Matrix{Int}} # valid/relevant permutations of the subgraph adj mat
    pinned         :: Bool
end
Base.length(subgraph_ps::SubGraphPermutations) = length(subgraph_ps.As)

# ---------------------------------------------------------------------------------------- #

struct BandGraphPermutations{D} <: AbstractVector{BandGraph{D}}
    partitions   :: Vector{Partition{D}}
    subgraphs_ps :: Vector{SubGraphPermutations{D}}
end

function Base.length(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    n = 1
    for subgraph_ps in subgraphs_ps
        n *= length(subgraph_ps)
    end
    return n
end
Base.size(bandgp :: BandGraphPermutations) = (length(bandgp),)

function permutation_counts(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    return [length(subgraph_ps) for subgraph_ps in subgraphs_ps]
end

function permutation_info(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    n = n_no_pins = 1
    foreach(subgraphs_ps) do subgraph_ps       
        multiples = subgraph_ps.p_nonmax.multiples 
        pss = collect.(permutations.(multiples))
        Np = prod(length, pss) # aggregate number of same-irrep permutations, across irreps
        if !subgraph_ps.pinned
            @assert Np == length(subgraph_ps)
            n *= Np
        end
        n_no_pins *= Np
        printstyled(stdout,
            subgraph_ps.p_max.klab, " → ", subgraph_ps.p_nonmax.klab, ": ",
            Np,
            subgraph_ps.pinned ? "\t(pinned)" : "";
            color=subgraph_ps.pinned ? :light_black : :normal)
        println(stdout)
    end
    println(stdout,     "Total permutations:          ", n)
    printstyled(stdout, "Total permutations w/o pins: ", n_no_pins, color=:light_black)
    println(stdout)
end

function Base.getindex(
            bandgp :: BandGraphPermutations{D},
            idx :: Integer
            ) where D
    subgraphs_ps = bandgp.subgraphs_ps
    partitions = bandgp.partitions
    subgraphs = Vector{SubGraph{D}}(undef, length(subgraphs_ps))

    # Conceptually, we now want to find the `idx`th permutation from the following set of
    # ranges:
    #       counts = permutation_counts(bandgp)
    #       (1:counts[1], 1:counts[2], ..., 1:counts[end])
    # This is equivalent to the following linear-index-to-subscript-index problem:
    #       CartesianIndices(Tuple(counts))[idx]
    # or, equivalently, to `collect(Base.product(counts...))[idx]`; both these variants are
    # problematic, however, since they are either type-unstable or allocates a lot.
    # Instead, we copy the logic from `Base.getindex` on `CartesianIndices`, specifically
    # from `_ind2sub_recurse`, converting also from a recursive implementation to a loop.
    @boundscheck 1 ≤ idx ≤ length(bandgp) || throw(BoundsError(bandgp, idx))
    idx = idx-1
    for (i, subgraph_ps) in enumerate(subgraphs_ps)
        r1 = length(subgraph_ps) # `= counts[i]`
        idx′ = div(idx, r1)
        sub_idx_i = idx - r1*idx′ + 1 # `i`th subscript in the linear-to-subscript problem

        subgraphs[i] = SubGraph{D}(subgraph_ps.p_max,
                                   subgraph_ps.p_nonmax,
                                   subgraph_ps.As[sub_idx_i],
                                   subgraph_ps.pinned)
        idx = idx′
    end
    return BandGraph{D}(subgraphs, partitions)
end
