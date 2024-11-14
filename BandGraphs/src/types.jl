# ---------------------------------------------------------------------------------------- #

"""
    SymmetryVector(nv, irlabs_nv, lgirsd)

Build a structured representation of a symmetry vector `nv` whose entries correspond to the
irrep labels in `irlabs_nv`, mapping these labels to the full irrep information provided in
`lgirsd`.
Note that the sorting of labels in `lgirsd` and (`nv`, `irlabs_nv`) is allowed to differ.
"""
function Crystalline.SymmetryVector(
            nv::AbstractVector{<:Integer},
            irlabs_nv::AbstractVector{String},
            lgirsd::Dict{String, <:AbstractVector{LGIrrep{D}}}) where D
    klabs = unique(klabel.(irlabs_nv))
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
Base.length(p::Partition) = length(p.lgirs)
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

struct SubGraph{D} <: AbstractSubGraph{D} # a 2-partite subgraph
    p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    A        :: Matrix{Int}  # a canonically ordered subgraph adjacency matrix
    monodromy_tie_via_Ω :: Bool # whether the subgraph is artificially tied to Ω to enforce
                                # energy-sameness between monodromy-related irreps; 
                                # indicate to avoid any permutations
end
function SubGraph(p_max::Partition, p_nonmax::Partition, A::AbstractMatrix)
    return SubGraph(p_max, p_nonmax, A, #=monodromy_tie_via_Ω=# false)
end

function Base.show(io :: IO, s :: SubGraph)
    print(io, "[")
    for i in 1:size(s.A, 1)
        print_identified_lgir_label(io, s.p_max.lgirs[i], s.p_max, i)
        print(io, " ↓ ")
        js = findall(!iszero, @view s.A[i,:])
        for (nⱼ,j) in enumerate(js)
            lgirⱼ = s.p_nonmax.lgirs[j]
            multiplicity = s.A[i,j] // Crystalline.irdim(lgirⱼ)
            isone(multiplicity) || print(io, multiplicity)
            print_identified_lgir_label(io, lgirⱼ, s.p_nonmax, j)
            nⱼ == length(js) || print(io, " + ")
        end
        i == size(s.A, 1) || print(io, ", ")
    end
    print(io, "]")
end

function print_identified_lgir_label(
            io :: IO, lgir :: LGIrrep{D}, p :: Partition{D}, idx :: Int
            ) where D
    # print a superscript index to identify subduction to nonmaximal irreps for which there
    # may be multiple copies; only print an identifier if there are multiple "copies", and
    # only if the subgraph is not `pinned`; the index is a local counter within the set of 
    # identical/"copied" nonmaximal irreps
    print(io, label(lgir))
    for m in p.multiples
        length(m) == 1 && continue # don't print a local index if no copies
        local_iridx = findfirst(==(idx), m)
        isnothing(local_iridx) && continue
        printstyled(io, Crystalline.supscriptify(string(local_iridx)); color=:light_black)
    end
end
function Base.show(io :: IO, ::MIME"text/plain", s :: SubGraph)
    summary(io, s)
    println(io, ": ", s)
    
    # max- and nonmax- labels
    max_labels = map(1:length(s.p_max)) do i
        io′ = IOBuffer()
        print_identified_lgir_label(io′, s.p_max.lgirs[i], s.p_max, i)
        String(take!(io′))
    end
    nonmax_labels = map(1:length(s.p_nonmax)) do i
        io′ = IOBuffer()
        print_identified_lgir_label(io′, s.p_nonmax.lgirs[i], s.p_nonmax, i)
        String(take!(io′))
    end
    pretty_table(io, 
        s.A;
        row_labels = max_labels,
        header = nonmax_labels,
        vlines = [1,],
        hlines = [:begin, 1, :end],
        row_label_alignment = :l,
        alignment = :c,
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        highlighters = (Highlighter((data,i,j) -> iszero(data[i,j]), crayon"dark_gray"), )
        )
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
Graphs.nv(bandg::BandGraph) = last(bandg.partitions[end].iridxs)

function Crystalline.occupation(bandg::BandGraph)
    partitions = bandg.partitions
    μ = sum(irdim, first(partitions).lgirs)
    for p in @view partitions[2:end]
        if sum(irdim, p.lgirs) ≠ μ
            error("uneqal number of bands in different partitions!")
        end
    end
    return μ
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

"""
    VectorProductIterator{T}(vs::Vector{Vector{T}})

An iterator over the the cartesian product of the vectors in `vs`. 
Pratically equivalent to `ProductIterator(vs...)` but avoids the splatting of the input and
associated type-instabilities.

Used internally to represent the combinatorial product of column-permutations of subgraphs
in `SubGraphPermutations`.
"""
struct VectorProductIterator{T} <: AbstractVector{Vector{T}}
    vs :: Vector{Vector{T}}
end
Base.size(p::VectorProductIterator) = (prod(length, p.vs; init=1),)

function init_product_iterator_state(p::VectorProductIterator)
    state = ones(Int, length(p.vs))
    state[1] = 0
    return state # `state[j]` indexes into p.vs[j]
end

function Base.iterate(p::VectorProductIterator{T}) where T
    any(isempty, p.vs) && return nothing # check for empty vectors; no iterants
    return iterate(p, init_product_iterator_state(p))
end

function Base.iterate(p::VectorProductIterator{T}, state::Vector{Int}) where T
    finished = increment_state!(p, state, 1)
    finished && return nothing
    v = [vⱼ[i] for (vⱼ, i) in zip(p.vs, state)]
    if T isa AbstractVector
        # Hack: for use in `SubGraphPermutations`; see TODO there.
        return reduce(vcat, v), state
    else
        return v, state
    end
end

function increment_state!(p::VectorProductIterator, state::Vector{Int}, j::Int)
    ip1 = state[j]+1
    if ip1 ≤ length(p.vs[j])
        state[j] = ip1
        return false # not finished
    else
        j == length(p.vs) && return true # finished
        state[j] = 1
        return increment_state!(p, state, j+1)
    end
end

function Base.getindex(p::VectorProductIterator{T}, idx::Integer) where T
    @boundscheck 1 ≤ idx ≤ length(p) || throw(BoundsError(p, idx))
    idx = idx - 1
    v = Vector{T}(undef, length(p.vs))
    @inbounds for (j, vⱼ) in enumerate(p.vs)
        r1 = length(vⱼ) # `= counts[i]`
        idx′ = div(idx, r1)
        sub_idx_j = idx - r1*idx′ + 1 # `j`th subscript in the linear-to-subscript problem
        idx = idx′
        @inbounds v[j] = vⱼ[sub_idx_j]
    end
    return reduce(vcat, v) # TODO: change back to just `v`; temporary change
end

# ---------------------------------------------------------------------------------------- #
@enum ColsOrRowsEnum::Int8 begin
    UNPERMUTED = 0 # neither row nor columns permuted
    ROWS = 1    
    COLS = 2
end

mutable struct SubGraphPermutations{D} <: AbstractSubGraph{D}
    const subgraph :: SubGraph{D}
    permutations   :: Union{Nothing, VectorProductIterator{Vector{Int}}}
    cols_or_rows   :: ColsOrRowsEnum # permutate over columns or rows of `subgraph.A`
    pinned         :: Bool

function Base.length(subgraph_ps :: SubGraphPermutations)
    permutations = subgraph_ps.permutations
    return isnothing(permutations) ? 1 : length(permutations)
end

function Base.iterate(subgraph_ps :: SubGraphPermutations)
    subgraph = subgraph_ps.subgraph
    A = subgraph.A
    if isnothing(subgraph_ps.permutations)
        return subgraph, Vector{Int}[]
    else
        permutation, state = iterate(something(subgraph_ps.permutations))
        A′ = if subgraph_ps.cols_or_rows == COLS
            A[:,permutation]
        elseif subgraph_ps.cols_or_rows == ROWS
            A[permutation,:]
        else # UNPERMUTED
            error("unexpected error: non-nothing permutation for unpermuted subgraph")
        end
        subgraph′ = SubGraph{D}(subgraph.p_max, subgraph.p_nonmax, A′, subgraph.monodromy_tie_via_Ω)
        return subgraph′, state
    end
end

function Base.iterate(subgraph_ps :: SubGraphPermutations, state::Vector{Int})
    subgraph = subgraph_ps.subgraph
    A = subgraph.A
    if isnothing(subgraph_ps.permutations)
        return nothing
    else
        permutation, state = iterate(something(subgraph_ps.permutations), state)
        A′ = if subgraph_ps.cols_or_rows == COLS
            A[:,permutation]
        elseif subgraph_ps.cols_or_rows == ROWS
            A[permutation,:]
        else # UNPERMUTED
            error("unexpected error: non-nothing permutation for unpermuted subgraph")
        end
        subgraph′ = SubGraph{D}(subgraph.p_max, subgraph.p_nonmax, A′, subgraph.monodromy_tie_via_Ω)
        return subgraph′, state
    end
end

function Base.getindex(subgraph_ps::SubGraphPermutations{D}, idx::Integer) where D
    @boundscheck 1 ≤ idx ≤ length(subgraph_ps) || throw(BoundsError(subgraph_ps, idx))
    subgraph = subgraph_ps.subgraph
    permutations = subgraph_ps.permutations
    A = if isnothing(permutations)
        subgraph.A
    else
        permutation = something(permutations)[idx]
        if subgraph_ps.cols_or_rows == COLS
            subgraph.A[:,permutation]
        elseif subgraph_ps.cols_or_rows == ROWS
            subgraph.A[permutation,:]
        else # UNPERMUTED
            error("unexpected error: non-nothing permutation for unpermuted subgraph")
        end
    end    
    return SubGraph{D}(subgraph.p_max, subgraph.p_nonmax, A, subgraph.monodromy_tie_via_Ω)
end
# ---------------------------------------------------------------------------------------- #

struct BandGraphPermutations{D} <: AbstractVector{BandGraph{D}}
    partitions   :: Vector{Partition{D}}
    subgraphs_ps :: Vector{SubGraphPermutations{D}}
end

# NB: Since the number of permutations potentially can be very, very large - larger than
# typemax(Int), `length` ought in principle to work with `BigInt`s. But this is very slow
# and since `length` is used ubiquitously, this in turns slows down iteration of 
# `BandGraphPermutations`. So, here, we use ordinary `Int` & throw if there is an overflow.
# To compute a definite, non-throwing value, use `safe_length(bandgp)`, which uses `BigInt`.
function Base.length(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    n = 1
    for subgraph_ps in subgraphs_ps
        n = Base.checked_mul(n, length(subgraph_ps)) # check for & throw on overflow
    end
    return n
end
Base.size(bandgp :: BandGraphPermutations) = (length(bandgp),)

function safe_length(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    n = one(BigInt)
    for subgraph_ps in subgraphs_ps
        n = Base.checked_mul(n, length(subgraph_ps)) # check for & throw on overflow
    end
    return n
end

function permutation_counts(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    return [length(subgraph_ps) for subgraph_ps in subgraphs_ps]
end

function permutation_info(bandgp :: BandGraphPermutations)
    subgraphs_ps = bandgp.subgraphs_ps
    n = BigInt(1)
    foreach(subgraphs_ps) do subgraph_ps
        subgraph = subgraph_ps.subgraph
        Np = length(subgraph_ps) # aggregate number of distinct same-irrep permutations, across irreps
        if subgraph_ps.cols_or_rows ≠ UNPERMUTED
            n *= Np
        end
        printstyled(stdout,
            subgraph.p_max.klab, " → ", subgraph.p_nonmax.klab, ": ",
            Np,
            subgraph_ps.pinned ? "\t(pinned)" : "";
            color=subgraph_ps.pinned ? :light_black : :normal)
        println(stdout)
    end
    println(stdout,     "Total permutations:          ", n)
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
    # problematic, however, since they are either type-unstable or allocate a lot.
    # Instead, we copy the logic from `Base.getindex` on `CartesianIndices`, specifically
    # from `_ind2sub_recurse`, converting also from a recursive implementation to a loop.
    @boundscheck 1 ≤ idx ≤ length(bandgp) || throw(BoundsError(bandgp, idx))
    idx = idx-1
    for (i, subgraph_ps) in enumerate(subgraphs_ps)
        r1 = length(subgraph_ps) # `= counts[i]`
        idx′ = div(idx, r1)
        sub_idx_i = idx - r1*idx′ + 1 # `i`th subscript in the linear-to-subscript problem
        idx = idx′

        subgraphs[i] = subgraph_ps[sub_idx_i]
    end
    return BandGraph{D}(subgraphs, partitions)
end