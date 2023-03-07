module BandGraphs

using Crystalline
using JLD2: jldopen
using PrettyTables
using LinearAlgebra
using Graphs
using MetaGraphsNext
using Combinatorics: permutations # for set_pinned_subgraphs!

# ---------------------------------------------------------------------------------------- #
# EXPORTS

export
    SymVector,
    Partition,
    # TODO


# ---------------------------------------------------------------------------------------- #

# TODO: much copying from `crawl_and_write_bandpaths.jl` below: fix
include("subduction-types.jl")

subductionsd = jldopen("data/connections/3d/subductions-tr.jld2", "r")["subductionsd"]

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
using BlockArrays

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
mutable struct SubGraph{D} # a 2-partite subgraph
    const p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    const p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    const A        :: Matrix{Int} # a canonically ordered subgraph adjacency matrix
    pinned         :: Bool
end
SubGraph(p_max::Partition, p_nonmax::Partition, A::AbstractMatrix) = SubGraph(p_max, p_nonmax, A, false)

mutable struct SubGraphPermutations{D}
    const p_max    :: Partition{D} # maximal k-manifold (high-symmetry point)
    const p_nonmax :: Partition{D} # non-maximal k-manifold (high-symmetry line/plane)
    const As       :: Vector{Matrix{Int}} # valid/relevant permutations of the subgraph adj mat
    pinned         :: Bool
end
Base.length(subgraphp::SubGraphPermutations) = length(subgraphp.As)
function build_partition(klab, nsₖ, lgirsₖ::Vector{LGIrrep{D}}, kidx, seen_irs) where D
    lgirs = LGIrrep{D}[]
    multiples = UnitRange{Int}[]
    max_local_iridx = 0
    for (nₖᵢ, lgirₖᵢ) in zip(nsₖ, lgirsₖ)
        nₖᵢ == 0 && continue # irrep not included, nothing to add
        append!(lgirs, Iterators.repeated(lgirₖᵢ, nₖᵢ))
        push!(multiples, (max_local_iridx+1):length(lgirs))
        max_local_iridx = length(lgirs)
    end
    global_iridxs = seen_irs .+ (1:max_local_iridx) # indices of irrep nodes in "global sorting"
    return Partition(klab, lgirs, multiples, true, kidx, global_iridxs)
end

function build_subgraphs(n::SymVector{D}, subductions,
                         lgirsd::Dict{String, Vector{LGIrrep{D}}}) where D

    seen_irs = 0 # current number of irreps aggregated across all partitions
    kidx = 0 # number of distinct k-points we've seen

    # maximal k-point partitions
    partitions_max = Partition{D}[]
    for (klab, nsₖ, lgirsₖ) in zip(n.klabs, n.mults, n.lgirsv)
        kidx += 1
        p_max = build_partition(klab, nsₖ, lgirsₖ, kidx, seen_irs)
        push!(partitions_max, p_max)
        seen_irs = last(p_max.iridxs)
    end

    # monodromy additions (same irrep population as their "non-monodromy"/"original" irrep
    # variants, but connect differently to non-maximal k-manifolds); PRE 96, 023310 (2017)
    for s in subductions
        s.monodromy || continue
        klab_max′ = string(s.c.kᴳ.label) # monodromy k-label (with a '′')
        klab_max = rstrip(klab_max′, '′') # "original" k-label (no '′')

        # skip if we already added info to `partitions_max`
        any(p->p.klab==klab_max′, partitions_max) && continue
        
        # now we have to find the entries of k′ ("monodromy" point) in k ("original" point)
        idx = findfirst(==(klab_max), n.klabs)
        if isnothing(idx)
            # user might have intentionally left out parts of a symmetry vector; then this
            # could fail to find the associated "original" irrep in which case it makes
            # best sense to leave it out; still, warn just in case since this should be very
            # rare
            @warn "failed to add monodromy irrep data at $klab_max′ since $klab_max is not included in symmetry data"
            continue
        else
            kidx += 1
        end
        @assert replace.(s.irlabsᴳ, Ref('′'=>"")) == label.(n.lgirsv[idx]) # verify assumption used below
        
        lgirs_max′ = map(enumerate(n.lgirsv[idx])) do (i,lgir)
            LGIrrep{D}(s.irlabsᴳ[i], lgir.g, lgir.matrices, lgir.translations, lgir.reality,
                       lgir.iscorep)
        end
        p_max = build_partition(klab_max′, n.mults[idx], lgirs_max′, kidx, seen_irs)
        push!(partitions_max, p_max)
        seen_irs = last(p_max.iridxs)
    end

    # non-maximal k-point partitions
    partitions_nonmax = Partition{D}[]
    for p_max in partitions_max # row index (maximal k)
        klab_max = p_max.klab
        idxs = findall(s -> string(s.c.kᴳ.label) == klab_max, subductions)
        subductions_from_max = @view subductions[idxs]
        for sᴳᴴ in subductions_from_max # subductions from maximal to non-maximal manifolds
            klab_nonmax = string(sᴳᴴ.c.kᴴ.label)
            if klab_nonmax ∈ getfield.(partitions_nonmax, Ref(:klab))
                continue # only need to add partition once
            else
                kidx += 1
            end

            # first, find and aggregate all the non-maximal irrep-labels as nodes, without
            # taking any care about the order of them in `irlabs`
            lgirs = LGIrrep{D}[]
            for lgir_max in p_max.lgirs
                irlab_max = label(lgir_max)
                idx_max = something(findfirst(==(irlab_max), sᴳᴴ.irlabsᴳ)) # in subduction irrep sorting
                for (nᴴ, irlabᴴ) in zip(sᴳᴴ.table[idx_max,:], sᴳᴴ.irlabsᴴ)
                    if nᴴ == 0
                        continue
                    else
                        lgirs_nonmax = lgirsd[klab_nonmax]
                        iridxᴴ = findfirst(lgir->label(lgir)==irlabᴴ, lgirs_nonmax)
                        if isnothing(iridxᴴ)
                            _throw_failed_to_find_irrep(irlabᴴ, lgirs_nonmax)
                        else
                            append!(lgirs,
                                    Iterators.repeated(lgirs_nonmax[something(iridxᴴ)], nᴴ))
                        end
                    end
                end
            end
            # now, the elements in `lgirs` are not sorted: restore sorting and compute
            # `multiples` and `global_iridxs` as well
            sort!(lgirs, by = label)
            local_iridx = 1
            multiples = UnitRange{Int}[]
            while local_iridx ≤ length(lgirs)
                local_irlab = label(lgirs[local_iridx])
                nₖᵢ = count(lgir -> label(lgir)==local_irlab, lgirs)
                push!(multiples, local_iridx:local_iridx+nₖᵢ-1)
                local_iridx += nₖᵢ
            end
            global_iridxs = seen_irs .+ (1:length(lgirs)) # indices of irrep nodes in "global sorting"
            seen_irs = last(global_iridxs)

            # now we have all the node info about the non-maximal irreps at `klab_nonmax`
            p_nonmax = Partition(klab_nonmax, lgirs, multiples, false, kidx, global_iridxs)
            push!(partitions_nonmax, p_nonmax)
        end
    end
    
    # TODO: sort `partitions_max` and `partitions_nonmax` so that we benefit from maximally
    #       from the pinning order when we consider permutations: basically, want to sort
    #       in such a way that subgraphs with most permutations get pinned
    
    # build all the subgraphs that make up the full graph
    subgraphs = SubGraph{D}[]
    for p_max in partitions_max # maximal k-manifolds
        i = p_max.kidx
        klab_max  = p_max.klab
        lgirs_max = p_max.lgirs
        for p_nonmax in partitions_nonmax # non-maximal k-manifolds
            j = p_nonmax.kidx
            klab_nonmax  = p_nonmax.klab
            lgirs_nonmax = p_nonmax.lgirs

            idx = findfirst(subductions) do s
                string(s.c.kᴳ.label) == klab_max && string(s.c.kᴴ.label) == klab_nonmax
            end
            if !isnothing(idx)
                Aᵢⱼ = zeros(Int, length(lgirs_max), length(lgirs_nonmax))
                # pick the "most diagonal" connection as our canonical, unpermuted matrix
                s = subductions[idx]
                for (q, lgir_max) in enumerate(lgirs_max)
                    q′ = something(findfirst(==(label(lgir_max)), s.irlabsᴳ))
                    for (r′,v) in enumerate(@view s.table[q′,:])
                        irlab_nonmax = s.irlabsᴴ[r′]
                        r = 0
                        for _ in 1:v
                            while true
                                r = something(findnext(lgir->label(lgir)==irlab_nonmax,
                                              lgirs_nonmax, r+1))
                                if iszero(Aᵢⱼ[:,r])
                                    Aᵢⱼ[q,r] = Crystalline.irdim(lgirs_nonmax[r])
                                    break
                                end
                            end
                        end
                    end
                end
                push!(subgraphs, SubGraph(p_max, p_nonmax, Aᵢⱼ))
            end
        end
    end

    return subgraphs, partitions_max, partitions_nonmax
end

# build adjacency matrix as a block matrix from subgraphs
function assemble_adjacency(subgraphs, partitions_max, partitions_nonmax)
    Nirs_total  = last(partitions_nonmax[end].iridxs)
    Nirs_max    = [length(p.lgirs) for p in partitions_max]
    Nirs_nonmax = [length(p.lgirs) for p in partitions_nonmax]
    Nirs_all    = vcat(Nirs_max, Nirs_nonmax)

    A = BlockArray{Int}(zeros(Int, Nirs_total, Nirs_total), Nirs_all, Nirs_all)
    for subgraph in subgraphs
        i, j = subgraph.p_max.kidx, subgraph.p_nonmax.kidx
        # fill in block (i,j) of the adjacency matrix
        A[Block(i), Block(j)] .= subgraph.A
        A[Block(j), Block(i)] .= subgraph.A' # symmetric matrix (undirected graph)
    end

    return A
end

# While the graph associated with `A` in principle could be constructed from e.g.
function assemble_graph(subgraphs, partitions_max, partitions_nonmax)
    g = MetaGraph(Graph();
            Label = Tuple{String, Int},
            VertexData = @NamedTuple{lgir::LGIrrep, maximal::Bool},
            EdgeData   = @NamedTuple{weight::Int})
    # add vertices to graph
    for p in Iterators.flatten((partitions_max, partitions_nonmax))
        j = 0
        for codes′ in p.multiples
            for n in eachindex(codes′)
                j += 1
                lgir = p.lgirs[j]
                # to distinguish multiple identical irrep labels, we include a incrementing
                # number (1, 2, …) along with the irrep label as the overall vertex label,
                # see also `partitions_ids(::Partition)`
                vertex_label = (label(lgir), n)
                add_vertex!(g, vertex_label, (; lgir=lgir, maximal=p.maximal))
            end
        end
    end
    @assert Graphs.nv(g) == last(partitions_nonmax[end].iridxs) # check number of vertices

    # add edges to graph
    for subgraph in subgraphs
        max_iridxs, nonmax_iridxs = subgraph.p_max.iridxs, subgraph.p_nonmax.iridxs
        for (j,a) in enumerate(eachcol(subgraph.A)) # j: local column index in block/subgraph
            j′ = nonmax_iridxs[j] # global column index in overall graph
            i = something(findfirst(≠(0), a)) # local row index in block/subgraph
            i′ = max_iridxs[i] # global column index in overall graph
            weight = a[i]
            add_edge!(g, label_for(g, i′), label_for(g, j′), (; weight=weight))
        end
    end
    return g
end

function assemble_degree(subgraphs, partitions_max, partitions_nonmax)
    # Diagonal matrix D, whose entries Dᵢᵢ are equal to the degree of the ith node (i.e.,
    # the number of edges incident upon it). We define the degree of node nᵢ as on p. 9 of
    # Vergniory's 2017 PRB, i.e. as `Pₖ * irdim(lgirₖᵢ)` with `Pₖ` denoting the number of
    # distinct k-manifolds connected to the k-point of the partition that nᵢ belongs to and
    # `lgirₖᵢ` denoting the irrep associated with nᵢ.
    # In practice, this is equal to the number of edges incident upon each node, if each
    # edge is counted with the dimensionality of the connecting non-maximal irrep.
    # These choices are e.g., motivated by the desire to have the rows of the Laplacian
    # matrix sum to zero.

    v = vcat(mapreduce(irdims, vcat, partitions_max),    # irrep dimension for each node
             mapreduce(irdims, vcat, partitions_nonmax))
    manifold_multiplicity_max = mapreduce(vcat, partitions_max) do p
        c = count(subgraphs) do subgraph
            subgraph.p_max == p
        end
        fill(c, length(p))
    end
    manifold_multiplicity_nonmax = mapreduce(vcat, partitions_nonmax) do p
        c = count(subgraphs) do subgraph
            subgraph.p_nonmax == p
        end
        fill(c, length(p))
    end
    manifold_multiplicity = vcat(manifold_multiplicity_max, manifold_multiplicity_nonmax)
    return Diagonal(v .* manifold_multiplicity)
end

function assemble_laplacian(subgraphs, partitions_max, partitions_nonmax)
    A = assemble_adjacency(subgraphs, partitions_max, partitions_nonmax)
    D = assemble_degree(subgraphs, partitions_max, partitions_nonmax)
    return D - A
end

# ---------------------------------------------------------------------------------------- #


function partition_graph(subgraphs, partitions)
    g = MetaGraph(Graph(); Label=String,
                           VertexData = @NamedTuple{maximal::Bool})
    for (i,p) in enumerate(partitions)
        @assert i == p.kidx
        add_vertex!(g, p.klab, (;maximal=p.maximal))
    end
    for subgraph in subgraphs
        src_code,  dst_code  = subgraph.p_max.kidx, subgraph.p_nonmax.kidx
        src_label, dst_label = subgraph.p_max.klab, subgraph.p_nonmax.klab
        @assert code_for(g, src_label) == src_code && code_for(g, dst_label) == dst_code
        add_edge!(g, src_label, dst_label, nothing)
    end

    return g
end

# modify the graph `kg` such that any nonmaximal node is only ever entered and exited by one
# edge; i.e., ensure the degree of nonmaximal nodes is 2. When the degree `n` is higher, we
# split the node into ∑ⱼ₌₁ⁿ⁻¹ j = (n-1)n/2 duplicated nodes.
function split_nonmaximal_nodes(kg::MetaGraph)
    kg′ = MetaGraph(Graph();
                  Label=Tuple{String, Int},
                  VertexData=@NamedTuple{klab::String, code::Int, maximal::Bool})
    for i in vertices(kg)
        klab = label_for(kg, i)
        kg[klab].maximal || break # assume that all non-maximal k-points come in sequence
        add_vertex!(kg′, (klab, 1), (; klab=klab, code=i, maximal=true))
    end
    for i in vertices(kg)
        klab_nonmax = label_for(kg, i)
        kg[klab_nonmax].maximal && continue
        neighbor_codes = neighbors(kg, i)
        neighbor_labels = label_for.(Ref(kg), neighbor_codes)
        j = 0
        for (n, klabⁿ) in enumerate(neighbor_labels)
            for m in n+1:length(neighbor_labels)
                j += 1
                klabᵐ = neighbor_labels[m]
                nonmax_vertex_label = (klab_nonmax, j)
                add_vertex!(kg′, nonmax_vertex_label,
                                 (; klab=klab_nonmax, code=i, maximal=false))
                add_edge!(kg′, (klabⁿ, 1), nonmax_vertex_label, nothing)
                add_edge!(kg′, (klabᵐ, 1), nonmax_vertex_label, nothing)
            end
        end
    end
    return kg′
end

# ---------------------------------------------------------------------------------------- #

function set_pinned_subgraphs!(subgraphs)
    # we assume and exploit a sorting of subgraphs here, where entries are sorted
    # lexicographically in the sense that
    #       [(g.p_max.kidx, g.p_nonmax.kidx) for g in subgraphs]
    # increases first in its second index (nonmax manifold), then its first (max manifold).

    seen_kidxs_nonmax = Set{Int}()
    for subgraph in subgraphs
        kidx_nonmax = subgraph.p_nonmax.kidx
        if kidx_nonmax ∈ seen_kidxs_nonmax
            # a subgraph has already been pinned with this (nonmax) kidx, so we cannot pin
            # another block in the column associated with this kidx
            subgraph.pinned = false
            continue
        else
            push!(seen_kidxs_nonmax, kidx_nonmax)
            subgraph.pinned = true
        end
    end

    # pin additional subgraphs via filter-criteria
    multiple_irrep_filter!.(subgraphs)
    fake_weyl_filter!.(subgraphs)

    return subgraphs
end
function subgraph_permutations(subgraph)
    # if the subgraph is pinned, we do not generate any permutations
    if subgraph.pinned
        As = [subgraph.A]
        return SubGraphPermutations(subgraph.p_max, subgraph.p_nonmax, As, #=pinned=# true)
    end
    
    # for each set of same-irrep nodes, we can perform all possible permutations of the
    # columns in the adjacency matrix associated with those nodes; if there are multiple
    # same-irrep node sets, we must perform all permutations across all sets
    multiples = subgraph.p_nonmax.multiples   
    pss = collect.(permutations.(multiples))
    Np = prod(length, pss) # aggregate number of same-irrep permutations, across irreps
    As = Vector{Matrix{Int}}(undef, Np)
    for (j,is) in enumerate(CartesianIndices(Tuple(Base.OneTo.(length.(pss)))))
        cols = reduce(vcat, [pss[k][i] for (k,i) in enumerate(Tuple(is))])
        As[j] = subgraph.A[:,cols]
    end
    return SubGraphPermutations(subgraph.p_max, subgraph.p_nonmax, As, #=pinned=# false)
end

function multiple_irrep_filter!(subgraph)
    # if there's only a single type of irrep included in the considered maximal k-manifold,
    # there is no point in considering variations from it as they are indistinguishable;
    # similarly if there is only a single type of nonmaximal k-manifold
    if length(subgraph.p_max.multiples) == 1 || length(subgraph.p_nonmax.multiples) == 1
        subgraph.pinned = true
    end
end
function fake_weyl_filter!(subgraph)
    # TODO
end

function permute_subgraphs(subgraphs)
    # find pinned subgraphs: only include a single, canonical, trivial permutation then
    set_pinned_subgraphs!(subgraphs)

    # generate all subgraph permutations of non-pinned subgraphs
    return subgraph_permutations.(subgraphs)
end
# TODO: change output structure of the above so that it becomes easier to generate the
#       different related graphs and also implement the graph filters. Then do that.

# TODO: some entries in `subductionsd[n]` contain trivial connections via the general point
#       Ω: we should probably just remove any nodes associated with Ω. In fact, there are
#       other examples of this where Bilbao/ISOTROPY include "special" k-points which
#       nevertheless have a trivial little group (presumably because it has the composition
#       of inversion and time-reversal?)

## ----------------------------------------------------------------------------------------
## Unfolding a band graph `g` via `kg_trail`
using LayeredLayouts # must use the equal_layers branch at github.com/thchr/LayeredLayouts.jl

function partition_ids(p :: Partition)
    # return labels of the form `(irlab, n)` where `n` is an incrementing counter from 1,
    # which increments for every `irlab` that is repeated.
    ids = [(label(p.lgirs[idxs[n]]), n) for idxs in p.multiples for n in eachindex(idxs)]
    return ids
end
function unfold_bandgraph(
            subgraphs,
            partitions,
            kg #= info about connections between k-partitions, from `partition_graph` =#
            )
    kg′ = split_nonmaximal_nodes(kg)
    kg_trail = chinese_postman(kg′)

    partitionsd = Dict(p.klab => p for p in partitions)
    k_ids = Dict(p.klab => 0 for p in partitions)
    g_trail = MetaGraph(DiGraph();
                Label = Tuple{String, Int, Int}, # (irrep label, multiplicity identifier, path-identifier)
                VertexData = @NamedTuple{lgir::LGIrrep, maximal::Bool, trailidx::Int},
                EdgeData   = @NamedTuple{weight::Int})
    # --- add vertices ---
    for (trailidx, kidx′) in enumerate(kg_trail.trail)
        klab = label_for(kg′, kidx′)[1]
        p = partitionsd[klab]
        lgir_ids = partition_ids(p)
        k_id = (k_ids[klab] += 1)
        for (idx, lgir_id) in enumerate(lgir_ids)
            vertex_id = (lgir_id..., k_id)
            add_vertex!(g_trail, vertex_id, (; lgir=p.lgirs[idx],
                                               maximal=p.maximal,
                                               trailidx=trailidx))
        end
    end

    # --- add edges ---
    subgraphsd  = Dict((s.p_max.klab, s.p_nonmax.klab) => s for s in subgraphs)
    foreach(klab->k_ids[klab]=0, keys(k_ids)) # reset counters
    # peel one iteration
    kidx′ = first(kg_trail.trail)
    p_prev = partitionsd[label_for(kg′, kidx′)[1]]
    k_id = (k_ids[p_prev.klab] += 1)
    for kidx′ in @view kg_trail.trail[2:end]
        vertex_ids_prev = [(lgir_id..., k_id) for lgir_id in partition_ids(p_prev)]

        p_curr = partitionsd[label_for(kg′, kidx′)[1]]
        k_id = (k_ids[p_curr.klab] += 1)
        vertex_ids_curr = [(lgir_id..., k_id) for lgir_id in partition_ids(p_curr)]

        ismax2nonmax = p_prev.maximal # we always go from max → nonmax or nonmax → max
        klabs = ismax2nonmax ? (p_prev.klab, p_curr.klab) : (p_curr.klab, p_prev.klab)
        subgraph = subgraphsd[klabs]

        vertex_ids_max, vertex_ids_nonmax = ismax2nonmax ? (vertex_ids_prev, vertex_ids_curr) : (vertex_ids_curr, vertex_ids_prev)
        for (j,a) in enumerate(eachcol(subgraph.A)) # j: local column index (nonmax) in block/subgraph
            i = something(findfirst(≠(0), a)) # local row index (max) in block/subgraph
            weight = a[i]
            vertex_id_max, vertex_id_nonmax = vertex_ids_max[i], vertex_ids_nonmax[j]
            vertex_id_prev, vertex_id_curr = ismax2nonmax ? (vertex_id_max, vertex_id_nonmax) : (vertex_id_nonmax, vertex_id_max)
            add_edge!(g_trail, vertex_id_prev, vertex_id_curr, (; weight=weight))
        end

        p_prev = p_curr
        vertex_ids_curr = vertex_ids_prev
    end

    klabs_trail = [kg′[label_for(kg′, code)].klab for code in kg_trail.trail]
    return g_trail, klabs_trail
end


end # module