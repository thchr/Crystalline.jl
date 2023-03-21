# ---------------------------------------------------------------------------------------- #

# While the graph associated with `A` in principle could be constructed from the `Graph(A)`
# this loses a lot of relevant information. Instead, we construct it directly from
# the subgraph and partition structures, and return additional metadata on vertices & edges
function assemble_graph(subgraphs, partitions)
    g = MetaGraph(Graph();
            label_type       = Tuple{String, Int},
            vertex_data_type = @NamedTuple{lgir::LGIrrep, maximal::Bool},
            edge_data_type   = @NamedTuple{weight::Int})
    # add vertices to graph
    for p in partitions
        j = 0
        for codes′ in p.multiples
            for n in eachindex(codes′)
                j += 1
                lgir = p.lgirs[j]
                # to distinguish multiple identical irrep labels, we include a incrementing
                # number (1, 2, …) along with the irrep label as the overall vertex label,
                # see also `partition_ids(::Partition)`
                vertex_label = (label(lgir), n)
                add_vertex!(g, vertex_label, (; lgir=lgir, maximal=p.maximal))
            end
        end
    end
    @assert Graphs.nv(g) == last(partitions[end].iridxs) # check number of vertices

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

# ---------------------------------------------------------------------------------------- #

function partition_graph(subgraphs, partitions)
    g = MetaGraph(Graph();
            label_type       = String,
            vertex_data_type = @NamedTuple{maximal::Bool})
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

# ---------------------------------------------------------------------------------------- #

# compute a modified graph from `kg` such that any nonmaximal node is only ever entered and
# exited by one edge; i.e., ensure the degree of nonmaximal nodes is 2. When the degree `n`
# is higher, we split the node into ∑ⱼ₌₁ⁿ⁻¹ j = (n-1)n/2 duplicated nodes. The resulting 
# graph will necessarily have a chinese postman solution (possibly for each disconnected
# component of `kg`) that can be computed without needing a multigraph implementation;
# effectively, this is an alternative to multigraphs.
function split_nonmaximal_nodes(kg::MetaGraph)
    kg′ = MetaGraph(Graph();
                    label_type       = Tuple{String, Int},
                    vertex_data_type = @NamedTuple{klab::String, code::Int, maximal::Bool})
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