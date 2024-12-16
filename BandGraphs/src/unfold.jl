## Unfolding a band graph `g` via `kg_trail`

function partition_ids(p :: Partition)
    # return labels of the form `(irlab, n)` where `n` is an incrementing counter from 1,
    # which increments for every `irlab` that is repeated.
    ids = [(label(p.lgirs[idxs[n]]), n) for idxs in p.multiples for n in eachindex(idxs)]
    return ids
end

unfold_bandgraph(bandg::BandGraph, kg) = unfold_bandgraph(bandg.subgraphs, bandg.partitions, kg)
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
        label_type       = Tuple{String, Int, Int}, # (irrep label, multiplicity identifier, path-identifier)
        vertex_data_type = @NamedTuple{lgir::LGIrrep, maximal::Bool, trailidx::Int},
        edge_data_type   = @NamedTuple{weight::Int})
    
    # --- add vertices ---
    for (trailidx, kidx′) in enumerate(kg_trail.trail)
        if kidx′ == -1
            # we've must have more than one disconnected component in the graph; so there
            # is not a single disconnected trail we can follow between the vertices; that is
            # okay though, there is just no vertex to add to `g_trail` here - and we need to
            # keep the absence of a point here in mind when we get to plotting later
            continue
        end
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
    trailidx_begin = firstindex(kg_trail.trail)
    # peel one iteration
    @label CONNECTED_COMPONENT_START
    kidx′ = kg_trail.trail[trailidx_begin]
    p_prev = partitionsd[label_for(kg′, kidx′)[1]]
    k_id = (k_ids[p_prev.klab] += 1)
    for trailidx in trailidx_begin+1:lastindex(kg_trail.trail)
        kidx′ = kg_trail.trail[trailidx]
        if kidx′ == -1
            # nothing to add for trail transitions between disconnected components; restart
            # the loop at the next trailidx
            trailidx_begin = trailidx + 1
            @goto CONNECTED_COMPONENT_START
        end
        vertex_ids_prev = [(lgir_id..., k_id) for lgir_id in partition_ids(p_prev)]

        p_curr = partitionsd[label_for(kg′, kidx′)[1]]
        k_id = (k_ids[p_curr.klab] += 1)
        vertex_ids_curr = [(lgir_id..., k_id) for lgir_id in partition_ids(p_curr)]

        ismax2nonmax = p_prev.maximal # we always go from max → nonmax or nonmax → max
        klabs = ismax2nonmax ? (p_prev.klab, p_curr.klab) : (p_curr.klab, p_prev.klab)
        subgraph = subgraphsd[klabs]

        vertex_ids_max, vertex_ids_nonmax = ismax2nonmax ? (vertex_ids_prev, vertex_ids_curr) : (vertex_ids_curr, vertex_ids_prev)
        for (j,a) in enumerate(eachcol(subgraph.A)) # j: local column index (nonmax) in block/subgraph
            i = 0 # local row index (max) in block/subgraph
            while (i = findnext(≠(0), a, i+1); !isnothing(i))
                weight = a[something(i)]
                vertex_id_max, vertex_id_nonmax = vertex_ids_max[i], vertex_ids_nonmax[j]
                vertex_id_prev, vertex_id_curr = ismax2nonmax ? (vertex_id_max, vertex_id_nonmax) : (vertex_id_nonmax, vertex_id_max)
                add_edge!(g_trail, vertex_id_prev, vertex_id_curr, (; weight=weight))
            end
            i == 0 && error("unexpectedly encountered zero-column: nonmax irrep unconnected to irrep in max partition")
        end

        p_prev = p_curr
        vertex_ids_curr = vertex_ids_prev
    end

    klabs_trail = [code == -1 ? "" : kg′[label_for(kg′, code)].klab for code in kg_trail.trail]
    return g_trail, klabs_trail
end