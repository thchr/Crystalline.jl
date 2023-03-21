## Unfolding a band graph `g` via `kg_trail`

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
        label_type       = Tuple{String, Int, Int}, # (irrep label, multiplicity identifier, path-identifier)
        vertex_data_type = @NamedTuple{lgir::LGIrrep, maximal::Bool, trailidx::Int},
        edge_data_type   = @NamedTuple{weight::Int})
    
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


# WIP below: aim was to remove the nonmaximal vertices and replace them with edge info; but
#            stalled due to problem in "NB" below.
function edgify_nonmax_vertices(g_trail)
    # FIXME: This does not work in general since the general case would require us to have a
    #        multigraph structure. To see this, consider e.g. the case that we have an irrep
    #        Γ₁ which splits into Δ₁+Δ₂ and then recombines to X₁. When we remove the edges
    #        due to Λ₁ and Λ₂ and seek to make them direct edges between Γ₁ and X₁, we end
    #        up needing _two_ edges between the Γ₁ and X₁ vertices

    g_decimated = MetaGraph(DiGraph();
        label_type       = Tuple{String, Int, Int}, # (irrep label, multiplicity identifier, path-identifier)
        vertex_data_type = @NamedTuple{lgir::LGIrrep, maximal::Bool, trailidx::Int},
        edge_data_type   = @NamedTuple{lgir::LGIrrep, weight::Int}
        )
    original_idxs = Int[]

    for v in vertices(g_trail)
        lᵥ = label_for(g_trail, v)
        vertex_data = g_trail[lᵥ]
        vertex_data.maximal || continue
        add_vertex!(g_decimated, lᵥ, vertex_data)
        push!(original_idxs, v)
    end

    for v in original_idxs
        lᵥ = label_for(g_trail, v)
        ns = outneighbors(g_trail, v)
        for n in ns # over nonmaximal nodes
            lₙ = label_for(g_trail, n)
            m = only(outneighbors(g_trail, n)) # nonmax nodes can only have one exiting edge
            lₘ = label_for(g_trail, m)
            edge_data = (; lgir=g_trail[lₙ].lgir, g_trail[lᵥ,lₙ]...)
            add_edge!(g_decimated, lᵥ, lₘ, edge_data) || error("edge insertion failure")
        end
    end

    # original_idxs is useful to e.g., index into the `xy` struct
    return g_decimated, original_idxs
end