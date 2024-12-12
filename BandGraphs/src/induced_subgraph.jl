function Graphs.induced_subgraph(bandg :: BandGraph{D}, vlist) where D
    partitions = bandg.partitions
    issorted(vlist) || (vlist = sort(vlist))

    # determine which partitions correspond to the vertex indices in `vlist`
    # (this aligns with the vertex-indexing ordering from `assemble_graph`)
    kidx = 1
    irmax = 0
    code = 0
    partitions′ = empty(partitions)
    vmap = Int[]
    lgirs_idxs = Int[]
    for p in partitions
        for j in eachindex(p.lgirs)
            code += 1
            if insorted(code, vlist)
                push!(lgirs_idxs, j)
                push!(vmap, code)
            end
        end
        isempty(lgirs_idxs) && continue
        lgirs = p.lgirs[lgirs_idxs]
        multiples = _identify_multiples(lgirs, first(lgirs))
        p′ = Partition{D}(p.klab, lgirs, multiples, p.maximal, kidx, irmax+1:irmax+length(lgirs))
        push!(partitions′, p′)

        kidx += 1
        irmax = last(p′.iridxs)
        empty!(lgirs_idxs)
    end
    
    # determine the subgraphs connected by elements in vlist (exploit that this is a
    # bipartite graph)
    subgraphs = bandg.subgraphs
    subgraphs′ = empty(subgraphs)
    for subgraph in subgraphs
        klab_max    = subgraph.p_max.klab
        klab_nonmax = subgraph.p_nonmax.klab
        idx_max = findfirst(p-> p.klab==klab_max, partitions′)
        isnothing(idx_max) && continue
        idx_nonmax = findfirst(p -> p.klab==klab_nonmax, partitions′)
        isnothing(idx_nonmax) && continue

        # determine include parts of the subgraph adjacency matrix
        p_max′    = partitions′[idx_max]
        p_nonmax′ = partitions′[idx_nonmax]
        A = subgraph.A
        A′ = Matrix{Int}(undef, length(p_max′.lgirs), length(p_nonmax′.lgirs))
        for (local_idx_max′, iridx_max′) in enumerate(p_max′.iridxs)
            for (local_idx_nonmax′, iridx_nonmax′) in enumerate(p_nonmax′.iridxs)
                iridx_max = vmap[iridx_max′]
                iridx_nonmax = vmap[iridx_nonmax′]
                local_idx_max = something(findfirst(==(iridx_max), subgraph.p_max.iridxs))
                local_idx_nonmax = something(findfirst(==(iridx_nonmax), subgraph.p_nonmax.iridxs))
                A′[local_idx_max′, local_idx_nonmax′] = A[local_idx_max, local_idx_nonmax]
            end
        end
        subgraph′ = SubGraph(p_max′, p_nonmax′, A′)
        push!(subgraphs′, subgraph′)
    end
    
    bandg′ = BandGraph(subgraphs′, partitions′)
    return bandg′, vmap
end

function _identify_multiples!(multiples :: Vector{UnitRange{Int}}, lgirs, lgir_target, start)
    # we assume that `lgirs` is sorted in the sense that identical `lgirs` are adjacent
    stop = start
    for j in start+1:length(lgirs)
        lgir = lgirs[j]
        if label(lgir) == label(lgir_target)
            stop = j
        else
            push!(multiples, start:stop)
            return _identify_multiples!(multiples, lgirs, lgirs[j], j)
        end
    end
    push!(multiples, start:stop)
    return multiples
end
function _identify_multiples(lgirs, lgir_target, start::Int=1)
    _identify_multiples!(UnitRange{Int}[], lgirs, lgir_target, start)
end

# obtain the induced graph that contains all vertices in the partitions associated with each
# k-label in klabs; think of this as cutting out the parts of a band graph that contain the
# vertices in klabs
function partition_subset_induced_subgraph(bandg :: BandGraph{D}, klabs) where D
    # determine the vertex indices to include from the klabel-indices
    vlist = Int[]
    for klab in klabs
        for p in bandg.partitions
            if p.klab == klab
                append!(vlist, p.iridxs)
            end
        end
    end

    # extract associated subgraph
    return induced_subgraph(bandg, vlist)[1]
end

function subgraph_subset_induced_subgraph(
    bandg_or_bandgp :: Union{BandGraph{D}, BandGraphPermutations{D}},
    subgraph_idxs :: Vector{Int}
    ) where D
    @assert allunique(subgraph_idxs)
    partitions = bandg_or_bandgp.partitions

    partitions′ = empty(partitions)
    p_max_kidxs_from_subgraph_idxs = Vector{Int}(undef, length(subgraph_idxs))
    p_nonmax_kidxs_from_subgraph_idxs = Vector{Int}(undef, length(subgraph_idxs))
    # extract partitions featured in `subgraph_idxs`; max ones first, then nonmax ones
    kidx = 1
    irmax = 0
    seen_kidxs = Set{Int}()
    for (j, s_idx) in enumerate(subgraph_idxs) # max partitions
        s = get_subgraph(bandg_or_bandgp, s_idx)
        if s.p_max.kidx ∈ seen_kidxs
            prev_kidx = findfirst(p -> p.klab == s.p_max.klab, @view partitions′[1:kidx-1])
            p_max_kidxs_from_subgraph_idxs[j] = prev_kidx
            continue
        end
        idx = findfirst(partitions) do p
            s.p_max.kidx == p.kidx
        end
        p = partitions[something(idx)]
        p′ = Partition{D}(p.klab, p.lgirs, p.multiples, p.maximal, kidx, irmax+1:irmax+length(p.lgirs))

        push!(partitions′, p′)
        push!(seen_kidxs, p.kidx)
        p_max_kidxs_from_subgraph_idxs[j] = kidx

        kidx += 1
        irmax = last(p′.iridxs)
    end
    N = length(partitions′) # number of maximal partitions in new graph

    for (j, s_idx) in enumerate(subgraph_idxs) # nonmax partitions
        s = get_subgraph(bandg_or_bandgp, s_idx)
        if s.p_nonmax.kidx ∈ seen_kidxs
            prev_kidx = findfirst(p -> p.klab == s.p_nonmax.klab, @view partitions′[N+1:kidx-1])
            p_nonmax_kidxs_from_subgraph_idxs[j] = N+prev_kidx
            continue
        end
        idx = findfirst(partitions) do p
            s.p_nonmax.kidx == p.kidx
        end
        p = partitions[something(idx)]
        p′ = Partition{D}(p.klab, p.lgirs, p.multiples, p.maximal, kidx, irmax+1:irmax+length(p.lgirs))

        push!(partitions′, p′)
        push!(seen_kidxs, p.kidx)
        p_nonmax_kidxs_from_subgraph_idxs[j] = kidx

        kidx += 1
        irmax = last(p′.iridxs)
    end

    # build new subgraphs or permuted subgraphs corresponding to those featured in
    # `subgraph_idxs`; important thing is to get the new `iridxs` and updated partitions
    T = bandg_or_bandgp isa BandGraph ? SubGraph{D} : #=BandGraphPermutations{D}}=# SubGraphPermutations{D}
    subgraphs′ = Vector{T}(undef, length(subgraph_idxs))
    for (j, s_idx) in enumerate(subgraph_idxs)
        s = get_subgraph(bandg_or_bandgp, s_idx)
        p_max = partitions′[p_max_kidxs_from_subgraph_idxs[j]]
        p_nonmax = partitions′[p_nonmax_kidxs_from_subgraph_idxs[j]]
        s′ = SubGraph{D}(p_max, p_nonmax, s.A, s.monodromy_tie_via_Ω)
        if bandg_or_bandgp isa BandGraphPermutations
            s_ps = bandg_or_bandgp.subgraphs_ps[s_idx]
            # need to update `chains_others` indices to new indices in `subgraphs′`
            chains_others = s_ps.chains_others
            chains_others′ = if isnothing(chains_others)
                chains_others
            else 
                map(s_idx′->something(findfirst(==(s_idx′), subgraph_idxs)), chains_others)
            end
            s_ps′ = SubGraphPermutations{D}(s′, s_ps.permutations, s_ps.cols_or_rows,
                                         s_ps.pinned, s_ps.is_chained, chains_others′)
            subgraphs′[j] = s_ps′
        else
            subgraphs′[j] = s′
        end
        
    end
    #=
    for (j, s_idx) in enumerate(subgraph_idxs)
        println(get_subgraph(bandg_or_bandgp, s_idx).p_max.klab => 
                get_subgraph(bandg_or_bandgp, s_idx).p_nonmax.klab, ": ",
                partitions′[p_max_kidxs_from_subgraph_idxs[j]].klab =>
                partitions′[p_nonmax_kidxs_from_subgraph_idxs[j]].klab)
    end
    =#

    if bandg_or_bandgp isa BandGraph
        return BandGraph{D}(subgraphs′, partitions′)
    elseif bandg_or_bandgp isa BandGraphPermutations
        return BandGraphPermutations{D}(partitions′, subgraphs′)
    else
        error("unreachable")
    end
end

get_subgraph(bandg :: BandGraph, i::Int) = bandg.subgraphs[i]
get_subgraph(bandgp :: BandGraphPermutations, i::Int) = bandgp.subgraphs_ps[i].subgraph