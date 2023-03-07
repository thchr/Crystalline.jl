# The method `chinese_postman` provides a solution to the Chinese postman problem for
# undirected graphs (https://en.wikipedia.org/wiki/Chinese_postman_problem). Unless a
# Eulerian cycle/trail exists, the trail will generally include overlapping parts.
# The solution is not provided as a cycle in general, but simply as the minimum-length trail
# (i.e., start and end points will in general differ).
# Implementation follows https://www.datacamp.com/tutorial/networkx-python-graph-tutorial.
# NB: At one point in the implementation, I exploit that we will never create a multigraph
#     when adding "fake" paths; this is using a particular property of the kinds of graphs
#     we use here. This should be generalized if the method is intended as a general-purpose
#     solver of the Chinese postman problem.
# NB: The solver presently assumes unit-length distances for each edge; that too could be
#     generalized.

using Graphs
using Graphs: AbstractSimpleGraph
using GraphsMatching
include("eulerian.jl")

struct Trail{T}
    nv :: T
    trail :: Vector{T}
end

function chinese_postman(
            g::AbstractSimpleGraph{T},
            check_connectedness::Bool=true
            ) where T
    # if disconnected, split into connected components and run on each component separately,
    # concatenating trails and indicating disconnections between by -1 elements
    if check_connectedness && !is_connected(g)
        return _chinese_postman_split_disconnected_components(g)
    end

    # --- Step 2 ---
    # Step 2.1. Find all vertices of odd order
    odd_verts = find_odd_vertices(g)

    # Stage 2.1*. Short-circuit to plain Eulerian trail if we have zero or two odd vertices
    if isempty(odd_verts)
        return tour_from_cycle(eulerian(g, 1, 1))
    elseif length(odd_verts) == 2
        return tour_from_cycle(eulerian(g, odd_verts[1], odd_verts[2]))
    end

    # Step 2.2. Build all pairwise connections between odd-order vertices & compute shortest
    #           distances between these pairs
    pairs = build_all_pairs(odd_verts) # edges between odd pairs
    pairs_shortest_paths = [a_star(g, e.src, e.dst) for e in pairs]

    # Step 2.3-2.4. Compute a minimal matching of the complete graph associated with
    #               these pairs (returning a set of pairs that cover all odd-order vertices,
    #               while minimizing the overall distance between pairs)
    idxs_in_shortest_paths, pairs_dists = find_minimal_matching(pairs_shortest_paths)

    # Step 2.4*. Of the above additional pair paths, we will, if possible, remove the
    #            longest-length pair path, which will leave two odd degree nodes remaining
    #            that will then be the start and end points of the Eulerian trail
    rm_pair_idx = argmax(pairs_dists)
    longest_pair_path = pairs_shortest_paths[idxs_in_shortest_paths[rm_pair_idx]]
    start, stop = src(first(longest_pair_path)), dst(last(longest_pair_path))
    deleteat!(idxs_in_shortest_paths, rm_pair_idx)
    
    # Step 2.4**. Build a dictionary of start/stop points in added pairs and their actual
    #             "real" paths
    pair_paths_d = Dict{Edge{T}, Vector{Edge{T}}}()
    for idx in idxs_in_shortest_paths
        path = pairs_shortest_paths[idx]
        e = Edge(src(first(path)), dst(last(path)))
        pair_paths_d[e] = path
    end

    # Step 2.5. Augment original graph by including new paths from above. Note that in our
    #           graph construction, the addition of new paths will never make it a multi
    #           graph since the intermediate k-points will always be bypassed
    g′, added_edge_d = augment_graph(g, pairs_shortest_paths, idxs_in_shortest_paths)

    # --- Step 3 ---
    # Step 3.1. Compute Eulerian cycle
    cycle = eulerian(g′, start, stop)

    # Step 3.2. Reconstruct "real" path by substituting out added edges
    return reconstruct_real_tour(cycle, added_edge_d, pair_paths_d)
end

function _chinese_postman_split_disconnected_components(g::AbstractSimpleGraph{T}) where T
    components = connected_components(g)
    filter!(length(verts)==1, components) # remove single-vertex components

    trail = T[]
    for (n, verts) in enumerate(components)
        # copy into new graph
        g′ = SimpleGraph{T}(length(verts))
        idxmap_lin2orig = verts
        idxmap_orig2lin = Dict{T,T}(v => i for (i,v) in enumerate(verts))
        for e in edges(g)
            if src(e) ∈ verts # don't need to check `dst`, cf. connectedness
                s, d = src(e), dst(e)
                e′ = Edge(idxmap_orig2lin[s], idxmap_orig2lin[d])
                add_edge!(g′, e′) || error("failed to add edge in copying")
            end
        end
        t′ = chinese_postman(g′, false)
        append!(trail, idxmap_lin2orig[t′.trail])
        # use `-1` as sentinel to indicate new disconnected trail
        n ≠ length(components) && push!(trail, -1)
    end
    return Trail{T}(nv(g), trail)
end

find_odd_vertices(g) = findall(v -> isodd(degree(g, v)), vertices(g))

build_all_pairs(vs) = [Edge(vs[i], vs[j]) for i in eachindex(vs) for j in i+1:length(vs)]

function find_minimal_matching(paths)
    # build complete graph and associated edge weights
    unique_path_verts = sort!(unique(vcat(src.(first.(paths)), dst.(last.(paths)))))
    vert2idx = Dict(v => idx for (idx, v) in enumerate(unique_path_verts))
    idx2vert = Dict(v => idx for (idx, v) in vert2idx)
    
    es = [Edge(vert2idx[first(path).src], vert2idx[last(path).dst]) for path in paths]
    cg = SimpleGraph(es)

    w = Dict{Edge{eltype(first(es))}, Float64}() # edge weights (distances)
    for (e, path) in zip(es, paths)
        d = length(path)
        w[e] = d
    end

    # solve minimum weight perfect matching problem (i.e., find list of pairs of the
    # vertices in `cg` that cover all vertices and whose edge weights sum to the smallest
    # possible value, when summed over the pairs)
    mates = if all(==(first(values(w))), values(w))
        # Work around bug in GraphsMatching.jl for 2-node or same-weight graphs
        # (https://github.com/JuliaGraphs/GraphsMatching.jl/issues/5)
        collect(reverse(vertices(cg)))
    else
        m = minimum_weight_perfect_matching(cg, w)
        m.mate
    end

    # extract the set of `nv(cg)/2` pairs from the matching
    pairs = unique([Edge(min(src, dst), max(src, dst)) for (src, dst) in enumerate(mates)])
    
    # find location of each pair in `paths`
    pair_idxs_in_paths = map(pairs) do e
        idx = findfirst(paths) do path
            v₁, v₂ = idx2vert[src(e)], idx2vert[dst(e)]
            p₁, p₂ = src(first(path)), dst(last(path))
            (v₁ == p₁ && v₂ == p₂) || (v₂ == p₁ && v₁ == p₂)
        end
        something(idx)
    end

    # find length of above paths
    pair_path_lengths = map(idx -> length(paths[idx]), pair_idxs_in_paths)

    return pair_idxs_in_paths, pair_path_lengths
end

function augment_graph(g::AbstractSimpleGraph{T}, paths, idxs) where T
    g′ = SimpleGraph{T}(nv(g))
    added_edge_d = Dict{Edge{T}, Bool}() # whether edges are original (0) or fictitious (1)
    for e in edges(g)
        add_edge!(g′, e) || error("edge inclusion was unsuccesful")
        added_edge_d[e] = false
    end
    for idx in idxs
        path = paths[idx]
        e = Edge(src(first(path)), dst(last(path)))
        add_edge!(g′, e) || error("edge inclusion was unsuccesful")
        added_edge_d[e] = true
    end
    return g′, added_edge_d
end

function reconstruct_real_tour(cycle, added_edge_d::Dict, pair_paths_d::Dict)
    # start by computing number of valid vertices in graph; could've used `g`, but nice not
    # needing to pass it around
    nv = maximum(added_edge_d) do (e, was_added)
        was_added ? 0 : max(src(e), dst(e))
    end

    # now we follow `cycle` through its edges, substituting out any fictituous path by its
    # actual path through the original graph
    tmp = iterate(cycle)
    isnothing(tmp) && error("empty cycle")
    previous, state = tmp
    trail = [previous]
    while true
        tmp = iterate(cycle, state)
        isnothing(tmp) && break
        current, state = tmp
        e = Edge(previous, current)
        was_added = haskey(added_edge_d, e) ? added_edge_d[e] : added_edge_d[reverse(e)]
        if !was_added
            push!(trail, current)
        else
            path = if haskey(pair_paths_d, e)
                pair_paths_d[e]
            else
                path′ = pair_paths_d[reverse(e)] # path goes in wrong direction
                [reverse(e) for e in Iterators.reverse(path′)]
            end
            for e′ in path
                # `path` is a sequential list of edges; add just the `dst` at each step
                push!(trail, dst(e′))
            end
        end
        previous = current
    end

    return Trail(nv, trail)
end

function tour_from_cycle(cycle)
    nv = maximum(cycle)
    return Trail(nv, cycle)
end


function Graphs.DiGraph(t::Trail)
    dg = DiGraph(t.nv)
    trail = t.trail
    return _digraph!(dg, trail)
end

function _digraph!(dg::DiGraph, trail, state=firstindex(trail))
    tmp = iterate(trail, state)
    isnothing(trail) && error("trail is empty")
    previous, state = tmp
    while true
        tmp = iterate(trail, state)
        isnothing(tmp) && break
        current, state = tmp
        if current == -1
            # there is a boundary between disconnected paths; start on new path
            return _digraph!(dg, trail, state)
        end
        add_edge!(dg, Edge(previous, current))
        previous = current
    end
    return dg
end