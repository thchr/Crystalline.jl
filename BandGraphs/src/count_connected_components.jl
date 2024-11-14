# Copy of https://github.com/JuliaGraphs/Graphs.jl/pull/407, until it is merged:
# implements `count_connected_components(g, [label, search_queue])`.
# TODO: remove when merged

"""
    _connected_components!(label, g, [search_queue])

Fill `label` with the `id` of the connected component in the undirected graph
`g` to which it belongs. Return a vector representing the component assigned
to each vertex. The component value is the smallest vertex ID in the component.

A `search_queue`, an empty `Vector{eltype(edgetype(g))}`, can be provided to reduce
allocations if `_connected_components!` is intended to be called multiple times sequentially.
If not provided, it is automatically instantiated.

### Performance
This algorithm is linear in the number of edges of the graph.
"""
function _connected_components!(
    label::AbstractVector, g::AbstractGraph{T}, search_queue::Vector{T}=Vector{T}()
) where {T}
    isempty(search_queue) || error("provided `search_queue` is not empty")
    for u in vertices(g)
        label[u] != zero(T) && continue
        label[u] = u
        push!(search_queue, u)
        while !isempty(search_queue)
            src = popfirst!(search_queue)
            for vertex in all_neighbors(g, src)
                if label[vertex] == zero(T)
                    push!(search_queue, vertex)
                    label[vertex] = u
                end
            end
        end
    end
    return label
end

"""
    count_connected_components( g, [label, search_queue]; reset_label::Bool=false)

Return the number of connected components in `g`.

Equivalent to `length(connected_components(g))` but uses fewer allocations by not
materializing the component vectors explicitly. Additionally, mutated work-arrays `label`
and `search_queue` can be provided to reduce allocations further (see
[`_connected_components!`](@ref)).

## Keyword arguments
- `reset_label :: Bool` (default, `false`): if `true`, `label` is reset to zero before
  returning.

## Example
```
julia> using Graphs

julia> g = Graph(Edge.([1=>2, 2=>3, 3=>1, 4=>5, 5=>6, 6=>4, 7=>8]));

length> connected_components(g)
3-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [4, 5, 6]
 [7, 8]

julia> count_connected_components(g)
3
```
"""
function count_connected_components(
    g::AbstractGraph{T},
    label::AbstractVector=zeros(T, nv(g)),
    search_queue::Vector{T}=Vector{T}();
    reset_label::Bool=false
) where T
    _connected_components!(label, g, search_queue)
    c = count_unique(label)
    reset_label && fill!(label, zero(eltype(label)))
    return c
end

function count_unique(label::Vector{T}) where T
    seen = Set{T}()
    c = 0
    for l in label
        if l ∉ seen
            push!(seen, l)
            c += 1
        end
    end
    return c
end