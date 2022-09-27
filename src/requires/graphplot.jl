using .Graphs, .StaticArrays, .LinearAlgebra
import .GraphMakie
import .GraphMakie: Makie
import LayeredLayouts
#using NetworkLayout: SFDP

## -------------------------------------------------------------------------------------- ##
# Visualize graph of maximal subgroups a particular group (descendants of `g`)

# -----------------
function map_numbers_to_oneto(xs; rev::Bool=false)
    uxs = sort!(unique(xs); rev)
    idxss = [findall(==(x), xs) for x in uxs]
    ys = Vector{Int}(undef, length(xs))
    for (y,idxs) in enumerate(idxss)
        ys[idxs] .= y
    end
    return ys
end

# TODO: generalize to groups besides `SpaceGroup`?
function layout_by_order(gr::GroupRelationGraph{D,SpaceGroup{D}};
                         ymap=map_numbers_to_oneto) where D
    orders = SG_PRIMITIVE_ORDERs[D][gr.nums]
    N = length(orders)
    yposs = ymap(orders)
    xposs = Vector{Float64}(undef, N)
    maxwidth = maximum([length(findall(==(o), orders)) for o in unique(orders)])
    for o in unique(orders)
        idxs = findall(==(o), orders)
        xposs[idxs] = ((1:length(idxs)) .- (1 + (length(idxs)-1)/2)) .* (maxwidth / length(idxs))
    end
    return Makie.Point2f.(xposs, yposs), orders
end
# # TODO: requires variant of https://github.com/JuliaGraphs/NetworkLayout.jl/pull/52 to be
# #       merged
# function layout_by_order_spring(gr::GroupRelationGraph{D,SpaceGroup{D}};
#                                 ymap=map_numbers_to_oneto) where D
#     xy, orders = layout_by_order(gr; ymap)
#     orderslim = extrema(orders)
#     
#     P = eltype(xy)
#     T = eltype(P)
#     rng = Xoshiro(1)
#     initialpos = [P(p[1] + (o∉orderslim) * (rand(rng, T)-0.5)*0.05, p[2]) for (p,o) in zip(xy,orders)]
#     pin = [Point{2,Bool}(false, true) for _ in xy]
#     
#     algo = NetworkLayout.SFDPpin(; initialpos, pin, K=.001, C=500, tol=1e-1)
#     xy′ = algo(gr)
# 
#     return xy′, orders
# end
function layout_by_minimal_crossings(gr::GroupRelationGraph{D,SpaceGroup{D}};
                                     force_layer_bool=true) where D
    orders = SG_PRIMITIVE_ORDERs[D][gr.nums]
    is_supergroup_relations = gr.direction==SUPERGROUP # if true, must flip y-ordering twice
    force_layer = if force_layer_bool
        layers′ = map_numbers_to_oneto(orders; rev = is_supergroup_relations)
        maxlayer = maximum(layers′)
        1:length(orders) .=> (maxlayer+1) .- layers′
    else
        Pair{Int, Int}[]
    end
    # `solve_positions` implicitly assumes `SimpleGraph`/`SimpleDiGraph` type graphs, so we
    # must first convert `gr` to `gr′::SimpleDiGraph`
    gr′ = SimpleDiGraphFromIterator(edges(gr))
    xs, ys, _ = LayeredLayouts.solve_positions(LayeredLayouts.Zarate(), gr′; force_layer)

    is_supergroup_relations && (xs .= -xs) # re-reverse the ordering
    
    return Makie.Point2f.(ys, -xs), orders
end
# -----------------

function prettystring(x::Rational)
    io = IOBuffer()
    show(io, numerator(x))
    isone(denominator(x)) && return String(take!(io))
    print(io, "/")
    show(io, denominator(x))
    return String(take!(io))
end
function rationalizedstring(x)
    s = sprint(show, prettystring.(rationalize.(x)))
    return replace(s, "\""=>"")
end
function stringify_relations(gr::GroupRelationGraph, e::Edge)
    s, d = src(e), dst(e)
    snum, dnum = gr.nums[s], gr.nums[d]
    info = gr.infos[snum]
    didx = findfirst(child->child.num==(dnum), info.children)
    isnothing(didx) && error("edge destination not contained in graph")
    Pps = info[didx]
    if length(Pps.classes) > 1
        return "multiple conjugacy classes"
    else
        Pp = only(Pps.classes)
        P, p = Pp.P, Pp.p
        if isone(P)
            if iszero(p)
                return "identity"                                    # trivial mapping
            else
                return "p = "*rationalizedstring(p)                  # translation mapping
            end
        else
            return "P = "*rationalizedstring(P) *                    # rotation mapping
                   (iszero(p) ? "" : ", p = "*rationalizedstring(p))
        end
    end
end

"""
    plot!(ax, gr::GroupRelationGraph{D}; layout_mode)

Visualize the graph structure of `gr`.
If the Makie backend is GLMakie, the resulting plot will be interactive (on hover and drag).

See also `Makie.plot`.

## Keyword arguments
- `layout_mode`: `:fancy` or `:quick`. By default, `:fancy` is chosen for graph with less
   than 40 edges; otherwise `:quick`. The `:fancy` algorithm uses a Sugiyama-style DAG
   layout algorithm: it performs well for relatively few edges, but very poorly otherwise.

## Example

```jl
using Crystalline
using GraphMakie, GLMakie
gr = maximal_subgroups(202)
plot(gr)
```

Note that `plot(gr)` is conditionally defined only upon loading of GraphMakie (via
Requires); additionally, a backend for Makie must be explicitly loaded (here, GLMakie).
"""
function Makie.plot!(
            ax::Makie.Axis,
            gr::GroupRelationGraph{D};
            layout_mode::Symbol=(Graphs.ne(gr) ≤ 40 && Graphs.nv(gr) ≤ 20) ? :fancy : :quick
            ) where D

    mapping_strings = stringify_relations.(Ref(gr), edges(gr))

    # --- compute layout of vertices ---
    xy, orders = layout_mode == :fancy ? layout_by_minimal_crossings(gr; force_layer_bool=true) :
                 layout_mode == :quick ? layout_by_order(gr) :
                 layout_mode == :quickspring ? layout_by_order_spring(gr) :
                 error("invalid `layout_mode` value")

    # --- graphplot ---
    p = GraphMakie.graphplot!(ax, gr; 
                nlabels = string.(gr.nums),
                nlabels_distance = 6,
                nlabels_align = (:center, :bottom),
                nlabels_textsize = 20,
                nlabels_attr=(;strokecolor=:white, strokewidth=2.5),
                node_color = Makie.RGBf(.2, .5, .3),
                edge_color = fill(Makie.RGBf(.5, .5, .5), ne(gr)),
                edge_width = fill(2.0, ne(gr)),
                node_size = fill(10, nv(gr)),
                arrow_size = 0,
                elabels = fill("", ne(gr)),
                elabels_textsize = 16,
                elabels_distance = 15,
                )    
    p[:node_pos][] = xy # update layout of vertices
    
    # --- group "order" y-axis ---
    Makie.hidespines!(ax, :r, :b, :t) # remove other axes spines
    Makie.hidexdecorations!(ax, grid = true)
    unique_yidxs = unique(i->xy[i][2], eachindex(xy))
    unique_yvals = [xy[i][2] for i in unique_yidxs]
    unique_yorders = [string(orders[i]) for i in unique_yidxs]

    ax.yticks = (unique_yvals, unique_yorders)
    ax.ylabel = "group order (primitive)"
    ax.ylabelsize = 18
    ax.ygridvisible = false
    ax.ytrimspine = true

    # --- enable interactions ---
    # (following http://juliaplots.org/GraphMakie.jl/stable/generated/depgraph/)
    Makie.deregister_interaction!(ax, :rectanglezoom)
    Makie.register_interaction!(ax, :edrag, GraphMakie.EdgeDrag(p))
    Makie.register_interaction!(ax, :ndrag, GraphMakie.NodeDrag(p))
    
    es = collect(edges(gr))
    function node_hover_action(state, idx, event, axis)
        num = gr.nums[idx]
        iuclabel = iuc(num, D)
        p.nlabels[][idx] = string(num) * (state ?  " (" * iuclabel * ")" : "")
        p.node_size[][idx] = state ? 17 : 10
        p.nlabels[] = p.nlabels[]; p.node_size[] = p.node_size[]; # trigger observables

        # highlight edge descendants of highlighted nodes
        nidxs = Int[idx]
        i = 1
        while i ≤ length(nidxs)
            s = nidxs[i]
            for d in gr.fadjlist[s]
                e = Edge(s, d)
                eidx = something(findfirst(==(e), es))
                p.edge_width[][eidx] = state ? 2.75 : 2.0
                p.edge_color[][eidx] = state ? Makie.RGBf(.35, .35, .7) : Makie.RGBf(.5, .5, .5) 
                d ∈ nidxs || push!(nidxs, d)
            end
            i += 1
        end
        p.edge_width[] = p.edge_width[]; p.edge_color[] = p.edge_color[]
    end
    nodehover = GraphMakie.NodeHoverHandler(node_hover_action)
    Makie.register_interaction!(ax, :nodehover, nodehover)

    function edge_hover_action(state, idx, event, axis)
        p.elabels[][idx] = state ? mapping_strings[idx] : ""       # edge text
        p.edge_width[][idx] = state ? 4.0 : 2.0                    # edge width
        p.elabels[] = p.elabels[]; p.edge_width[] = p.edge_width[] # trigger observable
    end
    edgehover = GraphMakie.EdgeHoverHandler(edge_hover_action)
    Makie.register_interaction!(ax, :edgehover, edgehover)

    return p
end

function Makie.plot(gr::GroupRelationGraph;
                    axis = NamedTuple(), figure = NamedTuple(), kws...)
    f = Makie.Figure(resolution=(800,900); figure...)
    f[1,1] = ax = Makie.Axis(f; axis...)
    p = Makie.plot!(ax, gr; kws...)
    return Makie.FigureAxisPlot(f, ax, p)
end