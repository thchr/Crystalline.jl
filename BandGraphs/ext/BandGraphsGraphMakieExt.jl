module BandGraphsGraphMakieExt

# ---------------------------------------------------------------------------------------- #

if isdefined(Base, :get_extension)
    using GraphMakie: graphplot!, NodeDragHandler
    using GraphMakie: Makie
else
    using ..GraphMakie: graphplot!, NodeDragHandler
    using ..GraphMakie: Makie
end
using BandGraphs

using LayeredLayouts: solve_positions, Zarate # must use https://github.com/thchr/LayeredLayouts.jl#equal_layers
using Crystalline
using Graphs
using MetaGraphsNext

import BandGraphs: plot_flattened_bandgraph, make_vertices_dragable! # exported in BandGraphs

# ---------------------------------------------------------------------------------------- #

## Unfolding a band graph `g` via `kg_trail`
function plot_flattened_bandgraph(
            n::SymVector{D},
            lgirsd::Dict{String, Vector{LGIrrep{D}}}
            ) where D
    
    sgnum = num(group(first(first(n.lgirsv))))
    @assert sgnum == num(group(first(first(values(lgirsd)))))

    # compute subgraphs and partitions
    subgraphs, partitions = build_subgraphs(n, SUBDUCTIONSD[sgnum], lgirsd)

    plot_flattened_bandgraph(subgraphs, partitions)
end

function plot_flattened_bandgraph(
            subgraphs  :: AbstractVector{SubGraph{D}},
            partitions :: AbstractVector{Partition{D}};
            xys = nothing
            ) where D

    # compute a path `kg` that connects all required maximal and nonmaximal k manifolds
    kg = partition_graph(subgraphs, partitions)

    plot_flattened_bandgraph(subgraphs, partitions, kg; xys)
end

function plot_flattened_bandgraph(
            subgraphs  :: AbstractVector{SubGraph{D}},
            partitions :: AbstractVector{Partition{D}},
            kg         :: MetaGraph; # ::(!IsDirected)
            xys = nothing
            ) where D

    # compute a directed acyclic graph `g_trail` that "unfolds" or "flattens" the original
    # graph representation associated with `subgraphs` into a frequency-momentum/dispersion
    # diagram; to unfold it, we allow ourselves to pass through certain k-manifolds multiple
    # times if necessary
    g_trail, klabs_trail = unfold_bandgraph(subgraphs, partitions, kg)

    # find layers in `g_trail` with the same maximal k-points and force them to have the
    # same frequency/y-position for the irreps; we do not require the same of nonmaximal
    # irrep nodes, since a notion of fixed frequency-position does not apply to those
    # manifolds (lines/planes, etc)
    if isnothing(xys)
        force_equal_layers = find_equal_maximal_layers(klabs_trail, partitions)
        xs, ys, _ = solve_positions(Zarate(), g_trail; force_equal_layers)
    else
        xs, ys = xys
    end
    
    plot_flattened_bandgraph(g_trail, klabs_trail, xs, ys)
end

function plot_flattened_bandgraph(
            g_trail     :: MetaGraph, # ::IsDirected
            klabs_trail :: AbstractVector{<:AbstractString},
            xs          :: AbstractVector{<:Real},
            ys          :: AbstractVector{<:Real}
            )

    xy = [Makie.Point(x,y) for (x,y) in zip(xs, ys)]
    miny, maxy = extrema(ys)
    Δy = maxy - miny
    miny -= Δy/25; maxy += Δy/25

    # find the maximum number of manifolds visited in the trail; generally, we will always
    # go from a maximal (odd indices) to a non-maximal (even indices) manifold as we
    # traverse; use this below
    trailidxs = [g_trail[label_for(g_trail, i)].trailidx for i in vertices(g_trail)]
    n_visited = maximum(trailidxs)

    # --- figure and axis object ---
    f = Makie.Figure()
    ax = Makie.Axis(f[1,1])

    # --- k-manifold lines ---
    max_idxs    = 1:2:n_visited
    nonmax_idxs = 2:2:n_visited-1
    Makie.linesegments!(ax,
            [(Makie.Point(idx, miny), Makie.Point(idx, maxy)) for idx in max_idxs];
            color=Makie.RGBf(.65,.65,.65), linewidth=1) # maximal
    #Makie.linesegments!(ax,
    #       [(Point(idx, miny), Point(idx, maxy)) for idx in nonmax_idxs];
    #       color=Makie.RGBf(.9,.9,.9), linewidth=1) # non-maximal
    Makie.text!(ax, max_idxs, fill(miny-Δy/100, length(max_idxs));
        text=klabs_trail[max_idxs], align=(:center, :top), color=Makie.RGBf(.35,.35,.35))
    Makie.text!(ax, nonmax_idxs, fill(miny-Δy/100, length(nonmax_idxs));
        text=klabs_trail[nonmax_idxs], align=(:center, :top), color=Makie.RGBf(.65,.65,.65))

    # --- graph ---
    maximality = [g_trail[label_for(g_trail, i)].maximal for i in vertices(g_trail)]
    nlabels = map(vertices(g_trail)) do i
        l = label_for(g_trail, i)
        l[1]
    end
    graphcol = Makie.RGBf(.1,.2,.6) # color of main elements in graph
    p = graphplot!(ax, g_trail; 
            layout = _->xy,
            arrow_size = 0,
            nlabels = nlabels,
            nlabels_color = [m ? graphcol : Makie.RGBf(.55,.55,.55) for m in maximality],
            nlabels_align = (:center,:bottom),
            nlabels_distance = 6.0,
            nlabels_attr = (;strokewidth=2, strokecolor=:white),
            node_color = graphcol, 
            edge_color = graphcol,
            node_size = [m ? 15 : 0 for m in maximality],
            node_attr = (; strokewidth=3, strokecolor=:white)
            )
    
    # --- prettification ---
    Makie.hidedecorations!(ax)  # hides ticks, grid and labels
    Makie.hidespines!(ax)

    faxp = Makie.FigureAxisPlot(f, ax, p)
    return faxp, (; g_trail, klabs_trail, xs, ys)
end

function find_equal_maximal_layers(klabs_trail, partitions)
    # NB: we normalize monodromy-related k-labels, with the intention of forcing equal
    # frequency-position of same irreps at monodromy-related k-points. E.g., if we have
    # Γ and Γ′, we force irreps Γᵢ and Γᵢ′ to have the same frequency (... would be nice
    # to understand the physics/reasoning of this a bit better)
    klabs_trail_normalized = rstrip.(klabs_trail, '′')
    seen = Set{Int}()
    force_equal_layers = Vector{Pair{Int,Int}}()
    for (i, klab) in enumerate(klabs_trail_normalized)
        # don't need to process k-label if already discovered
        i ∈ seen && continue

        # also, check if this is a maximal manifold: if not, continue
        pidx = something(findfirst(p->p.klab == klab, partitions))
        partitions[pidx].maximal || continue
        push!(seen, i)

        # next, find all indices in `klabs_trail` with identical/monodromy-related maximal
        # k-points 
        idxs = findall(klabs_trail) do klab′
            klab′_normalized = rstrip(klab′, '′')
            klab′_normalized == klab
        end

        for i′ in @view idxs[2:end]
            push!(seen, i′)
            push!(force_equal_layers, i => i′)
        end
    end
    return force_equal_layers
end

# ---------------------------------------------------------------------------------------- #
# Plotting utilities

function make_vertices_dragable!(ax, p)
    Makie.deregister_interaction!(ax, :rectanglezoom)
    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    Makie.register_interaction!(ax, :ndrag, ndrag)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    nothing
end

end # module BandGraphsGraphMakieExt