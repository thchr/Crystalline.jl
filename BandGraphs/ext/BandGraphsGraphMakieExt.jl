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
using Crystalline: AbstractSymmetryVector
using Graphs
using MetaGraphsNext

import BandGraphs: plot_flattened_bandgraph, make_vertices_dragable! # exported in BandGraphs

# ---------------------------------------------------------------------------------------- #

## Unfolding a band graph `g` via `kg_trail`
function plot_flattened_bandgraph(
            n::AbstractSymmetryVector{D},
            lgirsd::Dict{String, <:AbstractVector{LGIrrep{D}}};
            timereversal=true
            ) where D
    # TODO: take a `SubductionTable` as input?
    
    sgnum = num(n)
    @assert sgnum == num(group(first(first(values(lgirsd)))))

    # compute subgraphs and partitions
    bandg = build_subgraphs(n, subduction_tables(sgnum; timereversal), lgirsd)

    plot_flattened_bandgraph(bandg)
end

function plot_flattened_bandgraph(
            bandg :: BandGraph;
            xys = nothing
            )

    # compute a path `kg` that connects all required maximal and nonmaximal k manifolds
    kg = partition_graph(bandg)

    plot_flattened_bandgraph(bandg, kg; xys)
end

function plot_flattened_bandgraph(
            bandg :: BandGraph,
            kg    :: MetaGraph; # ::(!IsDirected)
            xys = nothing
            )

    # compute a directed acyclic graph `g_trail` that "unfolds" or "flattens" the original
    # graph representation associated with `subgraphs` into a frequency-momentum/dispersion
    # diagram; to unfold it, we allow ourselves to pass through certain k-manifolds multiple
    # times if necessary
    g_trail, klabs_trail = unfold_bandgraph(bandg, kg)

    # find layers in `g_trail` with the same maximal k-points and force them to have the
    # same frequency/y-position for the irreps; we do not require the same of nonmaximal
    # irrep nodes, since a notion of fixed frequency-position does not apply to those
    # manifolds (lines/planes, etc)
    if isnothing(xys)
        force_equal_layers = find_equal_maximal_layers(klabs_trail, bandg.partitions)
        force_layer = assign_layer_indices(g_trail)
        xs, ys, _ = solve_positions(Zarate(), g_trail; force_layer, force_equal_layers)
        ys = linearize_nonmaximal_y_coordinates(g_trail, ys)
    else
        xs, ys = xys
    end
    
    plot_flattened_bandgraph(g_trail, klabs_trail, xs, ys)
end

function assign_layer_indices(g_trail)
    # the x-layer of the vertices should equal the trail index of the vertex
    vs = vertices(g_trail)
    return vs .=> [g_trail[label_for(g_trail, i)].trailidx for i in vs]
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
    graphcol = Makie.RGBAf(.1,.2,.6, 1.0) # color of main elements in graph
    p = graphplot!(ax, g_trail; 
            layout = _->xy,
            arrow_size = 0,
            nlabels = nlabels,
            nlabels_color = [m ? graphcol : Makie.RGBf(.55,.55,.55) for m in maximality],
            nlabels_align = (:center,:bottom),
            nlabels_distance = 6.0,
            nlabels_attr = (;strokewidth=2, strokecolor=:white),
            # the weird transparent coloring situation below is there to avoid a bug on
            # dragging transparent / zero-size vertices (the nonmaximal vertices)
            node_color = [m ? graphcol : Makie.RGBAf(1.0,1.0,1.0,0.001) for m in maximality], 
            edge_color = graphcol,
            node_size = fill(15, length(maximality)),
            node_attr = (; strokewidth=3, 
                           strokecolor=[(m ? Makie.RGBAf(1.0,1.0,1.0,1.0) :
                                             Makie.RGBAf(1.0,1.0,1.0,0.0))
                                             for m in maximality],)
            )
    
    # --- prettification ---
    Makie.hidedecorations!(ax)  # hides ticks, grid and labels
    Makie.hidespines!(ax)
    make_vertices_dragable_bandgraph!(ax, p, g_trail)

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

        # if this is a transition point between disconnected components (`kidx == -1` which
        # maps to `klabs_trail[trailidx] = klab = ""`), there's nothing to enforce and we
        # just go to the next vertex
        isempty(klab) && continue
        
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

# the y-coordinates returned by `solve_positions` often make the paths/edges through the
# nonmaximal irrep "jagged" unnecessarily (essentially, an artifact of us treating the 
# nonmaximal irreps as nodes in the `g_trail` graph, instead of having a proper multigraph);
# we fix that here, by "linearizing" the possibly jagged edges between maximal irreps.
# note that the "frequency" of the nonmaximal irrep is not really a fixed thing in the same
# way as it is for the maximal irreps (it spans a line/plane etc; not a fixed point), so it
# is fine for us to change it ad hoc.
function linearize_nonmaximal_y_coordinates(g_trail, ys::Vector)
    Δy = -(-)(extrema(ys)...)
    ys′ = similar(ys)
    track = Dict{Int, Vector{Int}}() # trailidx => indices in `ys`
    for i in vertices(g_trail)
        l = label_for(g_trail, i)
        v = g_trail[l]
        maximal = v.maximal
        if maximal # don't change maximal y-positions
            ys′[i] = ys[i]
            continue
        end

        if !haskey(track, v.trailidx)
            track[v.trailidx] = [i]
        else
            push!(track[v.trailidx], i)
        end
    end
    for is in values(track)
        # a vertex assoc. w/ a nonmaximal irrep will only ever have two neighbors (**): it
        # cannot have more than 2 irrep neighbors, as any connected maximal irrep must have
        # equal or higher irrep dimensions; it cannot have less than 2 since it will always
        # connect two different maximal irreps
        # (**) This is true for a standard band graph, but not for a band graph that has
        #      been "split" at a maximal vertex connected to a nonmaximal vertex of irrep
        #      dimension higher than 1; ignoring this case, one could assert certain `only`
        #      outcomes below, e.g., for `(in,out)neighbors(g_trail, i)`. To allow the
        #      just-mentioned split-configuration, however, we just take averages in general
        avg_ys = Vector{Float64}(undef, length(is))
        for (idx_in_is, i) in enumerate(is)
            in_i = inneighbors(g_trail, i)   # vertex index of ingoing maximal irrep
            out_i = outneighbors(g_trail, i) # vertex index of outgoing maximal irrep
            avg_ys[idx_in_is] = (sum(@view ys[in_i])/length(in_i) + 
                                 sum(@view ys[out_i])/length(out_i))/2
        end
        for (idx_in_is, i) in enumerate(is)
            y′ = avg_ys[idx_in_is]
            if (any(≈(y′; atol=5e-2), @view avg_ys[1:idx_in_is-1]) || 
                any(≈(y′; atol=5e-2), @view avg_ys[idx_in_is+1:end]))
                # if linearizing would "collapse" multiple edges to lie on top of eachother,
                # or simply overlay their y-labels, we keep the existing positions
                if length(is) == 2 # special treatment for 2-irrep case
                    # TODO: we could do something better/more general things here by looking
                    #       again at the neighbors of the "overlapping" nonmaximal verts
                    ys′[i] = ys[i] - sign(ys[i])*Δy/3
                else
                    ys′[i] = ys[i]
                end
            else
                ys′[i] = y′
            end
        end
    end
    return ys′
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

function make_vertices_dragable_bandgraph!(ax, p, g_trail)
    Makie.deregister_interaction!(ax, :rectanglezoom)
    function node_drag_action(state, idx, event, axis)
        mouse_r = event.data

        # find related irreps, their vertex indices, and change their y-coordinates as well
        vert_id = label_for(g_trail, idx)[1:2]
        irlab, irmul = vert_id
        irlab_normalized = replace(irlab, "′"=>"") # to also move monodromy-related verts
        related_vert_idxs = findall(g_trail.vertex_labels) do l
            irlab′, irmul′ = l
            irlab_normalized == replace(irlab′, "′"=>"") && irmul == irmul′
        end
        for idx′ in related_vert_idxs
            prev_r′ = p[:node_pos][][idx′]
            p[:node_pos][][idx′] = Makie.Point(prev_r′[1], mouse_r[2])
        end
        
        # update observable
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    Makie.register_interaction!(ax, :ndrag, ndrag)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    nothing
end

end # module BandGraphsGraphMakieExt