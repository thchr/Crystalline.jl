
using Crystalline
using BandGraphs
using Graphs

## ----------------------------------------------------------------------------------------
# GraphMakie utilities
using GLMakie, GraphMakie, Colors

function make_vertices_dragable!(ax)
    deregister_interaction!(ax, :rectanglezoom)
    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)
    nothing
end

## ----------------------------------------------------------------------------------------
## Testing & visualization of band graph
using SymmetryBases

timereversal = true
sgnum = 81
sb, brs = compatibility_basis(sgnum; timereversal)
lgirsd = lgirreps(sgnum)
timereversal && (lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd))
_nv = sb[1]
n = SymVector(_nv, brs.irlabs, lgirsd)
subgraphs, partitions_max, partitions_nonmax = build_subgraphs(n, subductionsd[sgnum], lgirsd)
A = assemble_adjacency(subgraphs, partitions_max, partitions_nonmax)
L = assemble_laplacian(subgraphs, partitions_max, partitions_nonmax)
partitions = vcat(partitions_max, partitions_nonmax)
nothing

g = assemble_graph(subgraphs, partitions_max, partitions_nonmax) # structured equiv of `Graph(A)`


maximality = [g[label_for(g, i)].maximal for i in vertices(g)]
node_colors = [max ? :red : :black for max in maximality]
f, ax, p = graphplot(
    g;
    nlabels=[label(g[label_for(g, i)].lgir) for i in vertices(g)],
    nlabels_distance=4,
    node_color = node_colors, nlabels_color=node_colors,
    node_attr = (; strokewidth=4, strokecolor=:white),
    edge_color = :gray,
    layout=GraphMakie.NetworkLayout.SFDP(; C=10, K=.2, dim=2)
    )

make_vertices_dragable!(ax)
f


## ----------------------------------------------------------------------------------------
# Testing & visualization of k-connectivity graphs

kg = partition_graph(subgraphs, partitions)
kg′ = split_nonmaximal_nodes(kg)
f,ax,p=graphplot(kg′;
    nlabels = [label_for(kg′, i)[1] * " #$i" for i in vertices(kg′)],
    nlabels_distance=4,
    node_color = [kg′[label_for(kg′, i)].maximal ? :red : :black for i in vertices(kg′)],
    nlabels_color = [kg′[label_for(kg′, i)].maximal ? :red : :black for i in vertices(kg′)], 
    node_attr = (; strokewidth=5, strokecolor=:white),
    )

include("graph_routing.jl")
kg_trail = chinese_postman(kg′)
dkg = DiGraph(kg_trail)
graphplot!(ax, dkg,
    layout=_->p[:node_pos][],
    edge_color=RGB(0.2,0.2,.9),
    edge_width = 2,
    node_size = 0,
    curve_distance_usage = true)
f

# FIXME: Repeated step in trail for SG 31

## ----------------------------------------------------------------------------------------
## Unfolding a band graph `g` via `kg_trail`

g_trail, klabs_trail = unfold_bandgraph(subgraphs, partitions, kg)
trailidxs = [g_trail[label_for(g_trail, i)].trailidx for i in vertices(g_trail)]
n_trailidxs = maximum(trailidxs)

# find layers that have the same maximal k-points and force them to have the same frequency/
# y-position for the irreps; we do not require the same of nonmaximal irreps, since a notion
# of fixed frequency-position does not apply to those manifolds (lines/planes, etc).
force_equal_layers = collect(Iterators.flatten(filter(!isempty, 
    # NB: here, we normalize monodromy-related k-labels, with the intention of forcing equal
    # frequency-position of same irreps at monodromy-related k-points. E.g., if we have Γ
    # and Γ′, we force irreps Γᵢ and Γᵢ′ to have the same frequency (... would be nice to
    # understand the physics/reasoning of this a bit better)
    map(unique(rstrip.(klabs_trail, '′'))) do klab
        # check if klab is associated with a maximal manifold; if not, return empty
        if !partitions[something(findfirst(p->p.klab == klab, partitions))].maximal
            return Pair{Int,Int}[]
        end
        # find indices in `klabs_trail` of all identical/monodromy-related maximal k-points
        idxs = findall(klabs_trail) do klab′
            klab′_normalized = rstrip(klab′, '′')
            klab′_normalized == klab
        end
        println(klab, " : ", idxs)
        idx₀ = first(idxs)
        [idx₀ => idx for idx in @view idxs[2:end]]
    end)))
# TODO: must also fix Γ′ etc. correctly: what is the rule?
# TODO: remove all the nonmaximal nodes and transfer information to edges  
xs, ys, _ = solve_positions(Zarate(), g_trail;
                force_equal_layers=force_equal_layers)
xy = [Point(x,y) for (x,y) in zip(xs, ys)]
miny, maxy = extrema(ys)
Δy = maxy - miny
miny -= Δy/25; maxy += Δy/25

f, ax, _ = linesegments([(Point(idx, miny), Point(idx,maxy)) for idx in 1:2:n_trailidxs];
                color = RGBf(.5,.5,.5), linewidth=1)
text!(ax, 1:2:n_trailidxs, fill(miny-Δy/100, div(n_trailidxs+1, 2));
      text=klabs_trail[1:2:n_trailidxs], align=(:center, :top),
      color = RGBf(.25,.25,.25))
# nonmaximal k-lines
linesegments!(ax, [(Point(idx, miny), Point(idx,maxy)) for idx in 2:2:n_trailidxs-1];
      color = RGBf(.85,.85,.85), linewidth=1)
text!(ax, 2:2:n_trailidxs-1, fill(miny-Δy/100, div(n_trailidxs, 2));
      text=klabs_trail[2:2:n_trailidxs-1], align=(:center, :top),
      color = RGBf(.45,.45,.45))

# graph
maximality = [g_trail[label_for(g_trail, i)].maximal for i in vertices(g_trail)]
nlabels = map(vertices(g_trail)) do i
    l=label_for(g_trail, i)
    if g_trail[l].maximal
        l[1]#*Crystalline.supscriptify(string(l[2]))
    else
        ""
    end
end
p=graphplot!(ax, g_trail; 
           nlabels, layout=_->xy, arrow_size=0,
           nlabels_align=(:center,:bottom), nlabels_distance=6.0,
           node_color = RGBf(.1,.2,.6), edge_color = RGBf(.1,.2,.6),
           node_attr = (; strokewidth=3, strokecolor=:white),
           node_size = [m ? 15 : 0 for m in maximality]
           )

hidedecorations!(ax)  # hides ticks, grid and lables
hidespines!(ax)
f
# TODO: Debug visualization of SG 230

## ----------------------------------------------------------------------------------------
# Understanding how many permutations remain
let sgnum = 214
    sb, brs = compatibility_basis(sgnum; timereversal)
    lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirreps(sgnum))

    for nv in sb
        n = SymVector(nv, brs.irlabs, lgirsd)
        subgraphs, partitions_max, partitions_nonmax = build_subgraphs(n, subductionsd[sgnum], lgirsd)
        subgraphsp = permute_subgraphs(subgraphs)
        println(length.(subgraphsp) => prod(length.(subgraphsp)))
    end
end

## ----------------------------------------------------------------------------------------
## prune non-maximal vertices while retaining connectivity
# TODO: quite broken, e.g.:
#       - fails if any nonmax manifolds is not a line (e.g., plane) (e.g., `sgnum = 3`);
#         similar, but physically unrelated, failure when a k-line connects three irreps
#         where one of the three irreps is a "monodromy"-repetition
#       - doesn't correctly deal with two max-irreps of high-degen being connected by
#         multiple low-degen nonmax-irreps (e.g., `sgnum = 147; brs[end]`);
#       - if a max-manifold is not connected via a nonmax manifold, should it be pruned
#         (e.g., `sgnum = 6`) or connected via Ω? The latter seems preferable for
#         connectivity analysis
nonmaximal_vertices = partitions_nonmax[1].iridxs[1]:partitions_nonmax[end].iridxs[end]
es = collect(edges(g))
new_edges = Edge{Int}[]
for vidx in nonmaximal_vertices
    e12 = findall(e->e.dst==vidx, es)
    length(e12) != 2 && error("found more than two connected edges to non-maximal k-point")
    # FIXME: this fails when the connecting manifold is not a line (e.g., via F in SG 3)
    e1, e2 = es[e12[1]], es[e12[2]]
    e′ = Edge{Int}(e1.src, e2.src)
    push!(new_edges, e′)
end
g′ = Graph(partitions_max[end].iridxs[end])
foreach(e -> add_edge!(g′, e), new_edges)
# TODO: would be cool to restore information about the edge irrep as edge metadata

f, ax, p = graphplot(
    g′;
    nlabels=mapreduce(p->label.(p.lgirs), vcat, partitions_max),
    )

#deregister_interaction!(ax, :rectanglezoom)
#register_interaction!(ax, :ndrag, ndrag)
f