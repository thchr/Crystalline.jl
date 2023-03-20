using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie

## ----------------------------------------------------------------------------------------
## Testing & visualization of band graph

timereversal = true
sgnum = 200
sb, brs = compatibility_basis(sgnum; timereversal)
lgirsd = lgirreps(sgnum)
timereversal && (lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd))
_n = brs[end]
n = SymVector(_n, brs.irlabs, lgirsd)

subgraphs, partitions = build_subgraphs(n, SUBDUCTIONSD[sgnum], lgirsd)

A = assemble_adjacency(subgraphs, partitions)
L = assemble_laplacian(subgraphs, partitions)
g = assemble_graph(subgraphs, partitions) # structured equiv of `Graph(A)`

node_colors = [g[label_for(g, i)].maximal ? :red : :black for i in vertices(g)]
f, ax, p = graphplot(
    g;
    nlabels=[label(g[label_for(g, i)].lgir) for i in vertices(g)],
    nlabels_distance=4,
    node_color = node_colors, nlabels_color=node_colors,
    node_attr = (; strokewidth=4, strokecolor=:white),
    edge_color = :gray,
    layout = GraphMakie.NetworkLayout.SFDP(; C=10, K=.2, dim=2)
    )

make_vertices_dragable!(ax)
f

## ----------------------------------------------------------------------------------------
# Testing & visualization of k-connectivity graphs

kg = partition_graph(subgraphs, partitions)
kg′ = split_nonmaximal_nodes(kg)
f,ax,p=graphplot(kg′;
    nlabels = [label_for(kg′, i)[1] for i in vertices(kg′)],
    nlabels_distance=4,
    node_color = [kg′[label_for(kg′, i)].maximal ? :red : :black for i in vertices(kg′)],
    nlabels_color = [kg′[label_for(kg′, i)].maximal ? :red : :black for i in vertices(kg′)], 
    node_attr = (; strokewidth=5, strokecolor=:white),
    )

kg_trail = chinese_postman(kg′)
dkg = DiGraph(kg_trail)
graphplot!(ax, dkg,
    layout=_->p[:node_pos][],
    edge_color=RGB(0.2,0.2,.9),
    edge_width = 2,
    node_size = 0,
    curve_distance_usage = true)
f

## ----------------------------------------------------------------------------------------


sgnum = 200
sb, brs = compatibility_basis(sgnum; timereversal)
lgirsd = lgirreps(sgnum)
timereversal && (lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd))
_n = brs[end-1]
n = SymVector(_n, brs.irlabs, lgirsd)

f1, ax, plot_data = plot_flattened_bandgraph(n, lgirsd)
f1
# TODO: Debug visualization of SG 230

## ----------------------------------------------------------------------------------------
# Create a permuted graph (a single subgraph permutation) and compare
subgraphs, partitions = build_subgraphs(n, SUBDUCTIONSD[sgnum], lgirsd)
psubgraphs = permute_subgraphs(subgraphs)
p_idx = 12 # pick the second possible permutation of `pidx` subgraph
psubgraph_pidx = psubgraphs[p_idx] # a 
permuted_subgraph = SubGraph(
    psubgraph_pidx.p_max, psubgraph_pidx.p_nonmax, psubgraph_pidx.As[2], false)
permuted_subgraphs = vcat(
    subgraphs[1:p_idx-1], 
    permuted_subgraph,
    subgraphs[p_idx+1:end])

f2, _ = plot_flattened_bandgraph(subgraphs, partitions)

display(GLMakie.Screen(), f1)
display(GLMakie.Screen(), f2) # open both
## ----------------------------------------------------------------------------------------
# Understanding how many permutations remain

let sgnum = 73
    sb, brs = compatibility_basis(sgnum; timereversal)
    lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirreps(sgnum))

    for nv in sb
        n = SymVector(nv, brs.irlabs, lgirsd)
        subgraphs, partitions = build_subgraphs(n, SUBDUCTIONSD[sgnum], lgirsd)
        subgraphsp = permute_subgraphs(subgraphs)
        subgraph_klabs = [" (" * s.p_max.klab * "→" * s.p_nonmax.klab .* ")" for s in subgraphs]
        join(stdout, string.(length.(subgraphsp)) .* subgraph_klabs, ", ")
        println(" ⇒ ", prod(length.(subgraphsp)))
    end
end

## ----------------------------------------------------------------------------------------
## TODO: prune non-maximal vertices while retaining connectivity
# tricky bits include:
#       - fails if any nonmax manifolds is not a line (e.g., plane) (e.g., `sgnum = 3`);
#         similar, but physically unrelated, failure when a k-line connects three irreps
#         where one of the three irreps is a "monodromy"-repetition.
#         can be fixed by using on output of `unfold_bandgraph` or 
#       - doesn't correctly deal with two max-irreps of high-degen being connected by
#         multiple low-degen nonmax-irreps (e.g., `sgnum = 147; brs[end]`);
#       - if a max-manifold is not connected via a nonmax manifold, should it be pruned
#         (e.g., `sgnum = 6`) or connected via Ω? The latter seems preferable for
#         connectivity analysis