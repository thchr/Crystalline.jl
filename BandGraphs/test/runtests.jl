using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie

## ----------------------------------------------------------------------------------------
## Testing & visualization of band graph

D = 2
timereversal = true
sgnum = 17
sb, brs = compatibility_basis(sgnum, D; timereversal)
lgirsd = lgirreps(sgnum, D)
timereversal && realify!(lgirsd)
subts = subduction_tables(sgnum, D; timereversal)
_n = sb[end]
#_n = brs[end-1]
n = SymmetryVector(_n, brs.irlabs, lgirsd)

bandg = build_subgraphs(n, subts, lgirsd)
subgraphs, partitions = bandg.subgraphs, bandg.partitions

A = assemble_adjacency(bandg)
L = assemble_laplacian(bandg)
g = assemble_graph(bandg) # structured equiv of `Graph(A)`

node_colors = [g[label_for(g, i)].maximal ? :red : :black for i in vertices(g)]
f, ax, p = graphplot(
    g;
    #nlabels=[(l = label_for(g, i); l[1] * Crystalline.supscriptify(string(l[2]))) for i in vertices(g)],
    nlabels=[(l = label_for(g, i); l[1]) for i in vertices(g)],
    nlabels_distance=4,
    node_color = node_colors, nlabels_color=node_colors,
    node_attr = (; strokewidth=4, strokecolor=:white),
    edge_color = :gray,
    layout = GraphMakie.NetworkLayout.SFDP(; C=10, K=.2, dim=2)
    )

make_vertices_dragable!(ax, p)
f

## ----------------------------------------------------------------------------------------
# Testing & visualization of k-connectivity graphs

kg = partition_graph(bandg)
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
    edge_color = RGBf(0.2,0.2,.9),
    edge_width = 2,
    node_size = 0,
    curve_distance_usage = true)
f

## ----------------------------------------------------------------------------------------
# Create a permuted graph (a single subgraph permutation) and compare
subgraphs_ps = permute_subgraphs(bandg.subgraphs);
bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps);
BandGraphs.permutation_info(bandgp)
length(bandgp) > 10 && @info("Be warned: many permutations ($(length(bandgp)))")

GLMakie.closeall()
xys = nothing
maxplot = 4
fs = Vector{Figure}(undef, min(length(bandgp), maxplot))
for (n, bandg′) in enumerate(bandgp)
    n > maxplot && break
    faxp, (; xs, ys) = plot_flattened_bandgraph(bandg′; xys=xys)
    fs[n], ax, p = faxp
    ax.title = "Permutation $(n)"
    display(GLMakie.Screen(), fs[n])
    global xys = (xs, ys)
end


## ----------------------------------------------------------------------------------------
# Understanding how many permutations remain

let sgnum = 96, D = 3 # 200
    sb, brs = compatibility_basis(sgnum; timereversal)
    #sb, brs = nontopological_basis(sgnum; timereversal)
    lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirreps(sgnum))
    subts = subduction_tables(sgnum, D; timereversal)

    for nv in sb
        n = SymmetryVector(nv, brs.irlabs, lgirsd);
        bandg = build_subgraphs(n, subts, lgirsd)
        subgraphs_ps = permute_subgraphs(bandg.subgraphs)
        bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps);
        println(n)
        BandGraphs.permutation_info(bandgp)
        println()
    end
end