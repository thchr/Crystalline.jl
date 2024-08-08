using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie

## --------------------------------------------------------------------------------------- #
include(joinpath((@__DIR__), "..", "src", "complete-split.jl"))

## --------------------------------------------------------------------------------------- #

D = 2
timereversal = true
sgnum = 17
sb, brs = compatibility_basis(sgnum, D; timereversal)
lgirsd = lgirreps(sgnum, D)
timereversal && realify!(lgirsd)
subts = subduction_tables(sgnum, D; timereversal)
_n = sb[end]
n = SymVector(_n, brs.irlabs, lgirsd)

bandg = build_subgraphs(n, subts, lgirsd)

## --------------------------------------------------------------------------------------- #

v = ("K₃", 1)
bandgs′ = something(complete_split(bandg, v));

## --------------------------------------------------------------------------------------- #

g = assemble_graph(bandgs′[1])

node_colors = [g[label_for(g, i)].maximal ? :red : :black for i in vertices(g)]
f, ax, p = graphplot(
    g;
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
# Create a permuted graph (a single subgraph permutation) and compare
subgraphs_ps = permute_subgraphs(bandg.subgraphs);
bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps);
BandGraphs.permutation_info(bandgp)
length(bandgp) > 10 && @info("Be warned: many permutations ($(length(bandgp)))")

v = ("K₃", 1)
bandgs′ = something(complete_split(bandgp[1], v));
bandg′ = bandgs′[2]

faxp, (; xs, ys) = plot_flattened_bandgraph(bandg′)
display(GLMakie.Screen(), faxp)
