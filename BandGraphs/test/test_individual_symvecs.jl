using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using Crystalline: irdim
using BandGraphs
using SymmetryBases
using GraphMakie
using GLMakie
const B = BandGraphs
##

D = 3
sgnum = 130
timereversal = true
criterion = (lgir) -> isspecial(lgir) && irdim(lgir) == 4
separable_degree = 2

lgirsd = timereversal ? realify!(lgirreps(sgnum, D)) : lgirreps(sgnum, D)
subts = subduction_tables(sgnum, D; timereversal)
sb, brs = compatibility_basis(sgnum, D; timereversal)

sb_ns = SymmetryVector.(sb, Ref(brs.irlabs), Ref(lgirsd))
ns = sb_ns[findall(n -> has_vertex(criterion, n), sb_ns)]
vs = filter(criterion, collect(Iterators.flatten(values(lgirsd))))

##

i = 30
n = ns[i]
bandg = build_subgraphs(n, subts, lgirsd)
#=
sep, bandg_splits = findall_separable_vertices(
    criterion, bandg; 
    max_permutations = 1e3,
    separable_degree = separable_degree,
    max_subblock_permutations = 1e3)
=#
##
faxp, _ = plot_flattened_bandgraph(bandg); faxp