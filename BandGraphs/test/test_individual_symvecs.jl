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
sgnum = 131
timereversal = true
criterion = (lgir) -> isspecial(lgir) && irdim(lgir) == 2
separable_degree = nothing

lgirsd = timereversal ? realify!(lgirreps(sgnum, D)) : lgirreps(sgnum, D)
subts = subduction_tables(sgnum, D; timereversal)
sb, brs = compatibility_basis(sgnum, D; timereversal)

sb_ns = SymmetryVector.(sb, Ref(brs.irlabs), Ref(lgirsd))
ns_vertex = sb_ns[findall(n -> has_vertex(criterion, n), sb_ns)]
vs = filter(criterion, collect(Iterators.flatten(values(lgirsd))))

##

n = infeasible[131][end][2]
bandg = build_subgraphs(n, subts, lgirsd)
subgraphs_ps = permute_subgraphs(bandg.subgraphs; with_weyl_filter=true)
bandgp = B.BandGraphPermutations(bandg.partitions, subgraphs_ps);
findall_separable_vertices(
    criterion, n, subts, lgirsd; 
    max_permutations = 1e6,
    separable_degree = separable_degree,
    max_subblock_permutations = 1e4,
    with_weyl_filter=false)
length(bandgp)

subgraphs_ps_w = permute_subgraphs(bandg.subgraphs; with_weyl_filter=true)
bandgp_w = B.BandGraphPermutations(bandg.partitions, subgraphs_ps_w);
length(bandgp_w)

#nothing
#=
sep, bandg_splits = findall_separable_vertices(
    criterion, bandg; 
    max_permutations = 1e3,
    separable_degree = separable_degree,
    max_subblock_permutations = 1e4)
=#
##
#faxp, _ = plot_flattened_bandgraph(bandg); faxp