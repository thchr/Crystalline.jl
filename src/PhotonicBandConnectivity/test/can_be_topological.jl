using Crystalline, LinearAlgebra, Test, JuMP, GLPK
import Crystalline: rotation, rotation_order_3d
if !isdefined(Main, :(PhotonicBandConnectivity))
    includet("../src/PhotonicBandConnectivity.jl")
    using Main.PhotonicBandConnectivity
end
if !isdefined(Main, :(SymmetryBases))
    includet("../../SymmetryBases/SymmetryBases.jl")
    using Main.SymmetryBases
end

rotation_order(op::SymOperation{3}) = (W=rotation(op); rotation_order_3d(det(W), tr(W)))

sgnums = collect(1:MAX_SGNUM[3])
# restrict to sgs that do not have roto-inversions
filter!(sgnums) do sgnum
    !any(op->rotation_order(op)∈(-1, -3, -4, -6), spacegroup(sgnum, Val(3)))
end
# restrict to sgs that have regular 2T representations
filter!(sgnums) do sgnum
    all(≥(0), PhotonicBandConnectivity.find_representation²ᵀ(sgnum, timereversal=true))
end
# of these (68) space groups, just 10 have nontrivial topological classification (all Z₂):
#filter!(sgnums) do sgnum
#    classification(bandreps(sgnum, timereversal=true))≠"Z₁"
#end

# determine if expansions of 2T must include any nontrivial Hilbert bases via symmetry bases
nontopo_data = nontopological_bases.(sgnums, timereversal=true)
nontopo_sbs  = getindex.(nontopo_data, 1)

data = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=true)
cⁱss = getindex.(data, 1)
νᵀs  = getindex.(data, 2)
sbs  = getindex.(data, 3)

# determine the trivial respectively fragile indices in nontopo_sbs
BRSs       = bandreps.(sgnums, timereversal=true)
Bs         = matrix.(BRSs, true)
fragile_data = split_fragiletrivial_bases.(nontopo_sbs, Bs)
trivial_idxs = getindex.(fragile_data, 1)
fragile_idxs = getindex.(fragile_data, 2)

# compute the actual symmetry vectors for each compatibility-constraint band solution
nsᵀs = [unique!(sort(sum_symbases.(Ref(sb), cⁱs))) for (sb, cⁱs) in zip(sbs,  cⁱss)]

# ---------------------------------------------------------------------------------------- #

# for each of the band-solution in `nsᵀs`, check whether or not it can be expanded solely 
# using the nontopological basis (then it is a trivial band-combination!); if not, it must
# be a nontrivial and topological band-combination
println("\n-------------------\n")
for (sgidx, sgnum) in enumerate(sgnums)
    print("SG ", sgnum, ", ν = ", νᵀs[sgidx])

    # prep-work to get matrices; to avoid creating them multiple times in the loop
    M         = matrix(sbs[sgidx])         # all bases
    nontopo_M = matrix(nontopo_sbs[sgidx]) # nontopological bases only
    can_be_fragile = !isempty(fragile_idxs[sgidx])
    trivial_M = can_be_fragile ? (@views nontopo_M[:,trivial_idxs[sgidx]]) : nothing # trivial bases only (excl. fragile bases)

    # some printing about global properties of the space group wrt. topology
    classification(BRSs[sgidx]) ≠ "Z₁" && print(" [+NONTRIVIAL]")
    !isempty(fragile_idxs[sgidx])      && print(" [+FRAGILE]")
    println()

    # compute & print topology of each solution
    for nᵀ in nsᵀs[sgidx]
        topology_kind = get_solution_topology(nᵀ, nontopo_M, trivial_M, M)

        print("   ", topology_kind, " ⇐  ")
        Crystalline.prettyprint_symmetryvector(stdout, nᵀ, irreplabels(BRSs[sgidx]))
        println()
        # TODO: Maybe store `topology_kind` in some vector for each space group?
    end
    println()
end
