if !isdefined(Main, :(PhotonicBandConnectivity))
    includet("../src/PhotonicBandConnectivity.jl")
    using Main.PhotonicBandConnectivity
end
if !isdefined(Main, :(SymmetryBases))
    includet("../../SymmetryBases/SymmetryBases.jl")
    using Main.SymmetryBases
end

using Test
#= 
This is a test script to check the topology of the 2T solutions in the cases where the 2T 
solutions is irregular, such that it requires inclusion of a 1L mode. The approach is to 
contrast the topology of the 1L mode choice/pick with the topology of the 2T+1L solution.
=#

sgnums = 1:230
has_tr = true

# ω=0 solutions
data = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=has_tr);
cⁱss   = getindex.(data, 1) # coefficients of expansions
νᵀs    = getindex.(data, 2) # fillings for tranverse branch
sbs    = getindex.(data, 3) # symmetry bases
idx¹ᴸs = getindex.(data, 4) # index for chosen 1L branch

# nontopological Hilbert bases
nontopo_sbs = getindex.(nontopological_bases.(sgnums, timereversal=has_tr), 1)

# check whether 1L pick is trivial
for (sgidx, sgnum) in enumerate(sgnums)
    idx¹ᴸ = idx¹ᴸs[sgidx]
    
    if idx¹ᴸ !== nothing #&& classification(bandreps(sgnum, timereversal=has_tr)) ≠ "Z₁"
        # TODO: Treat the simpler `idx¹ᴸ !== nothing` case as well within this script?
        sb, nontopo_sb = sbs[sgidx], nontopo_sbs[sgidx]
        nᴸ = sb[idx¹ᴸ]
        # check topology of 1L pick
        is_trivialᴸ = nᴸ ∈ nontopo_sb

        println("SG ", sgnum)

        # now check topology of 2T+1L solutions: the approach is the same as that in 
        # src/topology_as_2T1L_vs_1L_difference.jl (see `topology_from_2T1L_xor_1L`).
        ns = unique!(sort(sum_symbases.(Ref(sb), cⁱss[sgidx])))
        nontopo_M = matrix(nontopo_sb)
        println("   ", length(ns), " ω=0 solutions")
        # TODO: use `topology_from_2T1L_xor_1L` instead of this?
        trivial_count²ᵀ⁺¹ᴸ = nontrivial_count²ᵀ⁺¹ᴸ = 0
        for n in ns
            topology_kind = get_solution_topology(n, nontopo_M, nothing, nothing)
            if topology_kind == trivial
                trivial_count²ᵀ⁺¹ᴸ += 1
            elseif topology_kind == nontrivial
                nontrivial_count²ᵀ⁺¹ᴸ += 1
            end
            #Crystalline.prettyprint_symmetryvector(stdout, n - nᴸ, irreplabels(bandreps(230, timereversal=has_tr)))
            #println()
        end
        
        # print summary of stats
        println("   │  1L pick:         ", is_trivialᴸ ? "trivial" : "nontrivial")
        println("   │  2T+1L solutions: ", trivial_count²ᵀ⁺¹ᴸ, " trivial", " / ", 
                                           nontrivial_count²ᵀ⁺¹ᴸ, " nontrivial")
        println("   ╰─ 2T solutions:    ", 
                is_trivialᴸ ? trivial_count²ᵀ⁺¹ᴸ : nontrivial_count²ᵀ⁺¹ᴸ, " trivial", " / ",
                is_trivialᴸ ? nontrivial_count²ᵀ⁺¹ᴸ : trivial_count²ᵀ⁺¹ᴸ, " nontrivial")
        println()

        # NOTE: SGs 48, 50, 68, 86 are interesting (all 2T solutions are nontrivial)
    end
end

# TODO: we should check that we get the same 2T Z₂ index by seing whether or not we can 
#       expand the non-Γ parts of the 2T solutions in a strictly positive coefficient EBR 
#       expansion