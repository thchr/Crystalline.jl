using Crystalline
if !isdefined(Main, :SymmetryBases)
    include("SymmetryBases.jl")
end
using Main.SymmetryBases

test_nontopo = true
spinful = false
timereversal = true
algorithm = "DualMode" # DualMode or PrimalMode
for sgnum in 1:MAX_SGNUM[3]
    BRS = bandreps(sgnum, spinful=spinful, timereversal=timereversal)
    
    B = matrix(BRS, true)        # Matrix with columns of EBRs.
    
    F   = Crystalline.smith(B)   # Smith normal decomposition of B
    dᵇˢ = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    Nⁱʳʳ, Nᴱᴮᴿ = size(B)

    # Print some simple stuff early on, to indicate that a calculation is running
    println("\nSG", sgnum, ": ", classification(BRS), 
            " (", dᵇˢ, " \"band structure dimensions\"; ", Nⁱʳʳ, " inequalities)")

    # Compatibility Hilbert bases  
    sb, zsᴴ = compatibility_bases(F, BRS, algorithm=algorithm)
    nsᴴ = matrix(sb) 
    Nᴴ = length(sb) # Number of Hilbert bases

    # Nontopological Hilbert bases 
    sb_nontopo, ysᴴ_nontopo = nontopological_bases(F, BRS, algorithm=algorithm)
    nsᴴ_nontopo = matrix(sb_nontopo)
    Nᴴ_nontopo  = length(sb_nontopo)

    # Splitting into trivial and fragile Hilbert bases
    nsᴴ_trivial, nsᴴ_fragile = split_fragiletrivial_bases(sb_nontopo, B)

    # Write some stats about the obtained Hilbert bases
    println("   ", Nᴱᴮᴿ,       " EBRs")
    println("   ", Nᴴ,         " compatibility bases")
    println("   ", Nᴴ_nontopo, " non-topological bases")
    println("      ", size(nsᴴ_trivial, 2), " trivial bases")
    println("      ", size(nsᴴ_fragile, 2), " fragile bases")

    # Test consistency of bases, if requested
    if test_nontopo
        SymmetryBases._test_hilbert_bases_consistency(BRS, F, nsᴴ, nsᴴ_nontopo, zsᴴ)
    end
end