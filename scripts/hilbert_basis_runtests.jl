using SGOps
if !isdefined(Main, :SGOpsHilbertBases)
    include("hilbert_basis.jl")
end
using Main.SGOpsHilbertBases

test_nontopo = true
spinful = false
timereversal = true
algorithm = "DualMode" # DualMode or PrimalMode
for sgnum in 1:230
    BRS = bandreps(sgnum, spinful=spinful, timereversal=timereversal)
    
    B = collect(transpose(matrix(BRS, true))) # Matrix with columns of EBRs.
    
    F   = SGOps._smith′(B)      # Smith normal decomposition of B
    dᵇˢ = count(!iszero, F.SNF) # "Dimensionality" of band structure
    Nⁱʳʳ, Nᴱᴮᴿ = size(B)

    # Print some simple stuff early on, to indicate that a calculation is running
    println("\nSG", sgnum, ": ", classification(BRS), 
            " (", dᵇˢ, " \"band structure dimensions\"; ", Nⁱʳʳ, " inequalities)")

    # Compatibility Hilbert bases  
    nsᴴ, zsᴴ = compatibility_bases(F, algorithm=algorithm) 
    Nᴴ = size(nsᴴ, 2) # Number of Hilbert bases

    # Nontopological Hilbert bases 
    nsᴴ_nontopo, ysᴴ_nontopo = nontopological_bases(F, algorithm=algorithm)
    Nᴴ_nontopo  = size(nsᴴ_nontopo, 2)

    # Splitting into trivial and fragile Hilbert bases
    nsᴴ_trivial, nsᴴ_fragile = split_fragiletrivial_bases(nsᴴ_nontopo, B)

    # Write some stats about the obtained Hilbert bases
    println("   ", Nᴱᴮᴿ,       " EBRs")
    println("   ", Nᴴ,         " compatibility bases")
    println("   ", Nᴴ_nontopo, " non-topological bases")
    println("      ", size(nsᴴ_trivial, 2), " trivial bases")
    println("      ", size(nsᴴ_fragile, 2), " fragile bases")

    # Test consistency of bases, if requested
    if test_nontopo
        SGOpsHilbertBases._test_hilbert_bases_consistency(BRS, F, nsᴴ, nsᴴ_nontopo, zsᴴ)
    end
end