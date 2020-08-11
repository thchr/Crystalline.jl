using Test
using SymmetryBases
using Crystalline
using SmithNormalForm

"""
    _test_hilbert_bases_consistency(BRS::BandRepSet, F::SmithNormalForm.Smith,
                    nsᴴ::AbstractMatrix, nsᴴ_nontopo::AbstractMatrix, zsᴴ::AbstractMatrix)

Test that the obtained "non-topological" bases indeed obey some of the conditions that we
know they must. Prints a checkmark (✓) if succesful; throws `Test.FallbackTestSetException`
otherwise. Returns `nothing`.
"""
function _test_hilbert_bases_consistency(BRS::BandRepSet, F::SmithNormalForm.Smith,
                nsᴴ::AbstractMatrix, nsᴴ_nontopo::AbstractMatrix, zsᴴ::AbstractMatrix)

    print("   ... checking consistency of non-topological bases: ")

    dᵇˢ        = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    Nᴴ         = size(nsᴴ, 2)           # Number of compatibility Hilbert bases
    Nᴴ_nontopo = size(nsᴴ_nontopo, 2)   # Number of nontopological Hilbert bases

    # All zsᴴ that are not element divisible by Λ ≡ F.SNF correspond to "proper" topological
    # cases (i.e. not fragile). This is because those z ≡ Λy with yᵢ∈ℚ\ℤ are topological, 
    # whereas all those with yᵢ∈ℤ are either trivial or fragilely topological.
    # Note: The approach below is not valid in general: while it does find all the 
    #       non-topological elements among the Hilbert bases nsᴴ, it does not find a full
    #       Hilbert basis for all non-topological states. The easiest way to see this is to
    #       think in terms of a "unit cell" for the Hilbert bases, and then realize that 
    #       this unit cell may be far larger, when we take a subcone of the original cone.
    Λ = @view F.SNF[OneTo(dᵇˢ)] # Nonzero invariant factors of Smith normal decomposition
    nontopo_idxs_subset = findall(zᴴ -> all(zᴴΛᵢ -> mod(zᴴΛᵢ[1], zᴴΛᵢ[2]) == 0, zip(zᴴ, Λ)),
                           collect(eachcol(zsᴴ)))
    topo_idxs_subset   = findall(i -> i ∉ nontopo_idxs_subset, OneTo(Nᴴ))
    nsᴴ_nontopo_subset = @view nsᴴ[:, nontopo_idxs_subset]
    nsᴴ_topo_subset    = @view nsᴴ[:, topo_idxs_subset]

    # If classification is Z₁, nsᴴ and nsᴴ_nontopo must be equal
    if classification(BRS) == "Z₁"
        @test Set(eachcol(nsᴴ)) == Set(eachcol(nsᴴ_nontopo))
        println("✓ (trivially)")
    else 
        # Must be a superset of the naive extraction approach
        @test Set(eachcol(nsᴴ_nontopo_subset)) ⊆ Set(unique(eachcol(nsᴴ_nontopo)))

        # Check that every basis vector in nsᴴ_nontopo can be expanded in the compatibility 
        # basis nsᴴ using only non-negative integer coefficients
        if Nᴴ < Nᴴ_nontopo # no need to test unless there are more elements in nsᴴ_nontopo
            for (j,nsᴴ_nontopoʲ) in enumerate(eachcol(nsᴴ_nontopo))
                j ≠ 1 && print(stdout, "\b"^ndigits(j-1))
                print(stdout, j)
                flush(stdout)

                m = Model(GLPK.Optimizer)
                @variable(m, c[1:Nᴴ] >= 0, Int)
                @constraint(m, nsᴴ*c .== nsᴴ_nontopoʲ)
                optimize!(m)

                cvals = Int.(value.(c)) # extract optimized expansion coefficients
                @test nsᴴ*cvals == nsᴴ_nontopoʲ
                @test all(cvals .≥ 0)
            end
        end
        print("\b"^ndigits(Nᴴ_nontopo), "✓", " "^ndigits(Nᴴ_nontopo-1))
    end

    nothing
end