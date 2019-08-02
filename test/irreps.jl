using SGOps, Test
if true
IR = parseisoir(Complex);
@testset "1ˢᵗ orthogonality theorem" begin
    for sgnum in 1:230          # test for every space group
        @testset "SG$sgnum" begin 
            firstorthogacc = true
            for ir in IR[sgnum] # test for every CDML label, i.e. (k-points)∪(irreps)
                lgir    = littlegroupirrep(ir)
                lgops   = operations(lgir)
                chars   = characters(lgir)

                # 1ˢᵗ orthogonality theorem requires that ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾ for each 
                # irrep (a) with i running over the Nₒₚ elements of the little group 
                firstorthogacc &= sum(abs2, chars) ≈ length(lgops) 
            end
            @test firstorthogacc 
        end
    end
end
end

LGIR = parselittlegroupirreps.(IR)
if true
@testset "2ⁿᵈ orthogonality theorem" begin
    # test ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾* = δₐᵦNₒₚ⁽ᵃ⁾  for distinct little group irreps (a) ≠ (β)
    for lgirs in LGIR        # lgirs: vectors of little group irrep collections
        for lgir in lgirs    # lgir:  tuples of distinct little group irreps
            Nₒₚ = order(first(lgir))
            @test all(x->order(x)==Nₒₚ, lgir)     # test that the size of the little group is identical at fixed 𝐤      
            for i in eachindex(lgir) 
                charsᵢ = characters(lgir[i])
                for j in eachindex(lgir)
                    orthog2nd = charsᵢ'*characters(lgir[j])
                    @test isapprox(orthog2nd, (i==j)*Nₒₚ, atol = 1e-12)
                end
            end
        end
    end
end
end

