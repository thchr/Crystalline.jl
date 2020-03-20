using SGOps, Test, LinearAlgebra

if !isdefined(Main, :LGIRS′)
    LGIRS = get_all_lgirreps(Val(3))  # loaded from our saved .jld2 files
end

@testset "CharacterTable and orthogonality theorems" begin
    @testset "Little group irreps" begin
        for lgirsvec in LGIRS
            for (idx, lgirs) in enumerate(lgirsvec)
                
                ct=SGOps.CharacterTable(lgirs)
                χs = characters(ct) # matrix of characters; each row is a different representation
                Nₒₚ = length(operations(ct))

                # 1ˢᵗ orthogonality theorem (characters):    ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾
                @test all(n->isapprox.(n, Nₒₚ, atol=1e-15), sum(abs2, χs, dims=2))

                # 2ⁿᵈ orthogonality theorem (characters):    ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦNₒₚ⁽ᵃ⁾ 
                @test conj(χs)*transpose(χs) ≈ Nₒₚ*I
            end
        end
    end
end