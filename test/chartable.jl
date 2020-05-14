using Crystalline, Test, LinearAlgebra

if !isdefined(Main, :LGIRS′)
    LGIRS = get_all_lgirreps(Val(3))  # loaded from our saved .jld2 files
end

@testset "CharacterTable and orthogonality theorems" begin
    @testset "Little group irreps" begin
        for lgirsvec in LGIRS
            for (idx, lgirs) in enumerate(lgirsvec)
                
                ct= CharacterTable(lgirs)
                χs = characters(ct) # matrix of characters; each row is a different representation
                Nₒₚ = length(operations(ct))

                # 1ˢᵗ orthogonality theorem:    ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾
                @test all(n->isapprox.(n, Nₒₚ, atol=1e-14), sum(abs2, χs, dims=1))

                # 2ⁿᵈ orthogonality theorem:    ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦNₒₚ⁽ᵃ⁾ 
                @test χs'*χs ≈ Nₒₚ*I
            end
        end
    end

    @testset "Point group irreps" begin
        for D in 1:3
            for pgiuc in PGS_IUCs[D]
                pgirs = get_pgirreps(pgiuc, Val(D))
                ct = CharacterTable(pgirs)
                χs = characters(ct) # matrix of characters; each row is a different representation
                Nₒₚ = length(operations(ct))

                # 1ˢᵗ orthogonality theorem:    ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾
                @test all(n->isapprox.(n, Nₒₚ, atol=1e-14), sum(abs2, χs, dims=1))

                # 2ⁿᵈ orthogonality theorem:    ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦNₒₚ⁽ᵃ⁾ 
                @test χs'*χs ≈ Nₒₚ*I
            end
        end
    end
end