using Crystalline, Test, LinearAlgebra

@testset "Basis vectors" begin
    @testset "Conventional bases" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                Rs = directbasis(sgnum, D)
                Gs = reciprocalbasis(Rs)
                
                @test all(dot(Gs[i], Rs[i])≈2π for i in 1:D)
            end
        end
    end

    @testset "Primitive bases" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                cntr = centering(sgnum, D)
                Rs = directbasis(sgnum, D)
                Rs′ = primitivize(Rs, cntr)
                Gs′ᵃ = primitivize(reciprocalbasis(Rs), cntr)
                Gs′ᵇ = reciprocalbasis(Rs′)

                @test Gs′ᵃ ≈ Gs′ᵇ
                @test all(dot(Gs′ᵇ[i], Rs′[i])≈2π for i in 1:D)
            end
        end
    end
end