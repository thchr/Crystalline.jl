using Crystalline, Test, LinearAlgebra, StaticArrays

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

    @testset "Mixed inputs for constructors" begin
        Rs = ([1.0, 0.0, 0.0], [-0.5, sqrt(3)/2, 0.0],   [0, 0, 1.25])
        Rs′ = DirectBasis{3}(Rs)
        @test Rs′ == DirectBasis(Rs)                         # NTuple{3, Vector{3, Float64}}
        Rsˢˢ = convert(SVector{3,SVector{3,Float64}}, Rs)  # SVector{3, SVector{3, Float64}}
        @test Rs′ == DirectBasis{3}(Rsˢˢ) == DirectBasis(Rsˢˢ)
        Rsᵛ = [Rs...]                                              # Vector{Vector{Float64}}
        @test Rs′ == DirectBasis{3}(Rsᵛ) == DirectBasis(Rsᵛ)
        Rsˢ = SVector{3,Float64}.(Rs)                       # NTuple{3, SVector{3, Float64}}
        @test Rs′ == DirectBasis{3}(Rsˢ) == DirectBasis(Rsˢ)
        Rsᵗ = ntuple(i->ntuple(j->Rs[i][j], Val(3)), Val(3)) # NTuple{3, NTuple{3, Float64}}
        @test Rs′ == DirectBasis{3}(Rsᵗ) == DirectBasis(Rsᵗ)
        @test Rs′ == DirectBasis(Rs[1], Rs[2], Rs[3]) == DirectBasis{3}(Rs[1], Rs[2], Rs[3])
                                                  # ↑ Varargs / splatting of Vector{Float64}
        @test Rs′ == DirectBasis(Rsˢ[1], Rsˢ[2], Rsˢ[3]) == DirectBasis{3}(Rsˢ[1], Rsˢ[2], Rsˢ[3])
                                              # ↑ Varargs / splatting of SVector{3, Float64}
    end
end