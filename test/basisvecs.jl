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
    @test reciprocalbasis(([1,0.1], [0.1, 1])) == reciprocalbasis([[1, 0.1], [0.1, 1]])

    @testset "Primitive bases" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                cntr = centering(sgnum, D)
                Rs   = directbasis(sgnum, D)
                Rs′  = primitivize(Rs, cntr)
                Gs   = reciprocalbasis(Rs)
                Gs′ᵃ = primitivize(Gs, cntr)
                Gs′ᵇ = reciprocalbasis(Rs′)

                @test Gs′ᵃ ≈ Gs′ᵇ
                @test all(dot(Gs′ᵇ[i], Rs′[i]) ≈ 2π for i in 1:D)
                @test conventionalize(Rs′, cntr) ≈ Rs
                @test conventionalize(Gs′ᵃ, cntr) ≈ Gs
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

    @testset "AbstractPoint" begin
        t = (.1,.2,.3)
        v = [.1,.2,.3]
        s = @SVector [.1,.2,.3]
        r = ReciprocalPoint(.1,.2,.3)

        # conversion
        @test r == convert(ReciprocalPoint{3}, v)  # AbstractVector conversion
        @test r == convert(ReciprocalPoint{3}, s)  # StaticVector conversion

        # construction
        @test r == ReciprocalPoint(s) == ReciprocalPoint{3}(s)
        @test r == ReciprocalPoint(t) == ReciprocalPoint{3}(t)
        @test r == ReciprocalPoint(v) == ReciprocalPoint{3}(v)
    end
end