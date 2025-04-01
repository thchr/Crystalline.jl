using Crystalline, Test, LinearAlgebra, StaticArrays

@testset "Basis vectors" begin
    @testset "Conventional bases" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                Rs = directbasis(sgnum, D)
                Gs = dualbasis(Rs)
                
                @test all(dot(Gs[i], Rs[i])≈2π for i in 1:D)
                @test dualbasis(dualbasis(Rs)) ≈ Rs # involution property
            end
        end
    end
    @test dualbasis(([1,0.1], [0.1, 1])) == dualbasis([[1, 0.1], [0.1, 1]])

    @testset "Primitive bases" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                cntr = centering(sgnum, D)
                Rs   = directbasis(sgnum, D)
                Rs′  = primitivize(Rs, cntr)
                Gs   = dualbasis(Rs)
                Gs′ᵃ = primitivize(Gs, cntr)
                Gs′ᵇ = dualbasis(Rs′)

                @test Gs′ᵃ ≈ Gs′ᵇ
                @test all(dot(Gs′ᵇ[i], Rs′[i]) ≈ 2π for i in 1:D)
                @test conventionalize(Rs′, cntr) ≈ Rs
                @test conventionalize(Gs′ᵃ, cntr) ≈ Gs
            end
        end
    end

    @testset "crystal(...)" begin
        # zero-elements should be _exactly_ =0 for π/2 angles; not merely approximately 0
        Rs = crystal(1.0, 1.0, 1.0, π/2, π/2, π/2)
        @test all(iszero, (map(filter(R->abs(R)<0.1), Rs)))
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
        k = ReciprocalPoint(.1,.2,.3)
        r = DirectPoint(.1,.2,.3)

        # conversion
        @test k == convert(ReciprocalPoint{3}, v)  # AbstractVector conversion
        @test k == convert(ReciprocalPoint{3}, s)  # StaticVector conversion

        # construction
        @test k == ReciprocalPoint(s) == ReciprocalPoint{3}(s)
        @test k == ReciprocalPoint(t) == ReciprocalPoint{3}(t)
        @test k == ReciprocalPoint(v) == ReciprocalPoint{3}(v)

        # transformation
        Rs = directbasis(2)
        Gs = dualbasis(Rs)
        @test cartesianize(k, Gs) isa ReciprocalPoint{3}
        @test latticize(cartesianize(k, Gs), Gs) isa ReciprocalPoint{3}
        @test cartesianize(r, Rs) isa DirectPoint{3}
        @test latticize(cartesianize(r, Rs), Rs) isa DirectPoint{3}

        @test cartesianize(parent(k), Gs) isa typeof(parent(k))
        @test latticize(cartesianize(parent(k), Gs), Gs) isa typeof(parent(k))

        # arithmetic
        @test k + k isa ReciprocalPoint{3} && k + k == s + s
        @test k - k isa ReciprocalPoint{3} && k - k == zero(k)
        @test zero(k) isa ReciprocalPoint{3}
        @test -k isa ReciprocalPoint{3}
        @test 2.3k isa ReciprocalPoint{3} && 2.3k == 2.3s
        @test 2.3r - r + 3r isa DirectPoint{3} && 2.3r - r + 3r ≈ 4.3r ≈ 4.3s

        # isapprox on near-zero-difference
        a = 4.3r
        b = 2.3r - r + 3r
        @test isapprox(a, b)
        @test isapprox(RVec(a), RVec(b))
    end
end

using Unitful: ustrip, @u_str
@testset "Unitful `AbstractBasis`" begin
    Rs = DirectBasis([1.0, 0.0, 0.0], [-0.5, sqrt(3)/2, 0.0],   [0, 0, 1.25])
    Us = DirectBasis{3}(map(Rs) do v v * 1u"Å" end)

    # issue #56
    @test all(zip(dualbasis(Rs), dualbasis(Us))) do (R, U)
        R == ustrip(U)
    end
    @test all(zip(dualbasis(Rs.vs), dualbasis(Us.vs))) do (R, U)
        R == ustrip(U)
    end
end