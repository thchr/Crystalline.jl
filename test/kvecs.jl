using Test, Crystalline, StaticArrays, LinearAlgebra

@testset "KVec" begin

    @testset "Construction, parsing, and show" begin
        # parsing: primary functionality
        @test dim(KVec("u")) == 1
        @test dim(KVec("0,0.5")) == 2
        @test dim(KVec("0.3,0.2,.1")) == 3
        @test parts(KVec("u,v,w")) == ([0,0,0], [1 0 0; 0 1 0; 0 0 1])

        @test KVec{1}("-u") == KVec{1}("-α") == KVec(zeros(1), fill(-1,1,1)) == KVec(zero(SVector{1,Float64}), -one(SMatrix{1,1,Float64}))

        @test string(KVec{3}("u+1,0,0.5")) == string(KVec{3}("1.0+u,0,0.5")) == "[1+α, 0, 1/2]"
        @test string(KVec{3}("u+v,u,0"))   == string(KVec{3}("α+β,α,0"))     == "[α+β, α, 0]"
        @test string(KVec{3}("-0.3,0.2,-.1")) == "[-3/10, 1/5, -1/10]"

        @test string(KVec{3}("α, 2α, 0.0"))   == string(KVec{3}("α, 2.0α, 0.0")) == "[α, 2α, 0]"
        @test string(KVec{3}("α, 2α+2, 0.0")) == string(KVec{3}("α, 2.0α+2.0, 0.0")) == "[α, 2+2α, 0]"
        @test string(KVec{3}("β, -2/3, -0.75")) == string(KVec{3}("v, -0.66666666666666667, -3/4")) == "[β, -2/3, -3/4]"

        # parsing: allow some basic multiplication usage
        @test sprint(show, KVec("2+3*u, 1-2*v")) == "[2+3α, 1-2β]"
        @test sprint(show, KVec("-3*u, -v")) == "[-3α, -β]"
        @test sprint(show, KVec("3*u, v")) == "[3α, β]"
        @test sprint(show, KVec("-3*u+1, -β+1/2")) == "[1-3α, 1/2-β]"
        @test sprint(show, KVec("-1*u+1/3, 2+1*β")) == "[1/3-α, 2+β]"

        # functor-like usage
        @test KVec{3}("x,y,z")(0,1,2) == [0.0,1.0,2.0]
        @test KVec{2}("α,.5")() == KVec{2}("α,.5")(nothing) == [0,0.5]
        @test KVec{3}("u,v,0.5")([.7,.2,.3]) == [.7,.2,.5]
        @test KVec{3}("v,u,0.5")([.7,.2,.3]) == [.2,.7,.5]

        @test KVec([0.0,1//1,2]) == KVec{3}("0,1,2") == KVec{3}([0.0,1//1,2])
        @inferred KVec{3}([0.0,1//1,2])
        @test KVec{3}([0, 0.5, 0], [0 1 0; 0 1/2 1/2; -1/2 0 1/3]) == KVec{3}("β, 1/2+1/2β+1/2γ, -1/2α+1/3γ")

        @test zero(KVec{3}) == KVec{3}("0,0,0") == zero(KVec{3}("1,2,3"))
    end

    αβγ = Crystalline.TEST_αβγ
    @testset "Composition" begin
        @test SymOperation{2}("x,y") * KVec{2}("1,α") == KVec{2}("1,α")

        # 4⁺ rotation
        @test SymOperation{2}("-y,x") * KVec(.4,.3) == KVec(-.3,.4)

        # composition is associative for both `KVec`s and `RVec`s
        kv  = KVec("α, -α")
        kv′ = KVec("α, 1/2-2β")
        g   = S"-x+y,-x" # 3⁻
        h   = S"-x,-x+y" # m₂₁
        @test h * (g * kv)  == (h * g) * kv
        @test h * (g * kv′) == (h * g) * kv′
        
        rv  = RVec("α, -α")
        rv′ = RVec("1/4+α, 1/2-2β")
        g   = S"-x+y,-x" # 3⁻
        h   = S"-x,-x+y+1/3" # m₂₁
        @test h * (g * rv)  == (h * g) * rv
        @test h * (g * rv′) == (h * g) * rv′
    end

    @testset "Composition of SymmetryOperation and AbstractPoint" begin
        @test SymOperation{2}("x,y") * ReciprocalPoint(.2, .3) == ReciprocalPoint(.2, .3)
        @test SymOperation{2}("x,-y") * DirectPoint(.2, -.3) == DirectPoint(.2, .3)

        op = SymOperation{3}("-z,x,-y")
        k = ReciprocalPoint(.2, .15, .42)
        r = DirectPoint(.8, .23, -.37)
        @test k⋅r ≈ (op*k)⋅(op*r)

        op = SymOperation{3}("-z,x,-y+1/3")
        @test op*r ≈ DirectPoint(.37, .8, -.23+1/3)
    end
end