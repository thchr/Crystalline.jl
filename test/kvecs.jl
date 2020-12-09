using Test, Crystalline, StaticArrays

@testset "KVec" begin
    @test dim(KVec("u")) == 1
    @test dim(KVec("0,0.5")) == 2
    @test dim(KVec("0.3,0.2,.1")) == 3
    @test parts(KVec("u,v,w")) == ([0,0,0], [1 0 0; 0 1 0; 0 0 1])

    @test KVec{1}("-u") == KVec{1}("-α") == KVec(zeros(1), fill(-1,1,1)) == KVec(zero(SVector{1,Float64}), -one(SMatrix{1,1,Float64}))

    @test string(KVec{3}("u+1,0,0.5")) == string(KVec{3}("1.0+u,0,0.5")) == "[1.0+α, 0.0, 0.5]"
    @test string(KVec{3}("u+v,u,0"))   == string(KVec{3}("α+β,α,0"))     == "[α+β, α, 0.0]"
    @test string(KVec{3}("-0.3,0.2,-.1")) == "[-0.3, 0.2, -0.1]"

    @test string(KVec{3}("α, 2α, 0.0"))   == string(KVec{3}("α, 2.0α, 0.0")) == "[α, 2.0α, 0.0]"
    @test string(KVec{3}("α, 2α+2, 0.0")) == string(KVec{3}("α, 2.0α+2.0, 0.0")) == "[α, 2.0+2.0α, 0.0]"

    @test KVec{3}("x,y,z")(0,1,2) == [0.0,1.0,2.0]
    @test KVec{2}("α,.5")() == KVec{2}("α,.5")(nothing) == [0,0.5]
    @test KVec{3}("u,v,0.5")([.7,.2,.3]) == [.7,.2,.5]
    @test KVec{3}("v,u,0.5")([.7,.2,.3]) == [.2,.7,.5]
end