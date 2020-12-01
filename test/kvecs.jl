using Test, Crystalline

@testset "KVec" begin
    @test dim(KVec("u")) == 1
    @test dim(KVec("0,0.5")) == 2
    @test dim(KVec("0.3,0.2,.1")) == 3
    @test parts(KVec("u,v,w")) == ([0,0,0], [1 0 0; 0 1 0; 0 0 1])

    @test KVec("-u") == KVec("-α") == KVec(zeros(1), fill(-1,1,1))

    @test string(KVec("u+1,0,0.5")) == string(KVec("1.0+u,0,0.5")) == "[1.0+α, 0.0, 0.5]"
    @test string(KVec("u+v,u,0"))   == string(KVec("α+β,α,0"))     == "[α+β, α, 0.0]"
end