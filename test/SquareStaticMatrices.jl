using StaticArrays
using Crystalline
using Crystalline.SquareStaticMatrices
using Test

@testset "SquareStaticMatrices" begin
    A   = [1 2 3; 4 5 6; 7 8 9]
    SA  = SMatrix{3,3}(A)
    SSA = SqSMatrix(SA)
    
    @test (@inferred SqSMatrix(SA)) === SSA
    @test SqSMatrix{3}(A) === SSA
    @test SqSMatrix(A) === SSA

    @test SSA[1] == SA[1]
    @test SSA[1,1] == SA[1,1]
    @test SSA == SA
    @test SSA[1:end] == SA[1:end]
    @test convert(typeof(A), SSA) isa typeof(A)
    @test convert(SqSMatrix{3, Int}, eachcol(SSA)) === SSA
    @test SqSMatrix{3,Int}(A) === SqSMatrix{3}(A)
    @test SqSMatrix{3,Int}(A) == A
    @test (@inferred SMatrix(SSA)) === SA
    @test SMatrix{3,3}(SSA) === SA
    @test SMatrix{3,3,Int}(SSA) === SA
    @test SqSMatrix(SSA) === SSA === SqSMatrix{3,Int}(SSA)

    A′ = [1 2; 3 4]
    @test SqSMatrix(A′) === SqSMatrix{2}(A′)

    @test zero(SSA) == zero(A)
end

