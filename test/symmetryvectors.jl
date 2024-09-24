using Test
using Crystalline

@testset "(Abstract)SymmetryVectors" begin
    brs = calc_bandreps(221) # ::Collection{NewBandRep{3}}

    # addition of symmetry vectors
    @test brs[1] + brs[2] == SymmetryVector(brs[1]) + SymmetryVector(brs[2])
    @test Vector(brs[1] + brs[2]) == Vector(brs[1]) + Vector(brs[2])

    # multiplication of symmetry vectors
    @test brs[1] + brs[1] == 2brs[1]
    @test -brs[1] == 2brs[1] - 3brs[1]
    @test zero(brs[1]) == 0*brs[1]
    @test brs[3]*7 == 7*brs[3]

    # type-stable summation
    @test sum(brs[1:1]) isa SymmetryVector{3}
    @test sum(brs[1:2]) isa SymmetryVector{3}

    # printing SymmetryVector that have zero-contents
    n = SymmetryVector(brs[1])
    @test n*0 == zero(n)
    @test string(n) == "[0Mᵢ, 0Xᵢ, 0Γᵢ, 0Rᵢ] (0 bands)"
end # @testset "(Abstract)SymmetryVectors"