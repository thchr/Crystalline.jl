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
    @test string(zero(n)) == "[0Mᵢ, 0Xᵢ, 0Γᵢ, 0Rᵢ] (0 bands)"
end # @testset "(Abstract)SymmetryVectors"

@testset "CompositeBandRep & @composite" begin
    # basic example with `@composite`: agree with `SymmetryVector`
    brs = calc_bandreps(2, Val(3))

    cbr = @composite brs[1] - brs[2] + 2brs[5] - brs[3]*7 - brs[4]*(-8)
    n = brs[1] - brs[2] + 2brs[5] - brs[3]*7 - brs[4]*(-8)

    @test cbr isa CompositeBandRep{3}
    @test cbr == n
    @test SymmetryVector(cbr) == n

    @test size(cbr)        == size(n)
    @test occupation(cbr)  == occupation(n)
    @test irreps(cbr)      == irreps(n)      == irreps(brs)
    @test num(cbr)         == num(n)         == num(brs)
    @test klabels(cbr)     == klabels(n)     == klabels(brs)
    @test irreplabels(cbr) == irreplabels(n) == irreplabels(brs)

    @test cbr - brs[1] == n - brs[1]

    # more challenging examples with `@composite`
    @test (@composite -brs[1]) == -brs[1]
    @test (@composite -brs[end]) == -brs[end]
    @test (@composite -brs[1+end-2]) == -brs[1+end-2]
    @test (@composite brs[1]) == brs[1] == (@composite brs[begin])
    @test ((@composite -brs[1] + brs[2]*3 - brs[7]*(-2) + brs[end] - (-3)*brs[end-2]) == 
           (-brs[1] + brs[2]*3 - brs[7]*(-2) + brs[end] - (-3)*brs[end-2]))
end