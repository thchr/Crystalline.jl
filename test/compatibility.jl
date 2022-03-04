using Test, Crystalline, StaticArrays

@testset "Compatibility" begin
    lgir218 = realify(lgirreps(218)["R"])[end]; # 6D irrep R₄R₅²¹⁸

    # test that complex corep subduces to itself
    @test subduction_count(lgir218, lgir218) == 1

    # test that complex corep subduces to its components
    lgirs218_notr = lgirreps(218)["R"]
    @test [subduction_count(lgir218, lgir) for lgir in lgirs218_notr] == [0,0,0,1,1]

    # test that two distinct coreps in different groups decompose consistenly
    lgir222 = realify(lgirreps(222)["R"])[end]; # 6D irrep R₄²²²
    lgir218′ = deepcopy(lgir218)                # 6D irrep R₄R₅²¹⁸
    
    # SG 218 is a subgroup of 222, but in a different setting - so convert the little group of 218 to that of 222
    lg_ops = operations(group(lgir218′)) 
    lg_ops .= transform.(lg_ops, Ref(one(SMatrix{3,3})), Ref((@SVector [1,1,1])/4))

    # check how the R₄²²² subduces into R₄R₅²¹⁸
    @test subduction_count(lgir222, lgir218′) == 1
end