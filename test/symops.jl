using SGOps, Test

@testset "Symmetry operations" begin
    @testset "Space group #1" begin
        sg = get_symops(1, 3; verbose=false)
        @test order(sg) == 1
        @test dim(sg) == 3
        op = operations(sg)[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test shorthand(op) == "x,y,z"
    end

    @testset "Space group #146" begin
        sg = get_symops(146, 3; verbose=false)
        @test order(sg) == 9
        @test dim(sg) == 3
        op = operations(sg)[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test shorthand(op) == "-x+y+1/3,-x+2/3,z+2/3"
    end

    @testset "Plane group 7" begin
        pg = get_symops(7, 2; verbose=false)
        @test order(pg) == 4
        @test dim(pg) == 2
        op = operations(pg)[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test shorthand(op) == "-x,-y"
    end
end