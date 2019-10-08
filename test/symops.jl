using SGOps, Test

@testset "Symmetry operations, Bilbao" begin
    @testset "Space group #1" begin
        sg = get_sgops(1, 3)
        @test order(sg) == 1
        @test dim(sg) == 3
        op = operations(sg)[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test xyzt(op) == "x,y,z"
    end

    @testset "Space group #146" begin
        sg = get_sgops(146, 3)
        @test order(sg) == 9
        @test dim(sg) == 3
        op = operations(sg)[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test xyzt(op) == "-x+y+1/3,-x+2/3,z+2/3"
    end

    @testset "Plane group #7" begin
        plg = get_sgops(7, 2)
        @test order(plg) == 4
        @test dim(plg) == 2
        op = operations(plg)[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test xyzt(op) == "-x,-y"
    end

    @testset "Conversion between xyzt and matrix forms" begin
        sgnum_span = Dict(2=>1:17, 3=>1:230)
        for dim = 2:3
            for sgnum = sgnum_span[dim]
                    @testset "SG$sgnum" begin
                    sg = get_sgops(sgnum, dim)
                    for op in operations(sg)
                        @test all(xyzt2matrix(xyzt(op))   == matrix(op)) # xyzt->matrix
                        @test all(matrix2xyzt(matrix(op)) == xyzt(op))   # matrix->xyzt
                    end
                end
            end
        end
    end
end