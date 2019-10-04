using SGOps, Test

@testset "Symmetry operations" begin
    @testset "Space group #1" begin
        sg = get_sgops(1, 3; verbose=false)
        @test order(sg) == 1
        @test dim(sg) == 3
        op = operations(sg)[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test xyzt(op) == "x,y,z"
    end

    @testset "Space group #146" begin
        sg = get_sgops(146, 3; verbose=false)
        @test order(sg) == 9
        @test dim(sg) == 3
        op = operations(sg)[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test xyzt(op) == "-x+y+1/3,-x+2/3,z+2/3"
    end

    @testset "Plane group #7" begin
        plg = get_sgops(7, 2; verbose=false)
        @test order(plg) == 4
        @test dim(plg) == 2
        op = operations(plg)[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test xyzt(op) == "-x,-y"
    end

    @testset "Conversion between xyzt and matrix forms" begin
        for sgnum = 1:230
                @testset "SG$sgnum" begin
                sg = get_sgops(sgnum; verbose=false)
                for op in operations(sg)
                    @test all(xyzt2matrix(xyzt(op))   == matrix(op)) # xyzt->matrix
                    @test all(matrix2xyzt(matrix(op)) == xyzt(op))   # matrix->xyzt
                end
            end
        end
    end
end