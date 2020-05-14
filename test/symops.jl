using Crystalline, Test

@testset "Symmetry operations, Bilbao" begin
    @testset "Space group #1" begin
        sg = spacegroup(1, 3)
        @test order(sg) == 1
        @test dim(sg) == 3
        op = sg[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test xyzt(op) == "x,y,z"
    end

    @testset "Space group #146" begin
        sg = spacegroup(146, 3)
        @test order(sg) == 9
        @test dim(sg) == 3
        op = sg[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test xyzt(op) == "-x+y+1/3,-x+2/3,z+2/3"
    end

    @testset "Plane group #7" begin
        plg = spacegroup(7, 2)
        @test order(plg) == 4
        @test dim(plg) == 2
        op = plg[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test xyzt(op) == "-x,-y"
    end

    @testset "Conversion between xyzt and matrix forms" begin
        for dim = 1:3
            for sgnum in 1:MAX_SGNUM[dim]
                    @testset "SG$sgnum ($(dim)D)" begin
                    sg = spacegroup(sgnum, dim)
                    for op in sg
                        @test all(xyzt2matrix(xyzt(op))   == matrix(op)) # xyzt->matrix
                        @test all(matrix2xyzt(matrix(op)) == xyzt(op))   # matrix->xyzt
                    end
                end
            end
        end
    end

    @testset "Composition" begin
        ops = spacegroup(230, 3) # random space group

        # test associativity (with and without modular arithmetic)
        g₁, g₂, g₃ = ops[5:7] 
        @test g₁∘(g₂∘g₃) == (g₁∘g₂)∘g₃
        @test compose(compose(g₁, g₂, false), g₃, false) == compose(g₁, compose(g₂, g₃, false), false)
    end

end