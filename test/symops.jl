using Crystalline, Test

@testset "Symmetry operations" begin
    @testset "Basics" begin
        # Space group #1
        sg = spacegroup(1, Val(3))
        @test order(sg) == 1
        @test dim(sg) == 3
        op = sg[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test xyzt(op) == "x,y,z"

        # Space group #146
        sg = spacegroup(146, Val(3))
        @test order(sg) == 9
        @test dim(sg) == 3
        op = sg[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test xyzt(op) == "-x+y+1/3,-x+2/3,z+2/3"

        # Plane group #7
        plg = spacegroup(7, 2) # keep as 2 (rather than Val(2)) intentionally, to test...
        @test order(plg) == 4
        @test dim(plg) == 2
        op = plg[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test xyzt(op) == "-x,-y"

        # Round-trippability of constructors
        op = SymOperation{3}("-y,x,z+1/2")
        @test op == S"-y,x,z+1/2"
        @test op == SymOperation(rotation(op), translation(op)) 
        @test op == SymOperation(matrix(op)) # SMatrix
        @test op == SymOperation(Matrix(op)) # Matrix
    end

    @testset "Conversion between xyzt and matrix forms" begin
        for dim = 1:3
            for sgnum in 1:MAX_SGNUM[dim]
                    @testset "SG$sgnum ($(dim)D)" begin
                    sg = spacegroup(sgnum, dim)
                    for op in sg
                        @test all(Crystalline.xyzt2matrix(xyzt(op)) == matrix(op)) # xyzt->matrix
                        @test all(Crystalline.matrix2xyzt(matrix(op)) == xyzt(op)) # matrix->xyzt
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

    @testset "Operators in differents bases" begin
        sgnum = 110 # bravais type "tI"
        cntr = centering(sgnum, 3) # 'I'

        csg = spacegroup(sgnum, Val(3))                             # conventional basis
        psg = SpaceGroup{3}(sgnum, primitivize.(csg, cntr, false))  # primitive basis

        # compute a random possible basis for tI bravais type
        cRs = directbasis(sgnum, Val(3))                            # conventional basis
        pRs = primitivize(cRs, cntr)                                # primitive basis

        # check that the cartesian representation of the space group operations, obtained
        # from `csg` & `cRs` versus `psg` and `pRs` agree (as it should)
        cartRs_from_cRs = Crystalline.cartesianize(csg, cRs)
        cartRs_from_pRs = Crystalline.cartesianize(psg, pRs)

        @test all(isapprox.(cartRs_from_cRs, cartRs_from_pRs, atol=1e-12))
    end

    @testset "Groups created from generators" begin # (`generate` default sorts by `seitz`)
        # generate plane group (17) p6mm from 6⁺ and m₁₀
        gens = SymOperation.(["x-y,x", "-x+y,y"]) 
        @test Set(generate(gens)) == Set(spacegroup(17,2))

        # generate site symmetry group of Wyckoff position 2b in p6mm
        ops  = SymOperation.(
                ["x,y","-y+1,x-y+1", "-x+y,-x+1",    # {1|0}, {3⁺|1.0,1.0}, {3⁻|0,1.0}, 
                "-y+1,-x+1", "-x+y,y", "x,x-y+1"])   # {m₁₁|1.0,1.0}, {m₁₀|0}, {m₀₁|0,1.0}
        gens = ops[[2,6]]
        @test Set(generate(gens, modτ=false)) == Set(ops)

        # generators do not specify a finite group under "non-modulo" composition
        @test_throws OverflowError (generate(SymOperation.(["x,y+1,z"]); modτ=false, Nmax=50))
    end

    @testset "Error types and domain checking" begin
        @test_throws DomainError spacegroup(231, 3)    
        @test_throws DomainError spacegroup(-1,  2)
        @test_throws DomainError spacegroup(2,   0)
        @test_throws DomainError spacegroup(41,  5)
    end

    @testset "Checking symmorphic space groups" begin
        # we do a memoized look-up for `issymmorph(::Integer, ::Integer)`: ensure that it 
        # agrees with explicit calculations
        for D in (1,2,3)
            for sgnum in 1:MAX_SGNUM[D]
                @test issymmorph(sgnum, D) == issymmorph(spacegroup(sgnum, D))
            end
        end
    end
end