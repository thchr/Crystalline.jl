using Crystalline, Test

@testset "Symmetry operations" begin
    @testset "Basics" begin
        # Space group ⋕1
        sg = spacegroup(1, Val(3))
        @test order(sg) == 1
        @test dim(sg) == 3
        op = sg[1]
        @test matrix(op) == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0]
        @test xyzt(op) == "x,y,z"

        # Space group ⋕146
        sg = spacegroup(146, Val(3))
        @test order(sg) == 9
        @test dim(sg) == 3
        op = sg[9]
        @test matrix(op) ≈ [-1.0  1.0  0.0  1/3; -1.0  0.0  0.0  2/3; 0.0  0.0  1.0  2/3]
        @test xyzt(op) == "-x+y+1/3,-x+2/3,z+2/3"

        # Plane group ⋕7
        sg = spacegroup(7, 2) # keep as 2 (rather than Val(2)) intentionally, to test...
        @test order(sg) == 4
        @test dim(sg) == 2
        op = sg[2]
        @test matrix(op) ≈ [-1.0 0.0 0.0; 0.0 -1.0 0.0]
        @test xyzt(op) == "-x,-y"

        # Round-trippability of constructors
        op = SymOperation{3}("-y,x,z+1/2")
        @test op == SymOperation("-y,x,z+1/2")
        @test op == S"-y,x,z+1/2"
        @test op == SymOperation(rotation(op), translation(op)) 
        @test op == SymOperation(matrix(op)) # SMatrix
        @test op == SymOperation(Matrix(op)) # Matrix

        # identity operation
        @test S"x,y,z" == one(S"y,z,x") == one(SymOperation{3})
    end

    @testset "Parsing cornercases" begin
        # allow some forms of whitespace in xyzt string parsing (issue #29)
        @test SymOperation("x,y,z") == SymOperation("x, y, z")
        @test SymOperation("x,-y,+z") == SymOperation(" x, - y, +z")
        @test SymOperation("x,-y+x,+z-1/2") == SymOperation("x, - y + x, +z - 1/2")
        @test SymOperation("  x,-y+x+1/3,+z-1/2") == SymOperation("   x, - y + x + 1/3, +z - 1/2")
    end

    @testset "Conversion between xyzt and matrix forms" begin
        for D = 1:3
            Dᵛ = Val(D)
            @testset "D = $(D)" begin
                for sgnum in 1:MAX_SGNUM[D]
                    @testset "⋕$sgnum" begin
                        sg = spacegroup(sgnum, Dᵛ)
                        for op in sg
                            @test Crystalline.xyzt2matrix(xyzt(op), Dᵛ) ≈ matrix(op) # xyzt->matrix (`≈` due to possible rounding errors and precision loss on round-trip)
                            @test Crystalline.matrix2xyzt(matrix(op))   == xyzt(op)  # matrix->xyzt
                        end
                    end
                end
            end
        end
    end

    @testset "Composition" begin
        sg = spacegroup(230, Val(3)) # random space group

        # test associativity (with and without modular arithmetic)
        g₁, g₂, g₃ = sg[5:7]
        @test g₁*(g₂*g₃) == (g₁*g₂)*g₃
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
        cartRs_from_cRs = cartesianize(csg, cRs)
        cartRs_from_pRs = cartesianize(psg, pRs)

        @test all(isapprox.(cartRs_from_cRs, cartRs_from_pRs, atol=1e-12))

        # `isapprox` with centering
        op  = S"-y+2/3,x-y+1/3,z+1/3" # {3₀₀₁⁺|⅔,⅓,⅓} (in a rhombohedral system)
        op′ = op * SymOperation{3}(Bravais.centeringtranslation('R',Val(3)))
        @test isapprox(op, op′, 'R', true)
        @test !isapprox(op, op′, 'P', true)

        # `isapprox` with `modw = false`
        op = S"x,-y,z"
        op′ = S"x-3,-y+10,z-5"
        @test !isapprox(op, op′, 'I', false)
        @test !isapprox(op, op′, 'P', false)
        @test isapprox(op, op′, 'P', true)
        @test isapprox(op, op′, 'I', true)
    end

    @testset "Groups created from generators" begin # (`generate` default sorts by `seitz`)
        # generate plane group (17) p6mm from 6⁺ and m₁₀
        gens = SymOperation.(["x-y,x", "-x+y,y"]) 
        sg = spacegroup(17, Val(2))
        @test sort!(generate(gens)) == sort!(sg)

        # generate site symmetry group of Wyckoff position 2b in p6mm
        ops  = SymOperation.(
                ["x,y","-y+1,x-y+1", "-x+y,-x+1",    # {1|0}, {3⁺|1.0,1.0}, {3⁻|0,1.0}, 
                "-y+1,-x+1", "-x+y,y", "x,x-y+1"])   # {m₁₁|1.0,1.0}, {m₁₀|0}, {m₀₁|0,1.0}
        gens = ops[[2,6]]
        @test sort!(generate(gens, modτ=false), by=xyzt) == sort!(ops, by=xyzt)

        # generate space group with nonsymmorphic operations
        sg = spacegroup(180, Val(3)) # P6₂22
        gens = sg[[6,8]]             # {6₀₀₁⁺|0,0,⅓}, 2₁₀₀
        @test sort!(generate(gens)) ≈ sort!(sg)

        # generators do not specify a finite group under "non-modulo" composition
        @test_throws OverflowError generate(SymOperation.(["x,y+1,z"]); modτ=false, Nmax=50)
    end

    @testset "Generators" begin
        for (Dᵛ, gtype) in ((Val(1), SpaceGroup{1}), (Val(2), SpaceGroup{2}), (Val(3), SpaceGroup{3}))
            D = typeof(Dᵛ).parameters[1]
            for sgnum in 1:MAX_SGNUM[D]
                ops1 = sort!(spacegroup(sgnum, Dᵛ))
                ops2 = sort!(generate(generators(sgnum, gtype)))
                @test ops1 ≈ ops2
            end
        end
    end

    @testset "Error types and domain checking" begin
        @test_throws DomainError spacegroup(231, 3)    
        @test_throws DomainError spacegroup(-1,  2)
        @test_throws DomainError spacegroup(2,   0)
        @test_throws DomainError spacegroup(41,  5)
        @test_throws ArgumentError SymOperation{2}("x,z")
        @test_throws ArgumentError SymOperation("x,z")
        @test_throws ArgumentError SymOperation("x,÷z")
        @test_throws ArgumentError SymOperation("x ,z") # don't allow spaces *after* entries
        @test_throws DimensionMismatch SymOperation{3}("x,y+z")
    end

    @testset "Checking symmorphic space groups" begin
        # we do a memoized look-up for `issymmorph(::Integer, ::Integer)`: ensure that it 
        # agrees with explicit calculations
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                @test issymmorph(sgnum, D) == issymmorph(spacegroup(sgnum, D))
            end
        end
    end
end

@testset "coset representatives" begin
    G = pointgroup("6mm")
    H = pointgroup("3")
    Q = cosets(G, H) 
    @test length(Q) == div(order(G), order(H))
    cG = reduce(vcat, (compose.(Ref(q), H) for q in Q))
    @test sort(cG, by=xyzt) == sort(G, by=xyzt)

    G = spacegroup(230)
    H  = spacegroup(206) # a maximal subgroup of G, same setting
    H′ = spacegroup(199) # a maximal subgroup of H, also same setting
    Q  = cosets(G, H)
    Q′ = cosets(G, H)
    cG  = reduce(vcat, (compose.(Ref(q), H) for q in Q))
    cG′ = reduce(vcat, (compose.(Ref(q), H) for q in Q′))
    @test sort(cG,  by=xyzt) == sort(G, by=xyzt)
    @test sort(cG′, by=xyzt) == sort(G, by=xyzt)
end