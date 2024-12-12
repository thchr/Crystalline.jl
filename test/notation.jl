using Crystalline, Test

@testset "Notation" begin
    @testset "Mulliken notation" begin
        for D in 1:3
            for pglab in PG_IUCs[D]
                pgirs′ = pgirreps(pglab, Val(D))
                g      = group(first(pgirs′))

                # check both "ordinary" irreps and physically real irreps (coreps)
                for irtype in (identity, realify)
                    pgirs = irtype(pgirs′)
                    mlabs = mulliken.(pgirs) # Mulliken symbols for point group irreps

                    # each irrep must be assigned a unique symbol
                    @test allunique(mlabs)

                    # check irrep dimension vs. label
                    irDs = irdim.(pgirs)
                    mlabs_main = getfield.(match.(Ref(r"A|B|¹E|²E|E|T"), mlabs), Ref(:match))
                    for (irD, mlab_main) in zip(irDs, mlabs_main)
                        if irD == 1
                            @test mlab_main ∈ ("A", "B", "¹E", "²E")
                        elseif irD == 2
                            @test mlab_main == "E"
                        elseif irD == 3
                            @test mlab_main == "T"
                        else
                            error("Unexpectedly got point group irrep of dimension higher than 3")
                        end
                    end
                    # check 'u' and 'g' subscripts for transformation under inversion for 3d
                    # point group irreps (doesn't make complete sense for 2D and 1D, where 
                    # inversion is really an embedded" 2-fold rotation (2D) or a mirror (1D))
                    if D == 3
                        if (idx⁻¹ = findfirst(op -> seitz(op) == "-1", g)) !== nothing
                            for (mlab, pgir) in zip(mlabs, pgirs)
                                irD = irdim(pgir)
                                χ⁻¹ = characters(pgir)[idx⁻¹]
                                if χ⁻¹ ≈ -irD # antisymmetric under inversion => ungerade
                                    @test occursin('ᵤ', mlab)
                                elseif χ⁻¹ ≈ +irD # symmetric under inversion => gerade
                                    @test occursin('g', mlab)
                                else
                                    throw(DomainError(χ_inversion, "inversion character must be equal to ±(irrep dimension)"))
                                end
                            end
                        else
                            @test all(mlab -> !occursin('ᵤ', mlab) && !occursin('g', mlab), mlabs)
                        end
                    end
                end
            end
        end
    end

    @testset "PGIrreps in Mulliken notation" begin
        # 2D & 3D
        for D in (2,3)
            pgirs          = pgirreps("6", Val(D); mulliken=false)
            pgirs_mulliken = pgirreps("6", Val(D); mulliken=true)
            @test label.(pgirs_mulliken) == mulliken.(pgirs)

            # we do some reductions of Mulliken labels during "realification"; check that
            # this works the same across `label.(realify(pgirreps(..; mulliken=true)))`
            # and `mulliken.(realify(pgirreps(..)))`
            pgirs′          = realify(pgirs)
            pgirs_mulliken′ = realify(pgirs_mulliken)
            @test label.(pgirs_mulliken′) == mulliken.(pgirs′)
        end
        # 1D
        @test label.(pgirreps("m", Val(1); mulliken=true)) == mulliken.(pgirreps("m", Val(1))) == ["A′", "A′′"]
    end

    @testset "Centering, Bravais types, & IUC" begin
        for D in 1:3
            for sgnum in 1:MAX_SGNUM[D]
                @test first(iuc(sgnum, D)) == centering(sgnum, D) == last(bravaistype(sgnum, D; normalize=false))
            end
        end
    end
end