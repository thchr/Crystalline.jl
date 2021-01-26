using Crystalline, Test

@testset "Notation" begin
    SqSMatrix = Crystalline.SquareStaticMatrices.SqSMatrix
    @testset "Mulliken notation" begin
        for D in 1:3
            for pglab in PGS_IUCs[D]
                pgirs = get_pgirreps(pglab, Val(D))
                g     = group(first(pgirs))
                mlabs = mulliken.(pgirs) # Mulliken symbols for point group irreps

                # each irrep must be assigned a unique symbol
                @test allunique(mlabs)

                # check irrep dimension vs. label
                irDs = Crystalline.irdim.(pgirs)
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
                            irD = Crystalline.irdim(pgir)
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