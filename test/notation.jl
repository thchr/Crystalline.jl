using Crystalline, Test


@testset "Notation" begin
    SqSMatrix = Crystalline.SquareStaticMatrices.SqSMatrix
    @testset "Mulliken notation" begin
        for D in (1,2,3)
            for pglab in PGS_IUCs[D]
                pgirs = get_pgirreps(pglab, Val(D))
                mlabs = mulliken.(pgirs) # Mulliken symbols for point group irreps
                
                # each irrep must be assigned a unique symbol
                @test allunique(mlabs)

                # check 'u' and 'g' subscripts for transformation under inversion
                g  = group(first(pgirs))
                if (idx⁻¹ = findfirst(op -> isone(-rotation(op)), g)) !== nothing
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

                # check agreement with Bilbao's conventions
                # TODO
            end
        end
    end
end