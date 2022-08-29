using Test, Crystalline

@testset "Cached orders:" begin
    @testset "Point groups" begin
        for (D,Dᵛ) in ((1,Val(1)), (2,Val(2)), (3,Val(3)))
            for iucs in Crystalline.PG_NUM2IUC[D]
                for iuc in iucs
                    @test Crystalline.PG_ORDERs[iuc] == order(pointgroup(iuc, Dᵛ))
                end
            end
        end
    end

    @testset "Space groups" begin
        for (D,Dᵛ) in ((1,Val(1)), (2,Val(2)), (3,Val(3)))
            for n in 1:MAX_SGNUM[D]
                sg = spacegroup(n, Dᵛ)
                psg = primitivize(sg)
                @test Crystalline.SG_ORDERs[D][n] == order(sg)
                @test Crystalline.SG_PRIMITIVE_ORDERs[D][n] == order(psg)
            end
        end
    end

    @testset "Subperiodic groups" begin
        subgind = Crystalline.subg_index
        for (D,P,Dᵛ,Pᵛ) in ((2,1,Val(2),Val(1)), (3,1,Val(3),Val(1)), (3,2,Val(3),Val(2)))
            for n in 1:MAX_SUBGNUM[(D,P)]
                subg = subperiodicgroup(n, Dᵛ, Pᵛ)
                psubg = primitivize(subg)
                @test Crystalline.SUBG_ORDERs[subgind(D,P)][n] == order(subg)
                @test Crystalline.SUBG_PRIMITIVE_ORDERs[subgind(D,P)][n] == order(psubg)
            end
        end
    end
end
