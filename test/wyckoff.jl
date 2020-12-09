using Test, Crystalline
using Crystalline: constant, free

@testset "SiteGroup" begin
    neg_error_tol = 1e-15
    for D in 1:3
        Dᵛ = Val(D)
        for sgnum in 1:MAX_SGNUM[D]
            sg  = spacegroup(sgnum, Dᵛ)
            wps = wyckpos(sgnum, Dᵛ)
            for wp in wps
                g = SiteGroup(sg, wp)
                @test MultTable(g).isgroup

                qv = qvec(wp)
                # test that ops in `g` leave the Wyckoff position representatives invariant
                for op in g
                    qv′ = op∘qv
                    @test isapprox(constant(qv), constant(qv′))
                    @test isapprox(free(qv),     free(qv′))

                    # test that all the constant parts of the positions in the Wyckoff orbit all
                    # have coordinates in [0,1)
                    wpreps = orbit(g, wp) # wyckoff position representatives
                    @test all(wpreps) do wp
                        all(xyz->xyz≥(-neg_error_tol) && xyz<1, constant(qvec(wp)))
                    end
                end
                
            end
        end
    end
end