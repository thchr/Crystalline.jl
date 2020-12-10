using Test, Crystalline
using Crystalline: constant, free

@testset "SiteGroup" begin
    neg_error_tol = 1e-15
    for D in 1:3
        Dᵛ = Val(D)
        for sgnum in 1:MAX_SGNUM[D]
            sg  = spacegroup(sgnum, Dᵛ)
            wps = get_wycks(sgnum, Dᵛ)
            for wp in wps
                g = SiteGroup(sg, wp)
                @test MultTable(g).isgroup

                qv = qvec(wp)
                # test that ops in `g` leave the Wyckoff position `wp` invariant
                for op in g
                    qv′ = op∘qv
                    @test isapprox(constant(qv), constant(qv′))
                    @test isapprox(free(qv),     free(qv′))
                end

                # test that all the constant parts of the positions in the Wyckoff orbit
                # all have coordinates in [0,1)
                wpreps = orbit(g, wp) # wyckoff position representatives
                @test all(wpreps) do wp
                    all(xyz->xyz≥(-neg_error_tol) && xyz<1, constant(qvec(wp)))
                end

                # test that `g` and `cosets(g)` furnishes a left-coset decomposition of `sg`
                ops = [opʰ∘opᵍ for opʰ in cosets(g) for opᵍ in g];
                @test sort!(seitz.(ops)) == sort!(seitz.(sg))
                # (the above test is analogous to checking `Set(ops) == Set(sg)`, but the
                #  `==`-check is not robust due to rounding errors here: unfortunately,
                #  there isn't currently a Base method for checking approximate equality
                #  of `Set`s)
            end
        end
    end
end