using Test, Crystalline
using Crystalline: constant, free

@testset "SiteGroup" begin
    neg_error_tol = 1e-15
    for D in 1:3
        Dᵛ = Val(D)
        for sgnum in 1:MAX_SGNUM[D]
            sg  = spacegroup(sgnum, Dᵛ)
            wps = wyckoffs(sgnum, Dᵛ)
            for wp in wps
                g = SiteGroup(sg, wp)
                @test g isa SiteGroup

                qv = parent(wp)
                # test that ops in `g` leave the Wyckoff position `wp` invariant
                for op in g
                    qv′ = op*qv
                    @test isapprox(qv, qv′, nothing, false)
                end

                # test that all the constant parts of the positions in the Wyckoff orbit
                # all have coordinates in [0,1)
                wpreps = orbit(g, wp) # wyckoff position representatives
                @test all(wpreps) do wp
                    all(xyz->xyz≥(-neg_error_tol) && xyz<1, constant(wp))
                end

                # test that `g` and `cosets(g)` furnishes a left-coset decomposition of `sg`
                ops = [opʰ*opᵍ for opʰ in cosets(g) for opᵍ in g];
                @test sort!(ops, by=xyzt) ≈ sort(sg, by=xyzt)
            end
        end
    end
end

@testset "Maximal Wyckoff positions" begin
    for sgnum in 1:MAX_SGNUM[3]
        sg  = spacegroup(sgnum, Val(3))
        wps = wyckoffs(sgnum, Val(3))
        sitegs = SiteGroup.(Ref(sg), wps)

        max_sitegs = findmaximal(sitegs)
        max_wps    = wyck.(max_sitegs)

        # the band representations should include all maximal wyckoff positions; 
        # check consistency against that
        max_wps_brs_str = getfield.(bandreps(sgnum, 3).bandreps, Ref(:wyckpos))
        @test sort(unique(max_wps_brs_str)) == sort(label.(max_wps))
    end

end