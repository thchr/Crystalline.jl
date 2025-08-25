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
                g = sitegroup(sg, wp)
                @test g isa SiteGroup

                rv = parent(wp)
                # test that ops in `g` leave the Wyckoff position `wp` invariant
                for op in g
                    rv′ = op*rv
                    @test isapprox(rv, rv′, nothing, false) # isapprox(::RVec, ::RVec)
                    wp′ = op*wp
                    @test isapprox(wp, wp′, nothing, false) # isapprox(::WyckoffPosition, ::WyckoffPosition)
                end

                # test that all the constant parts of the positions in the Wyckoff orbit
                # all have coordinates in [0,1)
                wpreps = orbit(g) # wyckoff position representatives
                @test all(wpreps) do wp
                    all(xyz->xyz≥(-neg_error_tol) && xyz<1, constant(wp))
                end

                # test that `g` and `cosets(g)` furnishes a left-coset decomposition of `sg`
                ops = [opʰ*opᵍ for opʰ in cosets(g) for opᵍ in g];
                @test sort!(ops, by=xyzt) ≈ sort(sg, by=xyzt)

                # test that primitivize works & has correct orbit/coset length
                N_cntr = Bravais.centering_volume_fraction(centering(sgnum, D), Val(D))
                @test length(orbit(primitivize(g))) == multiplicity(wp) ÷ N_cntr
            end
        end
    end
end

@testset "Maximal Wyckoff positions" begin
    for sgnum in 1:MAX_SGNUM[3]
        sg  = spacegroup(sgnum, Val(3))
        sitegs = sitegroups(sg)

        max_sitegs = findmaximal(sitegs)
        max_wps    = position.(max_sitegs)

        # the band representations should include all maximal wyckoff positions; 
        # check consistency against that
        brs = bandreps(sgnum, 3)
        max_wps_brs_str = map(_br -> _br.wyckpos, brs.bandreps)
        @test sort(unique(max_wps_brs_str)) == sort(label.(max_wps))
    end

    # type-stability
    @test (@inferred Vector{WyckoffPosition{1}} wyckoffs(1, Val(1))) isa Vector{WyckoffPosition{1}}
    @test (@inferred Vector{WyckoffPosition{2}} wyckoffs(1, Val(2))) isa Vector{WyckoffPosition{2}}
    @test (@inferred Vector{WyckoffPosition{3}} wyckoffs(1, Val(3))) isa Vector{WyckoffPosition{3}}
end