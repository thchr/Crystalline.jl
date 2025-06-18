using Crystalline
using Test

@testset "Isomorphic parent point group identification" begin
@testset "Site symmetry groups & parent point group identification" begin
    # here, we simply test that our identifications of the label associated with the site
    #symmetry group's isomorphic point group are consistent with those made in the
    # identification of the site symmetry group labels for band representations from the
    # Bilbao Crystallographic Server
    for sgnum in 1:MAX_SGNUM[3]
        brs = bandreps(sgnum, 3)
        sg  = spacegroup(sgnum, Val(3))
        sitegs = findmaximal(sitegroups(sg))
        for siteg in sitegs           
            # "our" label
            pg, _, _ = find_isomorphic_parent_pointgroup(siteg)
            pglab = label(pg)
            pglab = pglab ∈ ("312", "321")   ? "32"  : # normalize computed label (since
                    pglab ∈ ("-31m", "-3m1") ? "-3m" : # Bilbao doesn't distinguish 
                    pglab ∈ ("31m", "3m1")   ? "3m"  : # between e.g., 312 and 321)
                    pglab

            # Bilbao's label
            idx = something(findfirst(br->br.wyckpos == label(siteg.wp), brs))
            pglab_bilbao = brs[idx].sitesym

            # compare
            @test pglab == pglab_bilbao
        end
    end
end # @testset "Site symmetry groups & parent point group identification"

function Γ_littlegroup(brs::Collection{<:NewBandRep})
    lgirsv = irreps(brs)
    idx = something(findfirst(lgirs -> klabel(lgirs) == "Γ", lgirsv))
    return group(lgirsv[idx])
end
@testset "Corner cases (primitivized little group)" begin
    for sgnum in 1:230
        # we only test those space groups that are not "P" centered
        centering(sgnum, 3) == 'P' && continue
        
        if sgnum ∈ (225, 226, 227, 228, 229, 230) # times out
            timesout_token = false
            @test_broken timesout_token # calculation times out in isomorphism search
            continue
        end

        brs = calc_bandreps(sgnum, Val(3))
        lg = Γ_littlegroup(brs) # little group at Γ
        plg, _, _ = find_isomorphic_parent_pointgroup(lg)

        # now primitivize `brs` and check that we can still find the same parent point group
        brsᵖ = primitivize(brs)
        lgᵖ = Γ_littlegroup(brsᵖ) # little group at Γ
        
        if sgnum ∈ (196, 202, 203, 209, 210, 216, 219)
            # the "_sense_ + rotation" grouping leads us astray here and we don't manage
            # to find the parent point group that we ought to find
            errors_token = try
                find_isomorphic_parent_pointgroup(lgᵖ)
                true
            catch
                false
            end
            @test_broken errors_token
            continue
        end

        # works in the following casees, modulo indeterminacy around -42m and -4m2
        plgᵖ, _, _ = find_isomorphic_parent_pointgroup(lgᵖ)
        @test (label(plg) == label(plgᵖ) || (label(plg) ∈ ("-42m", "-4m2") &&
                                             label(plgᵖ) ∈ ("-42m", "-4m2")   ))
    end
end
end # @testset "Isomorphic parent point group identification"