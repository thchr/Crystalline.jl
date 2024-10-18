using Crystalline
using Test

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
end