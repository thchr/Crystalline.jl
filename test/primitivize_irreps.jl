using Crystalline, Test

@testset "`primitivize` for collections of irreps" begin
    # LGIrrep{3}
    sgnum = 227 # F-centered 
    lgirsd = lgirreps(sgnum, Val(3))
    lgirsd′ = primitivize(lgirsd) # primitivize of Dict
    for (klab, lgirs′) in lgirsd′
        klab ∈ ("Λ", "L") && continue # maps to itself under F-centering primitivize; so skip
        lgirs = lgirsd[klab]
        lgirs′′ = primitivize(lgirs) # primitivize of Collection, not Dict
        for (lgir, lgir′, lgir′′) in zip(lgirs, lgirs′, lgirs′′)
            @test group(lgir′′) == group(lgir′) && lgir′.matrices == lgir′′.matrices
            @test group(lgir) ≠ group(lgir′) # the underlying group should now be different
        end
    end

    sgnum = 2 # P-centered; already primitive
    lgirsd = lgirreps(sgnum, Val(3))
    lgirsd′ = primitivize(lgirsd) # primitivize of Dict
    for (klab, lgirs′) in lgirsd′
        lgirs = lgirsd[klab]
        lgirs′′ = primitivize(lgirs) # primitivize of Collection
        for (lgir, lgir′, lgir′′) in zip(lgirs, lgirs′, lgirs′′)
            @test group(lgir′′) == group(lgir′) && lgir′.matrices == lgir′′.matrices
            @test group(lgir) == group(lgir′) # P-centered, so should be unchanged
        end
    end

    # LGIrrep{2}
    sgnum = 5 # c-centered
    lgirsd = lgirreps(sgnum, Val(2))
    lgirsd′ = primitivize(lgirsd) # primitivize of Dict
    for (klab, lgirs′) in lgirsd′
        lgirs = lgirsd[klab]
        lgirs′′ = primitivize(lgirs) # primitivize of Collection, not Dict
        for (lgir, lgir′, lgir′′) in zip(lgirs, lgirs′, lgirs′′)
            @test group(lgir′′) == group(lgir′) && lgir′.matrices == lgir′′.matrices
            @test group(lgir) ≠ group(lgir′) # the underlying group should now be different
        end
    end

    # SiteIrrep{3}
    sgnum = 230 # I-centered
    sitegs = sitegroups(sgnum, Val(3))
    for siteg in sitegs
        siteg′ = primitivize(siteg)
        siteirs = siteirreps(siteg)
        siteirs′ = primitivize(siteirs)
        for siteir′ in siteirs′
            @test siteg′ == group(siteir′)
            # there is one site group here (16a) where the symmetry ops are sufficiently
            # simple that the same operations result before/after primitivization (think,
            # e.g., of just S"-x,-y,-z", which is always mapped to itself); skip these
            # for test below
            if label(position(siteg)) == "16a"
                continue
            end
            @test siteg ≠ group(siteir′)   # the underlying group should now be different
        end
    end
end