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

    # Correctly primitivizing also in the presence of a nonzero `lgir.translations` field,
    # where k-points are mapped. This only ever happens at nonspecial k-points; with e.g.
    # SG 214 at the Δ point being an example
    sgnum = 214
    lgirsd = lgirreps(sgnum, Val(3))
    lgirs = lgirsd["Δ"]
    lgirs′ = primitivize(lgirs)
    β = 0.314 # free parameter in Δ = (0, β, 0) [conv.] & Δ′ = [β, -β, β]/2 [primitive]
    αβγ = [0, β, 0]
    k = position(lgirs)(αβγ...)   # conventional basis
    k′ = position(lgirs′)(αβγ...) # primitive basis
    @test k == [0,β,0] && k′ ≈ [β,-β,β]/2
    for (lgir, lgir′) in zip(lgirs, lgirs′)
        # evaluate each at (α,β,γ) = (0,β,0) (specified in "basis" of free parameters, not 
        # in k-space per se)
        @test lgir(αβγ) ≈ lgir′(αβγ) # representation matrices should be the same

        # equivalently, the irrep matrices should _print_ the same, since both original and
        # primitivized versions refer to the same α, β, γ parameters: easiest to check by
        # comparing `Crystalline.prettyprint_irrep_matrix` calls
        matrix_sprint = (io, ir, i) -> Crystalline.prettyprint_irrep_matrix(io, ir, i)
        for i in eachindex(lgir.matrices)
            @test (sprint((io, ir) -> matrix_sprint(io, ir, i), lgir) == 
                   sprint((io, ir) -> matrix_sprint(io, ir, i), lgir′))
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