using Crystalline, Test

if !isdefined(Main, :LGIRSDIM)
    LGIRSDIM = Tuple(lgirreps.(1:MAX_SGNUM[D], Val(D)) for D in 1:3)
end

@testset "Little group order" begin
    # see e.g. Bradley & Cracknell p. 151(bottom)-152(top)
    @testset "Decomposition, order[star{k}]*order[G₀ᵏ] = order[G₀]" begin
    for (D, LGIRS) in enumerate(LGIRSDIM)
        for (sgnum, lgirsd) in enumerate(LGIRS)
            cntr = centering(sgnum, D)  
            # the "macroscopic" order is defined simply as the length of the 
            # point group associated with the space group
            sgops = operations(spacegroup(sgnum, D)) # from crawling Bilbao
            order_macroscopic = length(pointgroup(sgops))
            
            for lgirs in values(lgirsd)
                kv = kvec(first(lgirs))    
                for lgir in lgirs            
                    # number of k-vectors in the star of k
                    order_kstar = length(kstar(sgops, kv, cntr)) 
                    # number of operations in the little group of k
                    order_pointgroupofk = length(pointgroup(operations(lgir)))

                    # test that q⋅b=h, where
                    #   q ≡ order[star{k}]
                    #   b ≡ order[G₀ᵏ]      ≡ order of little co-group of k
                    #   h ≡ order[G₀]       ≡ macroscopic order of G
                    # and where G₀ denotes the point group derived from the space 
                    # group G, Gᵏ denotes the little group of k derived from G, 
                    # and G₀ᵏ denotes the point group derived from Gᵏ
                    @test order_kstar*order_pointgroupofk == order_macroscopic
                end
            end
        end
    end # for (D, LGIRS) in enumerate(LGIRSDIM)
    end # @testset "Decomposition, ..."

    @testset "Macroscopic order, Bilbao vs. ISOTROPY" begin
    for (D, LGIRS) in enumerate(LGIRSDIM)
        for (sgnum, lgirsd) in enumerate(LGIRS)
            sgops_bilbao = operations(spacegroup(sgnum, D)) # from crawling Bilbao
            sgops_isotropy = operations(first(lgirsd["Γ"])) # from ISOTROPY's Γ-point LittleGroup

            # note that we don't have to worry about whether the operations 
            # are represented in conventional or primitive bases, and whether 
            # there are "repeated translation sets" in the former; taking the 
            # point group automatically removes such sets
            order_macroscopic_bilbao = length(pointgroup(sgops_bilbao))
            order_macroscopic_isotropy = length(pointgroup(sgops_isotropy))

            @test order_macroscopic_bilbao == order_macroscopic_isotropy
        end
    end # for (D, LGIRS) in enumerate(LGIRSDIM)
    end # @testset "Macroscopic order, ..."
end

@testset "Little group data-consistency" begin
    for D in 1:3
        for sgnum in 1:MAX_SGNUM[D]
            (D == 2 && !issymmorph(sgnum, D)) && continue # 2D nonsymmorphic SGs not yet implemented
            lgs = littlegroups(sgnum, D)
            sg  = spacegroup(sgnum, Val(D))
            sg  = SpaceGroup{D}(sgnum, reduce_ops(sg))
            for (klab, lg) in lgs
                kv  = kvec(lg)
                lg′ = littlegroup(sg, kv)

                # we compare the primitive little groups; otherwise there can be spurious
                # differences relating to nonsymmorphic translations that "look" different
                # in the conventional basis but in fact are equal in the primitive basis.
                # an example is {m₀₁₀|0,½,½} (or S"x,-y+1/2,z+1/2") vs. {m₀₁₀|½,0,0} (or
                # "x+1/2,-y,z") in an 'I' centered space group; both are actually included
                # in `sg`, but when we reduce it, it may choose one or the other
                @test Set(xyzt.(primitivize(lg))) == Set(xyzt.(primitivize(lg′)))
            end
        end
    end
end