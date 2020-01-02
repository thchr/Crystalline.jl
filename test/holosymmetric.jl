using SGOps, Test


@testset "Representation vs basic domain BZ" begin

@testset "Identification of holosymmetric space groups" begin
    # Check that our cached values of the holosymmetric space group numbers
    # (stored in SGOps.HOLOSYMMETRIC_SGNUMS) agree with the results of 
    # SGOps._find_holosymmetric_sgnums(D)
    recalc_holosgnums = Tuple(tuple(SGOps._find_holosymmetric_sgnums(D)...) for D = 1:3)
    @test recalc_holosgnums == SGOps.HOLOSYMMETRIC_SGNUMS
    
    # Compare our calculation of the holosymmetric space groups with Table 3.10 
    # of CDML: that table lists the `sgnum`s of every _symmorphic_ nonholosymmetric
    # group in 3D, allowing us a sanity check on _find_holosymmetric_sgnums(D::Integer).
    holo_sgnums    = SGOps._find_holosymmetric_sgnums(3)
    nonholo_sgnums = filter(sgnum->sgnum∉holo_sgnums, 1:MAX_SGNUM[3])
    symmorph_nonholo_sgnums = filter(issymmorph, nonholo_sgnums)

    symmorph_nonholo_sgnums_CDML = (
        1,3,6,5,8,16,25,21,35,38,22,42,23,44,75,81,83,89,99,111,115,79,82,87,
        97,107,119,121,195,200,207,215,196,202,209,216,197,204,211,217,146,
        148,155,160,143,147,150,157,168,149,156,162,164,177,183,174,175,189,187)

    @test symmorph_nonholo_sgnums == sort([symmorph_nonholo_sgnums_CDML...])

    # Check that our cached values for ISOGONAL_PARENT_GROUPS remain computed correctly:
    isogonal_parent_groups′ = Tuple(tuple(getindex.(SGOps.get_arith_crystalclass_partner.(1:MAX_SGNUM[D], D), 2)...) for D in 1:3)
    @test SGOps.ISOGONAL_PARENT_GROUPS == isogonal_parent_groups′
end


# Test the corner cases noted for the method `_find_holosymmetric_sgnums`
@testset "Cornercases and transformations for _find_holosymmetric_sgnums" begin
    for (sgnum, info) in zip(keys(SGOps.CORNERCASES_SUBSUPER_NORMAL_SGS), 
                             values(SGOps.CORNERCASES_SUBSUPER_NORMAL_SGS))
        sgnum₀ = info[1] # supergroup number
        P, p = info[2:3] # transformation pair (P,p)
        cntr = centering(sgnum, 3)

        G  = operations(get_sgops(sgnum,  3))  # G in conventional basis of sgnum
        Gᵖ = reduce_ops(G, cntr, false)        # G in primitive basis of sgnum
        
        G₀ = operations(get_sgops(sgnum₀, 3))  # G₀ in conventional basis of sgnum₀    
        G₀′ = transform.(G₀, Ref(P), Ref(p))   # G₀ in conventional basis of sgnum
        G₀′ᵖ = reduce_ops(G₀′, cntr, false)    # G₀ in primitive basis of sgnum

        # verify that the G₀′ is indeed a normal, holosymmetric supergroup of G
        @test SGOps.issubgroup(G₀′ᵖ, Gᵖ)
        @test SGOps.isnormal(G₀′ᵖ, Gᵖ)
        @test SGOps.is_holosymmetric(sgnum₀, 3)
    end
end


@testset "Finding holosymmetric, normal supergroups (\"parent group\")" begin
    holosym_parents = SGOps.find_holosymmetric_parent.(1:MAX_SGNUM[3],3)
    
    # We know that there should be 24 orphans, as tabulated in CDML.
    # Verify that these are the cases "missed" by find_holosymmetric_parent
    holosym_parents_orphans = findall(isnothing, holosym_parents)
    @test holosym_parents_orphans == sort(collect(Iterators.flatten(SGOps.ORPHAN_SGNUMS)))
end


@testset "Test find_holosymmetric_superpointgroup" begin
    # the minimal super point group P₀ of a holosymmetric space group G₀ should equal its isogonal point group F₀
    for D in 1:3
        for sgnum₀ in SGOps.HOLOSYMMETRIC_SGNUMS[D]
            G₀ = get_sgops(sgnum₀, D)
            F₀ = pointgroup(G₀)                                            # isogonal point group of G₀
            P₀ = operations(SGOps.find_holosymmetric_superpointgroup(G₀))  # minimal holosymmetric super point group of G₀
            @test sort(xyzt.(F₀)) == sort(xyzt.(P₀))
        end
    end

    # every space group should have a holosymmetric super point group (i.e. never return nothing)
    for D = 1:3
        @test !any(sgnum->isnothing(SGOps.find_holosymmetric_superpointgroup(sgnum, D)), 1:MAX_SGNUM[D])
    end
end

@testset "Every space group maps to a holosymmetric group or to a group included in Φ-Ω?" begin
    D = 3
    # only need to do anything if sg is nonholosymmetric
    nonholo_sgnums = filter(sgnum->!SGOps.is_holosymmetric(sgnum,D), 1:MAX_SGNUM[D])
    # space group numbers of isogonal space groups that must exist as keys in SGOps.ΦNOTΩ_KVECS_AND_MAPS
    needed_sgnums_in_ΦnotΩ = unique(getindex.(Ref(SGOps.ISOGONAL_PARENT_GROUPS[D]), nonholo_sgnums))
    # actual keys in SGOps.ΦNOTΩ_KVECS_AND_MAPS
    sgsums_in_ΦnotΩ = keys(SGOps.ΦNOTΩ_KVECS_AND_MAPS)
    for needed_sgnum in needed_sgnums_in_ΦnotΩ
        if needed_sgnum ∉ sgsums_in_ΦnotΩ
            @test_broken needed_sgnum ∈ sgsums_in_ΦnotΩ
        end
        # TODO: This should not actually be an error! Some sgs just don't get "new" k-points
        # in Φ-Ω since they already exist in kstar or are equivalent via a primitive G-vec. 
        # We could verify that that is the case for the 6 "errors" here.
    end
end

end # file @testset