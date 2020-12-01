using Crystalline, Test, LinearAlgebra

# TODO: Cannot pass unless src/special_representation_domain_kpoints.jl is reincluded in 
#       src/Crystalline.jl. Currently disabled due to compilation time concerns.

@testset "Representation vs basic domain BZ" begin

# disable if src/special_representation_domain_kpoints.jl is reincluded
if (@test_broken false) isa Test.Broken
    println("Skipping tests of src/special_representation_domain_kpoints.jl")

else

@testset "Identification of holosymmetric space groups" begin
    # Check that our cached values of the holosymmetric space group numbers
    # (stored in Crystalline.HOLOSYMMETRIC_SGNUMS) agree with the results of 
    # Crystalline._find_holosymmetric_sgnums(D)
    recalc_holosgnums = Tuple(tuple(Crystalline._find_holosymmetric_sgnums(D)...) for D = 1:3)
    @test recalc_holosgnums == Crystalline.HOLOSYMMETRIC_SGNUMS
    
    # Compare our calculation of the holosymmetric space groups with Table 3.10 
    # of CDML: that table lists the `sgnum`s of every _symmorphic_ nonholosymmetric
    # group in 3D, allowing us a sanity check on _find_holosymmetric_sgnums(D::Integer).
    holo_sgnums    = Crystalline._find_holosymmetric_sgnums(3)
    nonholo_sgnums = filter(sgnum->sgnum∉holo_sgnums, 1:MAX_SGNUM[3])
    symmorph_nonholo_sgnums = filter(issymmorph, nonholo_sgnums)

    symmorph_nonholo_sgnums_CDML = (
        1,3,6,5,8,16,25,21,35,38,22,42,23,44,75,81,83,89,99,111,115,79,82,87,
        97,107,119,121,195,200,207,215,196,202,209,216,197,204,211,217,146,
        148,155,160,143,147,150,157,168,149,156,162,164,177,183,174,175,189,187)

    @test symmorph_nonholo_sgnums == sort([symmorph_nonholo_sgnums_CDML...])

    # Check that our cached values for ARITH_PARTNER_GROUPS remain computed correctly: (avoid bit-rotting)
    arith_partner_groups′ = Tuple(tuple(getindex.(Crystalline._find_arithmetic_partner.(1:MAX_SGNUM[D], D), 2)...) for D in 1:3)
    @test Crystalline.ARITH_PARTNER_GROUPS == arith_partner_groups′
end


# Test the corner cases noted for the method `_find_holosymmetric_sgnums`
@testset "Cornercases and transformations for _find_holosymmetric_sgnums" begin
    for (sgnum, info) in zip(keys(Crystalline.CORNERCASES_SUBSUPER_NORMAL_SGS), 
                             values(Crystalline.CORNERCASES_SUBSUPER_NORMAL_SGS))
        sgnum₀ = info[1] # supergroup number
        P, p = info[2:3] # transformation pair (P,p)
        cntr = centering(sgnum, 3)

        G  = operations(spacegroup(sgnum,  3))  # G in conventional basis of sgnum
        Gᵖ = reduce_ops(G, cntr, false)        # G in primitive basis of sgnum
        
        G₀ = operations(spacegroup(sgnum₀, 3))  # G₀ in conventional basis of sgnum₀    
        G₀′ = transform.(G₀, Ref(P), Ref(p))   # G₀ in conventional basis of sgnum
        G₀′ᵖ = reduce_ops(G₀′, cntr, false)    # G₀ in primitive basis of sgnum

        # verify that the G₀′ is indeed a normal, holosymmetric supergroup of G
        @test issubgroup(G₀′ᵖ, Gᵖ) && issubgroup(G₀′, G) # test for both conventional 
        @test isnormal(G₀′ᵖ, Gᵖ)   && isnormal(G₀′, G)   # **and** primitive bases
        @test Crystalline.is_holosymmetric(sgnum₀, 3)
    end
end


@testset "Finding holosymmetric, normal supergroups (\"parent group\")" begin
    holosym_parents = Crystalline.find_holosymmetric_parent.(1:MAX_SGNUM[3],3)
    
    # We know that there should be 24 orphans, as tabulated in CDML.
    # Verify that these are the cases "missed" by find_holosymmetric_parent
    holosym_parents_orphans = findall(isnothing, holosym_parents)
    @test holosym_parents_orphans == sort(collect(Iterators.flatten(Crystalline.ORPHAN_SGNUMS)))
end


@testset "Test find_holosymmetric_superpointgroup" begin
    # the minimal super point group P₀ of a holosymmetric space group G₀ should equal its isogonal/arithmetic point group F₀
    for D in 1:3
        for sgnum₀ in Crystalline.HOLOSYMMETRIC_SGNUMS[D]
            G₀ = spacegroup(sgnum₀, D)
            F₀ = pointgroup(G₀)                                            # isogonal/arithmetic point group of G₀
            P₀ = operations(Crystalline.find_holosymmetric_superpointgroup(G₀))  # minimal holosymmetric super point group of G₀
            @test sort(xyzt.(F₀)) == sort(xyzt.(P₀))
        end
    end

    # every space group should have a holosymmetric super point group (i.e. never return nothing)
    for D = 1:3
        @test !any(sgnum->isnothing(Crystalline.find_holosymmetric_superpointgroup(sgnum, D)), 1:MAX_SGNUM[D])
    end
end

@testset "Compare find_new_kvecs(sgnum, D) to \"new\" kvecs from CDML" begin
    D = 3
    # check that all holosymmetric sgs get no new k-vecs (trivial)
    @test all(sgnum-> Crystalline.find_new_kvecs(sgnum, 3) === nothing, Crystalline.HOLOSYMMETRIC_SGNUMS[D])

    # check that the new k-vecs from find_new_kvecs(...) agree with those of CDML 
    # for the case of nonholosymmetric sgs [TEST_BROKEN]
    verbose_debug = false # toggle to inspect disagreements between CDML and find_new_kvecs
    if verbose_debug; more, fewer = Int64[], Int64[]; end

    for sgnum in Base.OneTo(MAX_SGNUM[D])
        Crystalline.is_holosymmetric(sgnum, D) && continue # only check holosymmetric sgs

        # look for new kvecs (and mapping) in CDML data; if nothing is there, return an empty 
        # array of Crystalline.KVecMapping
        arith_partner_sgnum = Crystalline.find_arithmetic_partner(sgnum, D)
        kvmaps_CDML = get(Crystalline.ΦNOTΩ_KVECS_AND_MAPS, arith_partner_sgnum, Crystalline.KVecMapping[])
        newklabs_CDML = [getfield.(kvmaps_CDML, :kᴮlab)...]

        # look for new kvecs using our own search implementation
        _, _, _, newklabs = Crystalline.find_new_kvecs(sgnum, D)
        newklabs = collect(Iterators.flatten(newklabs)) # flatten array of array structure
        newklabs .= rstrip.(newklabs, '′') # strip "′" from newklabs

        # compare the set of labels
        sort!(newklabs_CDML)
        sort!(newklabs)
        @test_skip newklabs == newklabs_CDML # this is broken at present: toggle verbose_debug to inspect issue more clearly
        # TODO: Make these tests pass (low priority; find_new_kvecs is not used anywhere)

        if verbose_debug && (newklabs != newklabs_CDML)
            bt = bravaistype(sgnum, 3)
            if length(newklabs_CDML) > length(newklabs)
                @info "More KVecs in CDML:" sgnum bt setdiff(newklabs_CDML, newklabs)
                push!(more, sgnum)
            elseif length(newklabs) > length(newklabs_CDML)
                @info "Fewer KVecs in CDML:" sgnum bt setdiff(newklabs, newklabs_CDML)
                push!(fewer, sgnum)
            else
                @info "Distinct KVecs from CDML:" sgnum bt newklabs newklabs_CDML
            end
        end       
    end
    if verbose_debug && (!isempty(fewer) || !isempty(more))
        println("\n$(length(more)+length(fewer)) disagreements were found relative to CDML's listing of new KVecs:")
        print("   More KVecs in CDML:\n      ");  join(stdout, more,  ' '); println()
        print("   Fewer KVecs in CDML:\n      "); join(stdout, fewer, ' '); println()
    end

end

@testset "Tabulated orphan parent groups are indeed normal supergroups" begin
    Pᴵ = Matrix{Float64}(I, 3, 3) # trivial transformation rotation (only translation may be nontrivial here)
    for (sgnum, (parent_sgnum, p)) in pairs(Crystalline.ORPHAN_AB_SUPERPARENT_SGNUMS)
        G  = operations(spacegroup(sgnum, 3))
        G′ = operations(spacegroup(parent_sgnum, 3)) # basis setting may differ from G (sgs 151-154)

        # G′ in conventional basis of G
        if !iszero(p)
            G′ .= transform.(G′, Ref(Pᴵ), Ref(p))
            # Note: we don't have to bother with going to the primitive lattice since
            #       all the group/supergroup pairs have the same bravais type.
        end

        # test that G◁G′
        @test issubgroup(G′, G)
        @test isnormal(G′, G)
    end
end

end # @test_broken if
end # outer @testset