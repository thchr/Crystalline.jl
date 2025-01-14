using Crystalline, Test

@testset "subgroup traversal and transformations" begin
    # traversal and transformations via subgroups (G<H)
    for sgnumᴳ in 1:MAX_SGNUM[3]
        sgᴳ = spacegroup(sgnumᴳ)
        gr = maximal_subgroups(sgnumᴳ)
        for sgnumᴴ in gr.nums[2:end]
            classes = conjugacy_relations(gr, sgnumᴳ, sgnumᴴ)
            sgᴴ = reduce_ops(spacegroup(sgnumᴴ), centering(sgnumᴴ))
            for (i,t) in enumerate(classes)
                P, p = t.P, t.p
                sgᴳ′ = reduce_ops(transform.(sgᴳ, Ref(P), Ref(p)), centering(sgnumᴴ))
                @test issubgroup(sgᴳ′, sgᴴ)
            end
        end
    end
end # @testset

@testset "supergroup traversal and transformations" begin
    # traversal and transformations through supergroups (H<G)
    for sgnumᴴ in 1:MAX_SGNUM[3]
        sgᴴ = reduce_ops(spacegroup(sgnumᴴ), centering(sgnumᴴ))
        gr = minimal_supergroups(sgnumᴴ)
        for sgnumᴳ in gr.nums[2:end]
            classes = conjugacy_relations(gr, sgnumᴴ, sgnumᴳ)
            sgᴳ = spacegroup(sgnumᴳ)
            for (i,t) in enumerate(classes)
                P, p = t.P, t.p
                sgᴳ′ = reduce_ops(transform.(sgᴳ, Ref(P), Ref(p)), centering(sgnumᴴ))
                @test issubgroup(sgᴳ′, sgᴴ)
            end
        end
    end
end # @testset