using Test, Crystalline, StaticArrays

@testset "Compatibility" begin
@testset "Subduction" begin
    lgirsd = lgirreps(218)
    Rᵢ²¹⁸ = realify(lgirsd["R"]) # 2D, 4D, & 6D irreps R₁R₂, R₃R₃, & R₄R₅ at R = [½,½,½]
    Λᵢ²¹⁸ = realify(lgirsd["Λ"]) # 1D, 1D, & 2D irreps Λ₁, Λ₂, Λ₃ at Λ = [α,α,α]

    # --- point to line subduction (R → Λ) ---
    # R₁R₂ → Λ₁+Λ₂ | R₃R₃ → 2Λ₃ | R₄R₅  → Λ₁+Λ₂+2Λ₃
    @test [[subduction_count(Rᵢ, Λᵢ) for Λᵢ in Λᵢ²¹⁸] for Rᵢ in Rᵢ²¹⁸] == [[1,1,0],[0,0,2],[1,1,2]]

    # --- corep-specifics ---
    R₄R₅²¹⁸ = Rᵢ²¹⁸[end]
    # test that complex corep subduces to itself
    @test subduction_count(R₄R₅²¹⁸, R₄R₅²¹⁸) == 1

    # test that physically real corep subduces to its complex components
    Rᵢ²¹⁸_notr = lgirsd["R"]
    @test [subduction_count(R₄R₅²¹⁸, Rᵢ_notr) for Rᵢ_notr in Rᵢ²¹⁸_notr] == [0,0,0,1,1]

    # --- sub-groups ---
    # test that two distinct coreps in different groups decompose consistently
    R₄²²² = realify(lgirreps(222)["R"])[end]; # 6D irrep
    R₄R₅²¹⁸′ = deepcopy(R₄R₅²¹⁸)              # 6D irrep
    
    # SG 218 is a subgroup of 222, but in a different setting - so convert the little group
    # of 218 to that of 222
    lg_ops = operations(group(R₄R₅²¹⁸′)) 
    lg_ops .= transform.(lg_ops, Ref(one(SMatrix{3,3})), Ref((@SVector [1,1,1])/4))

    # check that R₄²²² subduces into R₄R₅²¹⁸
    @test subduction_count(R₄²²², R₄R₅²¹⁸′) == 1
end # @testset "Subduction"

@testset "`remap_to_kstar`" begin
    # TODO: very limited, incomplete testing - testing only a single bug encountered. This
    #       needs actual testing of the implementation in cases where it does something
    #       nontrivial

    # trivial case: input k-point is the same as the actual k-point of the `lgirs` - nothing
    # needs to be done, and the same exact set of irreps should just be returned (deepcopy,
    # not aliased at any field or indexing level)
    lgirs = lgirreps(225)["X"]
    kv = position(lgirs)
    lgirs_map_to_self = remap_to_kstar(lgirs, kv)
    @test lgirs_map_to_self ==  lgirs  # map to same k-star doesn't do anything
    @test lgirs_map_to_self !== lgirs # should not alias
    @test lgirs_map_to_self[1] ==  lgirs[1]
    @test lgirs_map_to_self[1] !== lgirs[1]
    
    # smoke-test: test it works/runs without error for non-trivial (different) `kv` input
    kv′ = kstar(first(lgirs))[2]
    @test kv ≠ kv′
    @test remap_to_kstar(lgirs, kv′) isa Collection{LGIrrep{3}}
end
end # @testset "Compatibility"