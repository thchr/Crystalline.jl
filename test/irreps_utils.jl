using Crystalline, Test, LinearAlgebra

@testset "Equality and hashing of irreps" begin
    # Julia defaults to the awkward ==(A,B) = A === B for user-defined types, which is not
    # which is usually impractical and awkward: we use StructEquality.jl to automatically
    # overload ==, to check field-based equality (same for hashing)

    # LGIrrep
    lgirsdA = lgirreps(221, Val(3))
    lgirsdB = lgirreps(221, Val(3))
    @test lgirsdA == lgirsdB
    @test lgirsdA !== lgirsdB # not `===`
    @test hash(lgirsdA) == hash(lgirsdB)

    # NewBandRep
    brsA = calc_bandreps(16, Val(3))
    brsB = calc_bandreps(16, Val(3))
    @test brsA == brsB
    @test brsA !== brsB # not `===`
    @test hash(brsA) == hash(brsB)

    # LGIrrep vectors from NewBandRep
    brsA_lgirsvA = irreps(brsA)
    brsB_lgirsvB = irreps(brsB)
    @test brsA_lgirsvA == brsB_lgirsvB
    @test brsA_lgirsvA !== brsB_lgirsvB # not `===`
    @test hash(brsA_lgirsvA) == hash(brsB_lgirsvB)

    # CompositeBandRep
    nA = @composite brsA[1]+brsA[2]
    nB = @composite brsB[1]+brsB[2]
    @test nA == nB
    @test nA !== nB # not `===`
    @test hash(nA) == hash(nB)

    # SymmetryVector
    _nA = SymmetryVector(nA)
    _nB = SymmetryVector(nB)
    @test _nA == _nB
    @test _nA !== _nB # not `===`
    @test hash(_nA) == hash(_nB)

    # PGIrrep
    pgirsA = pgirreps("m-3m", Val(3))
    pgirsB = pgirreps("m-3m", Val(3))
    @test pgirsA == pgirsB
    @test pgirsA !== pgirsB # not `===`
    @test hash(pgirsA) == hash(pgirsB)
    
    # SiteIrrep & SiteGroup
    sitegA = sitegroups(221, Val(3))[end]
    sitegB = sitegroups(221, Val(3))[end]
    @test sitegA == sitegB
    @test sitegA !== sitegB # not `===`
    @test hash(sitegA) == hash(sitegB)
    siteirsA = siteirreps(sitegA)
    siteirsB = siteirreps(sitegB)
    @test siteirsA == siteirsB
    @test siteirsA !== siteirsB # not `===`
    @test hash(siteirsA) == hash(siteirsB)
end