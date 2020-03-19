using SGOps, Test

@testset "Find a point group for every space group" begin
    for D = 1:3
        for sgnum in 1:MAX_SGNUM[D]
            G = spacegroup(sgnum)
            pg = SGOps.find_parent_pointgroup(G)
            @test pg !== nothing
        end
    end
end