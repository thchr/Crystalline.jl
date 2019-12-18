using SGOps, Test

@testset "Find a point group for every space group" begin
    for D = 1:3
        for sgnum in 1:MAX_SGNUM[D]
            G = get_sgops(sgnum)
            pg_iuclab = SGOps.find_parent_pointgroup(G)
            @test pg_iuclab !== nothing
        end
    end
end