@testset "Miscellaneous" begin
    @testset "`Collect` of non-Vector input" begin
        vs = [1,2,3]
        vs_range = 1:3
        c = Collection(vs)

        @test c == Collection(vs_range)
        @test c == Collection(c)
        @test c === Collection(c) # returned as-is (egal, ===), no copy of anything

        lgirs = lgirreps(22)["Γ"]
        lgirs_vs = lgirs[2:end-1]
        lgirs_view = @view lgirs[2:end-1] # a view type
        lgirs_vs_subset = Collection(lgirs_vs)
        @test lgirs_vs_subset == Collection(lgirs_view) # collect the `view`
        @test lgirs_vs_subset === Collection(lgirs_vs_subset) # returned as-is (egal, ===)
    end

    @testset "`uniquetol`: specialized on positional arguments" begin
        # positional arguments forwarded to `isapprox` (here, `cntr` & `modw`) must not make
        # each comparison a dynamic dispatch (previously, ~2 allocations per comparison)
        ops = operations(spacegroup(230, Val(3)))
        count_allocs(ops) = @allocations Crystalline.uniquetol(ops, nothing, false)
        count_allocs(ops) # compile
        @test count_allocs(ops) < length(ops) # only growth of the returned vector
        @test Crystalline.uniquetol(ops, nothing, false) == ops # ops are distinct
        @test Crystalline.isapproxin(ops[end], ops[1:end-1], nothing, false) == false
    end
end