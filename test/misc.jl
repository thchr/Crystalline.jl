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
end