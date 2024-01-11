using Crystalline, Test
@testset "`spacegroup` vs `mspacegroup` agreement" begin
    for sgnum in 1:MAX_SGNUM[3]
        bns_num = Crystalline.SG2MSG_NUMs_D[sgnum]

        sg = spacegroup(sgnum)
        msg = Crystalline.mspacegroup(bns_num[1], bns_num[2])

        sort_sg  = sort(sg,  by=seitz)
        sort_msg = sort(SymOperation.(msg), by=seitz)

        @test sort_sg ≈ sort_msg
    end
end