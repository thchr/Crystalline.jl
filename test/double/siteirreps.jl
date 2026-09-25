using Crystalline, Test, LinearAlgebra
using Crystalline: check_multtable_vs_ir

@testset "Double site groups and their irreps" begin
    for sgnum in 1:MAX_SGNUM[3]
        sg = spacegroup(sgnum, Val(3))
        dsg = spacegroup(sgnum, Val(3); spinful=Val(true))
        for wp in wyckoffs(sgnum, Val(3))
            siteg = sitegroup(sg, wp)
            dsiteg = doublegroup(siteg)
            dsiteg′ = sitegroup(dsg, wp)
            @test operations(dsiteg′) == operations(dsiteg)
            @test cosets(dsiteg′) == cosets(dsiteg)
            n = length(siteg)
            @test length(dsiteg) == 2n
            @test position(dsiteg) == position(siteg)
            @test length(cosets(dsiteg)) == length(cosets(siteg))

            siteirs = siteirreps(dsiteg)
            @test siteirs isa Collection{DSiteIrrep{3}}
            @test all(isspinful, siteirs)
            @test sum(ir -> irdim(ir)^2, siteirs) == n # double-valued irreps only

            mt = MultTable(dsiteg)
            χs = [characters(ir) for ir in siteirs]
            for (a, ir) in enumerate(siteirs)
                @test all(check_multtable_vs_ir(mt, ir))
                @test ir.matrices[n+1:2n] ≈ -ir.matrices[1:n] # D(Ēg) = -D(g)

                # orthogonality of characters, over the 2n operations of the double group
                for b in eachindex(siteirs)
                    @test dot(χs[b], χs[a]) / 2n ≈ (a == b) atol=1e-10
                end

                # reality, by the Frobenius-Schur indicator
                fs = sum(χs[a][mt.table[i,i]] for i in 1:2n) / 2n
                @test round(Int, real(fs)) == Int(reality(ir))
            end
        end
    end
end

# the site irrep labels are compared with Bilbao's in test/bandreps.jl (spinful EBRs)

