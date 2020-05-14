using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS  = parselittlegroupirreps() # parsed directly from ISOTROPY's files
end
if !isdefined(Main, :LGIRS′)
    LGIRS′ = get_all_lgirreps(Val(3))  # loaded from our saved .jld2 files
end

@testset "Test equivalence of parsed and loaded LGIrreps" begin
    for sgnum in 1:230
        lgirsvec  = LGIRS[sgnum]  # parsed variant
        lgirsvec′ = LGIRS′[sgnum] # loaded variant

        @test length(lgirsvec) == length(lgirsvec′)
        for (kidx, lgirs) in enumerate(lgirsvec)
            lgirs′ = lgirsvec′[kidx]
            @test length(lgirs) == length(lgirs′)
            for (iridx, lgir) in enumerate(lgirs)
                lgir′ = lgirs′[iridx]
                # test that labels agree
                @test label(lgir) == label(lgir′)
                # test that little groups agree
                @test isapprox(kvec(lgir), kvec(lgir′))
                @test all(operations(lgir) .== operations(lgir′))
                # test that irreps agree
                for kabc in (nothing, Crystalline.TEST_αβγ)
                    @test irreps(lgir) == irreps(lgir′)
                end

            end
        end
    end
end