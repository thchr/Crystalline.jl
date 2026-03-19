using Test
using Crystalline

@testset "isspecial" begin
    # special points
    @test isspecial(KVec("1/2,1/2,1/2")) == true
    @test isspecial(RVec("1/2,0,1/2")) == true
    @test isspecial(RVec("1/2,1/3")) == true

    # non-special points
    @test isspecial(KVec("1/2,α,1/2")) == false
    @test isspecial(RVec("1/2,α,1/2+γ")) == false
    @test isspecial(RVec("1/2,1/3+β")) == false

    # WyckoffPosition
    wps = wyckoffs(22)
    wp_4a  = wps[end] # 4a  [0, 0, 0] (special)
    wp_16k = wps[1]   # 16k [α, β, γ] (not special)
    @test isspecial(wp_4a) == true
    @test isspecial(wp_16k) == false

    # SiteGroup
    siteg_1a = sitegroup(22, wp_4a)
    siteg_16k = sitegroup(22, wp_16k)
    @test isspecial(siteg_1a) == true
    @test isspecial(siteg_16k) == false

    # SiteIrrep and Collection{<:SiteIrrep}
    @test isspecial(siteirreps(siteg_1a)) == true   # testing `Collection{<:SiteIrrep}`
    @test isspecial(siteirreps(siteg_16k)) == false # --||--
    @test all(isspecial, siteirreps(siteg_1a)) == true # testing ::SiteIrrep directly
    @test all(isspecial, siteirreps(siteg_16k)) == false

    # LittleGroup, Collection{<:LGIrrep}, LGIrrep
    lgs = littlegroups(230, Val(3))
    @test isspecial(lgs["Γ"]) == true
    @test isspecial(lgs["Ω"]) == false
    @test isspecial(lgs["Λ"]) == false

    lgirsd = lgirreps(230, Val(3))
    @test isspecial(lgirsd["Γ"]) == isspecial(first(lgirsd["Γ"])) == true
    @test isspecial(lgirsd["Σ"]) == isspecial(first(lgirsd["Σ"])) == false

    # Errors
    # ⇒ error("`nothing` is not a valid position input to `isspecial` [...]")
    @test_throws ErrorException isspecial(nothing)
    @test_throws ErrorException isspecial(spacegroup(22, Val(3)))
    @test_throws ErrorException isspecial(pgirreps("3", Val(3)))
end