using Crystalline, Test
using Crystalline: _su2_entry_string

plain(x) = sprint(show, MIME"text/plain"(), x)

@testset "`show` overloads for double groups" begin

@testset "SU(2) entries" begin
    @test _su2_entry_string(0)                     == "0"
    @test _su2_entry_string(-1)                    == "-1"
    @test _su2_entry_string(-im)                   == "-i"
    @test _su2_entry_string(0.5)                   == "1/2"
    @test _su2_entry_string(im*sqrt(2)/2)          == "i√2/2"
    @test _su2_entry_string(-im*sqrt(3)/2)         == "-i√3/2"
    @test _su2_entry_string((1-im)/2)              == "(1-i)/2"
    @test _su2_entry_string(-(1+im)*sqrt(2)/2)     == "-(1+i)√2/2"
    @test _su2_entry_string(0.5 - im*sqrt(3)/2)    == "1/2-i√3/2"
    @test _su2_entry_string(-sqrt(3)/2 + 0.5im)    == "-√3/2+i/2"
    @test _su2_entry_string(cis(0.1))              == "0.995+0.0998i" # not crystallographic

    # every SU(2) element of every crystallographic operation is written exactly
    @test all(1:MAX_SGNUM[3]) do sgnum
        all(spacegroup(sgnum, Val(3); spinful=Val(true))) do dop
            all(z -> !occursin('.', _su2_entry_string(z)), matrix(su2(dop)))
        end
    end
end

@testset "SU2" begin
    u = su2(pointgroup("m-3m", Val(3); spinful=Val(true))[5])
    @test plain(u) == """
        SU(2) element:
         (1-i)/2  -(1+i)/2
         (1-i)/2   (1+i)/2"""
    @test sprint(show, SU2(-im, 0)) == "SU2(0.0 - 1.0im, 0.0 + 0.0im)"
end

@testset "DSymOperation" begin
    g = pointgroup("m-3m", Val(3); spinful=Val(true))
    @test plain(g[2]) == """
        2₀₀₁ ─────────────────────────── (-x,-y,z)
         ┌ -1  0  0 ╷ 0 ┐
         │  0 -1  0 ┆ 0 │  ┌ -i  0 ┐
         └  0  0  1 ╵ 0 ┘, └  0  i ┘"""
    @test plain(g[end]) == """
        ᵈm₁₀₋₁ ─────────────────────────── (z,y,x)
         ┌ 0  0  1 ╷ 0 ┐
         │ 0  1  0 ┆ 0 │  ┌ -i√2/2  i√2/2 ┐
         └ 1  0  0 ╵ 0 ┘, └  i√2/2  i√2/2 ┘"""

    # hexagonal setting, and a nonzero translation part
    h = pointgroup("6/mmm", Val(3); spinful=Val(true))
    @test plain(h[7]) == """
        2₁₁₀ ──────────────────────────── (y,x,-z)
         ┌ 0  1  0 ╷ 0 ┐
         │ 1  0  0 ┆ 0 │  ┌         0  -1/2-i√3/2 ┐
         └ 0  0 -1 ╵ 0 ┘, └ 1/2-i√3/2           0 ┘"""
    @test plain(spacegroup(230, Val(3); spinful=Val(true))[30]) == """
        {-3₁₋₁₁⁻|½,½,0} ───────── (-z+1/2,x+1/2,y)
         ┌ 0  0 -1 ╷ 1/2 ┐
         │ 1  0  0 ┆ 1/2 │  ┌ (1+i)/2  -(1-i)/2 ┐
         └ 0  1  0 ╵   0 ┘, └ (1+i)/2   (1-i)/2 ┘"""

    # compact forms, as used when printing groups
    @test sprint(show, MIME"text/plain"(), g[end]; context = :compact => true) == "ᵈm₁₀₋₁"
    @test sprint(show, g[end]) == "ᵈm₁₀₋₁"
    @test plain(g[[1, 49]]) == """
        2-element Vector{DSymOperation{3}}:
         1
         ᵈ1"""
end

@testset "Double little and site groups: position labels" begin
    lg = littlegroups(221, Val(3); spinful=Val(true))["X"]
    @test startswith(sprint(show, MIME"text/plain"(), lg),
                     "DLittleGroup{3} ⋕221 (Pm-3m) at X = [0, 1/2, 0] with 32 operations:")
    wp = only(filter(wp -> label(wp) == "2p", wyckoffs(47, Val(3))))
    siteg = sitegroup(spacegroup(47, Val(3); spinful=Val(true)), wp)
    @test startswith(sprint(show, MIME"text/plain"(), siteg),
                     "DSiteGroup{3} ⋕47 (Pmmm) at 2p = [1/2, β, 1/2] with 8 operations:")
end

end # @testset "`show` overloads for double groups"
