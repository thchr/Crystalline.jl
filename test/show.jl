using Crystalline, Test
using LinearAlgebra: dot

# ---------------------------------------------------------------------------------------- #
# test print with nicely printed diff on failures (from https://github.com/invenia/PkgTemplates.jl/blob/master/test/runtests.jl)
using DeepDiffs: deepdiff
function print_diff(a, b)
    old = Base.have_color
    @eval Base have_color = true
    try
        println(deepdiff(a, b))
    finally
        @eval Base have_color = $old
    end
end

function test_show(expected::AbstractString, observed::AbstractString)
    if expected == observed
        @test true
    else
        print_diff(expected, observed)
        @test :expected == :observed
    end
end
test_tp_show(v, observed::AbstractString) = test_show(repr(MIME"text/plain"(), v), observed)

# ---------------------------------------------------------------------------------------- #

@testset "`show` overloads" begin
# -------------------------------
# DirectBasis
# -------------------------------
Rs = DirectBasis([1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]) # cubic
str = """
      DirectBasis{3} (cubic):
       [1.0, 0.0, 0.0]
       [0.0, 1.0, 0.0]
       [0.0, 0.0, 1.0]"""
test_tp_show(Rs, str)

Rs = DirectBasis([1,0,0], [0,1,0], [0,0,1]) # cubic with Int coordinate value type
str = """
      DirectBasis{3, Int64} (cubic):
       [1, 0, 0]
       [0, 1, 0]
       [0, 0, 1]"""
test_tp_show(Rs, str)

Rs = DirectBasis([1,0,0], [-0.5, √(3)/2, 0.0], [0, 0, 1.5]) # hexagonal
str = """
      DirectBasis{3} (hexagonal):
       [1.0, 0.0, 0.0]
       [-0.5, 0.8660254037844386, 0.0]
       [0.0, 0.0, 1.5]"""
test_tp_show(Rs, str)
Rs′ = directbasis(183, Val(3))
@test Rs[1] ≈ Rs′[1] && Rs[2] ≈ Rs′[2]
@test abs(dot(Rs[1], Rs′[3])) < 1e-14
@test abs(dot(Rs[2], Rs′[3])) < 1e-14
@test abs(dot(Rs[3], Rs′[3])) > 1e-1

Gs = dualbasis(Rs)
str = """
      ReciprocalBasis{3} (hexagonal):
       [6.283185307179586, 3.6275987284684357, -0.0]
       [0.0, 7.255197456936871, 0.0]
       [0.0, -0.0, 4.1887902047863905]"""
test_tp_show(Gs, str)

# test printing for unitful bases; don't clobber with big type before "[...]"
using Unitful
Rs = DirectBasis([1.0, 0, 0]*1u"Å" , [-0.5, sqrt(3)/2, 0]*1u"Å" ,   [0, 0, 1.25]*1u"Å")
str = """
      DirectBasis{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}} (hexagonal):
       [1.0 Å, 0.0 Å, 0.0 Å]
       [-0.5 Å, 0.8660254037844386 Å, 0.0 Å]
       [0.0 Å, 0.0 Å, 1.25 Å]"""
test_tp_show(Rs, str)

if !Sys.isapple()
    # don't test on Apple, since Unitful prints "^-1" as "⁻¹" there
    # (see https://github.com/PainterQubits/Unitful.jl/pull/446#issuecomment-2770636653)
    Gs = dualbasis(Rs)
    str = """
          ReciprocalBasis{3, Quantity{Float64, 𝐋^-1, Unitful.FreeUnits{(Å^-1,), 𝐋^-1, nothing}}} (hexagonal):
          [6.283185307179586 Å^-1, 3.6275987284684357 Å^-1, -0.0 Å^-1]
          [0.0 Å^-1, 7.255197456936871 Å^-1, 0.0 Å^-1]
          [0.0 Å^-1, -0.0 Å^-1, 5.026548245743669 Å^-1]"""
    test_tp_show(Gs, str)
end

# -------------------------------
# SymOperation
# -------------------------------
str = """
      1 ──────────────────────────────── (x,y,z)
       ┌ 1  0  0 ╷ 0 ┐
       │ 0  1  0 ┆ 0 │
       └ 0  0  1 ╵ 0 ┘"""
test_tp_show(S"x,y,z", str)

str = """
      {-3₋₁₋₁₁⁺|0,½,⅓} ──────── (z,-x+1/2,y+1/3)
       ┌  0  0  1 ╷   0 ┐
       │ -1  0  0 ┆ 1/2 │
       └  0  1  0 ╵ 1/3 ┘"""
test_tp_show(S"z,-x+1/2,y+1/3", str)

str = """
      3⁻ ───────────────────────────── (-x+y,-x)
       ┌ -1  1 ╷ 0 ┐
       └ -1  0 ╵ 0 ┘"""
test_tp_show(S"y-x,-x", str)

str = """
      3-element Vector{SymOperation{3}}:
       1
       2₀₁₁
       {3₁₁₁⁻|0,0,⅓}"""
test_tp_show([S"x,y,z", S"-x,z,y", S"y,z,x+1/3"], str)

str = """
      4×4 Matrix{SymOperation{3}}:
       1     2₀₀₁  2₀₁₀  2₁₀₀
       2₀₀₁  1     2₁₀₀  2₀₁₀
       2₀₁₀  2₁₀₀  1     2₀₀₁
       2₁₀₀  2₀₁₀  2₀₀₁  1"""
sg = spacegroup(16)
test_tp_show(sg .* permutedims(sg), str)

# -------------------------------
# MultTable
# -------------------------------
str = """
      6×6 MultTable{SymOperation{3}}:
      ───────┬──────────────────────────────────────────
             │     1  3₀₀₁⁺  3₀₀₁⁻   2₀₀₁  6₀₀₁⁻  6₀₀₁⁺
      ───────┼──────────────────────────────────────────
           1 │     1  3₀₀₁⁺  3₀₀₁⁻   2₀₀₁  6₀₀₁⁻  6₀₀₁⁺
       3₀₀₁⁺ │ 3₀₀₁⁺  3₀₀₁⁻      1  6₀₀₁⁻  6₀₀₁⁺   2₀₀₁
       3₀₀₁⁻ │ 3₀₀₁⁻      1  3₀₀₁⁺  6₀₀₁⁺   2₀₀₁  6₀₀₁⁻
        2₀₀₁ │  2₀₀₁  6₀₀₁⁻  6₀₀₁⁺      1  3₀₀₁⁺  3₀₀₁⁻
       6₀₀₁⁻ │ 6₀₀₁⁻  6₀₀₁⁺   2₀₀₁  3₀₀₁⁺  3₀₀₁⁻      1
       6₀₀₁⁺ │ 6₀₀₁⁺   2₀₀₁  6₀₀₁⁻  3₀₀₁⁻      1  3₀₀₁⁺
      ───────┴──────────────────────────────────────────
      """
test_tp_show(MultTable(pointgroup("6")), str)

str = """
      4×4 MultTable{SymOperation{3}}:
      ──────────────┬────────────────────────────────────────────────────────
                    │            1  {2₀₁₀|0,0,½}            -1  {m₀₁₀|0,0,½}
      ──────────────┼────────────────────────────────────────────────────────
                  1 │            1  {2₀₁₀|0,0,½}            -1  {m₀₁₀|0,0,½}
       {2₀₁₀|0,0,½} │ {2₀₁₀|0,0,½}             1  {m₀₁₀|0,0,½}            -1
                 -1 │           -1  {m₀₁₀|0,0,½}             1  {2₀₁₀|0,0,½}
       {m₀₁₀|0,0,½} │ {m₀₁₀|0,0,½}            -1  {2₀₁₀|0,0,½}             1
      ──────────────┴────────────────────────────────────────────────────────
      """
test_tp_show(MultTable(spacegroup(13)), str)

# -------------------------------
# KVec
# -------------------------------
for v in (KVec, RVec)
    test_tp_show(v("0,0,.5+u"), "[0, 0, 1/2+α]")
    test_tp_show(v("1/2+α,β+α,1/4"), "[1/2+α, α+β, 1/4]")
    test_tp_show(v("β,-α"), "[β, -α]")
    @test repr(MIME"text/plain"(), v("y,γ,u")) == repr(MIME"text/plain"(), v("β,w,x"))
end

# -------------------------------
# AbstractGroup
# -------------------------------
str = """
      PointGroup{3} ⋕21 (6) with 6 operations:
       1
       3₀₀₁⁺
       3₀₀₁⁻
       2₀₀₁
       6₀₀₁⁻
       6₀₀₁⁺"""
test_tp_show(pointgroup("6", Val(3)), str)

str = """
      SpaceGroup{3} ⋕213 (P4₁32) with 24 operations:
       1
       {2₀₀₁|½,0,½}
       {2₀₁₀|0,½,½}
       {2₁₀₀|½,½,0}
       3₁₁₁⁺
       {3₋₁₁₋₁⁺|½,½,0}
       {3₋₁₁₁⁻|½,0,½}
       {3₋₁₋₁₁⁺|0,½,½}
       3₁₁₁⁻
       {3₋₁₁₁⁺|0,½,½}
       {3₋₁₋₁₁⁻|½,½,0}
       {3₋₁₁₋₁⁻|½,0,½}
       {2₁₁₀|¾,¼,¼}
       {2₋₁₁₀|¾,¾,¾}
       {4₀₀₁⁻|¼,¼,¾}
       {4₀₀₁⁺|¼,¾,¼}
       {4₁₀₀⁻|¾,¼,¼}
       {2₀₁₁|¼,¾,¼}
       {2₀₋₁₁|¾,¾,¾}
       {4₁₀₀⁺|¼,¼,¾}
       {4₀₁₀⁺|¾,¼,¼}
       {2₁₀₁|¼,¼,¾}
       {4₀₁₀⁻|¼,¾,¼}
       {2₋₁₀₁|¾,¾,¾}"""
test_tp_show(spacegroup(213, Val(3)), str)

str = """
      SiteGroup{2} ⋕17 (p6mm) at 2b = [1/3, 2/3] with 6 operations:
       1
       {3⁺|1,1}
       {3⁻|0,1}
       {m₁₁|1,1}
       m₁₀
       {m₀₁|0,1}"""
sg = spacegroup(17,Val(2))
wps = wyckoffs(17, Val(2))
test_tp_show(sitegroup(sg, wps[end-1]), str)

# -------------------------------
# LGIrrep
# -------------------------------
str = """
4-element Collection{LGIrrep{3}} for ⋕16 (P222) at Γ = [0, 0, 0]:
Γ₁┌    1: 1
  ├ 2₁₀₀: 1
  ├ 2₀₁₀: 1
  └ 2₀₀₁: 1

Γ₂┌    1: 1
  ├ 2₁₀₀: -1
  ├ 2₀₁₀: -1
  └ 2₀₀₁: 1

Γ₃┌    1: 1
  ├ 2₁₀₀: 1
  ├ 2₀₁₀: -1
  └ 2₀₀₁: -1

Γ₄┌    1: 1
  ├ 2₁₀₀: -1
  ├ 2₀₁₀: 1
  └ 2₀₀₁: -1"""
test_tp_show(lgirreps(16)["Γ"], str)

str = """
4-element Collection{LGIrrep{3}}:
 #undef
 #undef
 #undef
 #undef"""
test_tp_show(similar(lgirreps(16)["Γ"]), str)

str = """
Γ₅┌             1: ⎡ 1  0 ⎤
  │                ⎣ 0  1 ⎦
  ├  {2₁₀₀|0,½,¼}: ⎡ 1   0 ⎤
  │                ⎣ 0  -1 ⎦
  ├  {2₀₁₀|0,½,¼}: ⎡ -1  0 ⎤
  │                ⎣  0  1 ⎦
  ├          2₀₀₁: ⎡ -1   0 ⎤
  │                ⎣  0  -1 ⎦
  ├        -4₀₀₁⁺: ⎡  0  1 ⎤
  │                ⎣ -1  0 ⎦
  ├        -4₀₀₁⁻: ⎡ 0  -1 ⎤
  │                ⎣ 1   0 ⎦
  ├  {m₁₁₀|0,½,¼}: ⎡  0  -1 ⎤
  │                ⎣ -1   0 ⎦
  ├ {m₋₁₁₀|0,½,¼}: ⎡ 0  1 ⎤
  └                ⎣ 1  0 ⎦"""
test_tp_show(lgirreps(122)["Γ"][end], str)

str = """
Γ₆┌     1: 1
  ├  2₀₀₁: -1
  ├ 3₀₀₁⁺: exp(-0.6667iπ)
  ├ 3₀₀₁⁻: exp(0.6667iπ)
  ├ 6₀₀₁⁺: exp(-0.3333iπ)
  └ 6₀₀₁⁻: exp(0.3333iπ)"""
test_tp_show(lgirreps(168)["Γ"][end], str)

@test repr(lgirreps(168)["Γ"]) == "[Γ₁, Γ₂, Γ₃, Γ₄, Γ₅, Γ₆] (Γ = [0, 0, 0])"

# -------------------------------
# PGIrrep
# -------------------------------
str = """
Γ₄Γ₆┌     1: ⎡ 1  0 ⎤
    │        ⎣ 0  1 ⎦
    ├ 3₀₀₁⁺: ⎡ -0.5+0.866im             0 ⎤
    │        ⎣            0  -0.5-0.866im ⎦
    ├ 3₀₀₁⁻: ⎡ -0.5-0.866im             0 ⎤
    │        ⎣            0  -0.5+0.866im ⎦
    ├  2₀₀₁: ⎡ -1   0 ⎤
    │        ⎣  0  -1 ⎦
    ├ 6₀₀₁⁻: ⎡ 0.5-0.866im            0 ⎤
    │        ⎣           0  0.5+0.866im ⎦
    ├ 6₀₀₁⁺: ⎡ 0.5+0.866im            0 ⎤
    └        ⎣           0  0.5-0.866im ⎦"""
pgirs = pgirreps("6")
pgirs′ = realify(pgirs)
test_tp_show(pgirs′[end], str)
@test summary(pgirs) == "6-element Collection{PGIrrep{3}}"
@test summary(pgirs′) == "4-element Collection{PGIrrep{3}}"

@test repr(realify(pgirreps("4/m"))) == "[Γ₁⁺, Γ₁⁻, Γ₂⁺, Γ₂⁻, Γ₃⁺Γ₄⁺, Γ₃⁻Γ₄⁻]"
@test repr(realify(pgirreps("4/m"; mulliken=true))) == "[Ag, Aᵤ, Bg, Bᵤ, Eg, Eᵤ]"

# -------------------------------
# CharacterTable
# -------------------------------
str = """
CharacterTable{3} for ⋕21 (6):
───────┬────────────────────
       │ Γ₁  Γ₂  Γ₃Γ₅  Γ₄Γ₆
───────┼────────────────────
     1 │  1   1     2     2
 3₀₀₁⁺ │  1   1    -1    -1
 3₀₀₁⁻ │  1   1    -1    -1
  2₀₀₁ │  1  -1     2    -2
 6₀₀₁⁻ │  1  -1    -1     1
 6₀₀₁⁺ │  1  -1    -1     1
───────┴────────────────────
"""
test_tp_show(characters(pgirs′), str)

str = """
CharacterTable{3} for ⋕230 (Ia-3d) at P = [1/2, 1/2, 1/2]:
─────────────────┬──────────────────────
                 │     P₁      P₂    P₃
─────────────────┼──────────────────────
               1 │      2       2     4
    {2₁₀₀|0,0,½} │      0       0     0
    {2₀₁₀|½,0,0} │      0       0     0
    {2₀₀₁|0,½,0} │      0       0     0
           3₁₁₁⁺ │     -1      -1     1
           3₁₁₁⁻ │     -1      -1     1
  {3₋₁₁₁⁺|½,0,0} │   -1im    -1im   1im
  {3₋₁₁₁⁻|0,½,0} │    1im     1im  -1im
 {3₋₁₁₋₁⁻|0,½,0} │   -1im    -1im   1im
 {3₋₁₁₋₁⁺|0,0,½} │    1im     1im  -1im
 {3₋₁₋₁₁⁻|0,0,½} │   -1im    -1im   1im
 {3₋₁₋₁₁⁺|½,0,0} │    1im     1im  -1im
  {-4₁₀₀⁺|¼,¼,¾} │ -1+1im   1-1im     0
  {-4₁₀₀⁻|¾,¼,¼} │  1-1im  -1+1im     0
  {-4₀₁₀⁺|¾,¼,¼} │ -1+1im   1-1im     0
  {-4₀₁₀⁻|¼,¾,¼} │  1-1im  -1+1im     0
  {-4₀₀₁⁺|¼,¾,¼} │ -1+1im   1-1im     0
  {-4₀₀₁⁻|¼,¼,¾} │  1-1im  -1+1im     0
    {m₁₁₀|¾,¼,¼} │      0       0     0
   {m₋₁₁₀|¼,¼,¼} │      0       0     0
    {m₁₀₁|¼,¼,¾} │      0       0     0
   {m₋₁₀₁|¼,¼,¼} │      0       0     0
    {m₀₁₁|¼,¾,¼} │      0       0     0
   {m₀₋₁₁|¼,¼,¼} │      0       0     0
─────────────────┴──────────────────────
"""
test_tp_show(characters(lgirreps(230)["P"]), str)

# -------------------------------
# BandRepSet and BandRep
# -------------------------------
brs = bandreps(42, 3)
str = """
BandRepSet (⋕42): 6 BandReps, sampling 17 LGIrreps (spin-1 w/ TR)
────┬────────────────────────
    │ 4a  4a  4a  4a  8b  8b
    │ A₁  A₂  B₁  B₂  A   B
────┼────────────────────────
 Γ₁ │ 1   ·   ·   ·   1   ·
 Γ₂ │ ·   1   ·   ·   1   ·
 Γ₃ │ ·   ·   ·   1   ·   1
 Γ₄ │ ·   ·   1   ·   ·   1
 T₁ │ 1   ·   ·   ·   ·   1
 T₂ │ ·   1   ·   ·   ·   1
 T₃ │ ·   ·   ·   1   1   ·
 T₄ │ ·   ·   1   ·   1   ·
 Y₁ │ 1   ·   ·   ·   ·   1
 Y₂ │ ·   1   ·   ·   ·   1
 Y₃ │ ·   ·   ·   1   1   ·
 Y₄ │ ·   ·   1   ·   1   ·
 Z₁ │ 1   ·   ·   ·   1   ·
 Z₂ │ ·   1   ·   ·   1   ·
 Z₃ │ ·   ·   ·   1   ·   1
 Z₄ │ ·   ·   1   ·   ·   1
 L₁ │ 1   1   1   1   2   2
────┼────────────────────────
 μ  │ 1   1   1   1   2   2
────┴────────────────────────
  KVecs: Γ, T, Y, Z, L"""
test_tp_show(brs, str)

test_tp_show(brs[1],   "1-band BandRep (A₁↑G at 4a):\n [Γ₁, T₁, Y₁, Z₁, L₁]")
test_tp_show(brs[end], "2-band BandRep (B↑G at 8b):\n [Γ₃+Γ₄, T₁+T₂, Y₁+Y₂, Z₃+Z₄, 2L₁]")

str = "[(8b|A), (8b|B), (4a|A₁), (4a|A₂), (4a|B₂), (4a|B₁)]"
@test repr(calc_bandreps(42, Val(3))) == str

end # @testset