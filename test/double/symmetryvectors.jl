using Crystalline, Test

datafile = joinpath(pkgdir(Crystalline), "data", "irreps", "lgs", "3d",
                    "irreps_data_spinful.jld2")
if !isfile(datafile)
    @warn "spinful irrep data not found; skipping spinful symmetry vector tests" datafile
else

@testset "Spinful symmetry vectors" begin
    lgirsd = lgirreps(221, Val(3); spinful=Val(true))
    lgirsv = [lgirsd[klab] for klab in ("Γ", "X", "M", "R")]

    # parsing: a 2-band symmetry vector over the double-valued irreps
    n = parse(SymmetryVector, "[Γ₆ˢ, X₆ˢ, M₆ˢ, R₆ˢ]", lgirsv)
    @test n isa SymmetryVector{3, DLGIrrep{3}}
    @test occupation(n) == 2
    @test irreps(n) === lgirsv
    @test sum(multiplicities(n)[1]) == 1

    # arithmetic stays within the same type
    @test 2n + n isa SymmetryVector{3, DLGIrrep{3}}
    @test 2n + n == 3n && occupation(3n) == 6
    @test iszero(n - n) && zero(n) == n - n

    # construction from a raw vector and labels, and printing (irrep type not shown)
    nv = Vector(n)
    @test SymmetryVector(nv, irreplabels(n), lgirsd) == n
    @test only(SymmetryVectors([nv], irreplabels(n), lgirsd)) == n
    @test startswith(sprint(show, MIME"text/plain"(), n),
                     "20-irrep SymmetryVector{3} (spinful):")

    # spinless and spinful symmetry vectors are distinct types
    br = first(calc_bandreps(221))
    @test !(SymmetryVector(br) isa SymmetryVector{3, DLGIrrep{3}})
    @test isspinful(n) && !isspinful(br) && !isspinful(SymmetryVector(br))
    @test isspinful(first(lgirsv[1])) && !isspinful(first(lgirreps(221)["Γ"]))
    @test startswith(sprint(show, MIME"text/plain"(), br),
                     "40-irrep NewBandRep{3} (spinless):")
end

@testset "Spinful irreps in symmetry eigenvalue analysis and primitivization" begin
    lgirsd = lgirreps(229, Val(3); spinful=Val(true)) # body-centered
    αβγ = [0.1, 0.2, 0.3]
    for lgirs in values(lgirsd)
        lgirs′ = primitivize(lgirs)
        @test lgirs′ isa Collection{DLGIrrep{3}}
        @test klabel(lgirs′) == klabel(lgirs)
        for (lgir, lgir′) in zip(lgirs, lgirs′)
            @test all(splat(≈), zip(lgir(αβγ), lgir′(αβγ)))
        end
    end
    @test primitivize(lgirsd) isa Dict{String, Collection{DLGIrrep{3}}}

    # the characters of each irrep, as "symmetry eigenvalues", are identified as that irrep
    lgirs = lgirsd["Γ"]
    annotations = collect_irrep_annotations([characters(lgir) for lgir in lgirs], lgirs)
    @test last.(annotations) == label.(lgirs)
end

end # if isfile(datafile)
