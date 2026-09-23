using Crystalline, Test, LinearAlgebra
using Crystalline: check_multtable_vs_ir, matrices, PG_IUCs

datafile = joinpath(pkgdir(Crystalline), "data", "irreps", "pgs", "3d",
                    "irreps_data_spinful.jld2")
if !isfile(datafile)
    @warn "spinful point group irrep data not found; skipping tests of `DPGIrrep`s" datafile
else

@testset "Double-valued point group irreps" begin

@test @inferred(pointgroup("4mm", Val(3); spinful=Val(true))) isa DPointGroup{3}
@test pointgroup("4mm", 3; spinful=true) == pointgroup(13, Val(3), 1; spinful=Val(true))
@test pointgroup("4mm", Val(3); spinful=Val(false)) == pointgroup("4mm", Val(3))
@test_throws DomainError pointgroup("4mm", Val(2); spinful=Val(true))

@test @inferred(pgirreps("4mm", Val(3); spinful=Val(true))) isa Collection{DPGIrrep{3}}
@test pgirreps("4mm", 3; spinful=true) == pgirreps(13, Val(3); spinful=Val(true))
@test pgirreps("4mm", Val(3); spinful=Val(false)) == pgirreps("4mm", Val(3))
@test_throws DomainError pgirreps("4mm", Val(2); spinful=Val(true))

for iuc in PG_IUCs[3]
    pgirs = pgirreps(iuc, Val(3); spinful=Val(true))
    pg = group(first(pgirs))
    n = order(pg) ÷ 2
    @test n == order(pointgroup(iuc, Val(3)))
    @test all(pgir -> endswith(label(pgir), 'ˢ'), pgirs)

    # the barred half of the double group is represented by `D(Ēg) = -D(g)`
    @test all(pgir -> all(i -> matrices(pgir)[i+n] == -matrices(pgir)[i], 1:n), pgirs)

    # the double-valued irreps are a complete set of the double group's irreps with
    # `D(Ē) = -𝟙`, and are orthogonal over the double group
    @test sum(pgir -> irdim(pgir)^2, pgirs) == n
    χs = characters.(pgirs)
    @test all(((a, b),) -> isapprox(dot(χs[a], χs[b]), a == b ? 2n : 0; atol=1e-10),
              Iterators.product(eachindex(χs), eachindex(χs)))

    # character tables, per operation and per conjugacy class (class-weighted orthogonality)
    @test characters(pgirs) isa CharacterTable{DSymOperation{3}}
    ct = classcharacters(pgirs)
    ws, X = length.(classes(ct)), matrix(ct)
    @test sum(ws) == 2n
    @test X' * (ws .* X) ≈ 2n * I

    # the irreps respect the double group's multiplication
    mt = MultTable(pg)
    @test all(pgir -> all(check_multtable_vs_ir(mt, pgir)), pgirs)

    # Bilbao's stated realities agree with the Frobenius-Schur criterion in the double group
    @test all(pgir -> calc_reality(pgir) == reality(pgir), pgirs)

    # time reversal: with T² = -1, Kramers degeneracy makes every co-representation
    # even-dimensional
    @test all(iseven ∘ irdim, realify(pgirs))
end

# At Γ, the double-valued little group irreps must be those of the double point group of the
# space group: compare characters, matching each operation by its rotation part and SU(2)
# element
@testset "Little group irreps at Γ vs. point group irreps" begin
datafile_lgs = joinpath(pkgdir(Crystalline), "data", "irreps", "lgs", "3d",
                        "irreps_data_spinful.jld2")
if isfile(datafile_lgs)
for sgnum in 1:MAX_SGNUM[3]
    lgirs = lgirreps(sgnum, Val(3); spinful=Val(true))["Γ"]
    lg = group(first(lgirs))
    iuc = label(Crystalline.find_parent_pointgroup(littlegroups(sgnum)["Γ"]))
    pgirs = pgirreps(iuc, Val(3); spinful=Val(true))
    pg = group(first(pgirs))
    idxs = map(operations(lg)) do op
        findfirst(pgop -> rotation(pgop) ≈ rotation(op) && SU2(pgop) ≈ SU2(op), pg)
    end
    @test all(!isnothing, idxs)
    @test sort(label.(lgirs)) == sort(label.(pgirs))
    for lgir in lgirs
        pgir = pgirs[findfirst(pgir -> label(pgir) == label(lgir), pgirs)]
        @test characters(lgir) ≈ characters(pgir)[idxs]
    end
end
end
end

@testset "Mulliken labels" begin
    for iuc in Crystalline.PG_IUCs[3]
        pgirs  = pgirreps(iuc, Val(3); spinful=Val(true))
        pgirsₘ = pgirreps(iuc, Val(3); spinful=Val(true), mulliken=true)
        @test label.(pgirsₘ) == mulliken.(pgirs)
        @test allunique(label.(pgirsₘ)) && all(l -> endswith(l, 'ˢ'), label.(pgirsₘ))
        # co-reps: `mulliken` of a co-rep agrees with `realify` of Mulliken-labelled irreps
        @test label.(realify(pgirsₘ)) == mulliken.(realify(pgirs))
    end
end

end # @testset "Double-valued point group irreps"
end # if isfile(datafile)
