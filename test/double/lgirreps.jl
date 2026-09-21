using Crystalline, Test, LinearAlgebra
using Crystalline: check_multtable_vs_ir, matrices, can_intersect, TEST_αβγs

datafile = joinpath(pkgdir(Crystalline), "data", "irreps", "lgs", "3d",
                    "irreps_data_spinful.jld2")
if !isfile(datafile)
    @warn "spinful irrep data not found; skipping tests of `DLGIrrep`s" datafile
else

@testset "Double-valued little group irreps" begin

@test @inferred(lgirreps(1, Val(3), Val(true))) isa Dict{String, Collection{DLGIrrep{3}}}
@test lgirreps(1, Val(3), Val(false)) == lgirreps(1, Val(3))
@test lgirreps(1, 3, true) == lgirreps(1, Val(3), Val(true))
@test lgirreps(1, 3) == lgirreps(1, 3, false)
@test @inferred(littlegroups(1, Val(3), Val(true))) isa Dict{String, DLittleGroup{3}}
@test littlegroups(1, 3, true) == littlegroups(1, Val(3), Val(true))
@test_throws DomainError littlegroups(1, Val(2), Val(true))
@test_throws DomainError lgirreps(1, Val(2), Val(true))

αβγ = TEST_αβγs[3]
for sgnum in 1:MAX_SGNUM[3]
    dlgirsd = lgirreps(sgnum, Val(3), Val(true))
    lgirsd  = lgirreps(sgnum, Val(3))
    @test Set(keys(dlgirsd)) == Set(keys(lgirsd))
    for (klab, dlgirs) in dlgirsd
        dlg = group(first(dlgirs))
        n = order(dlg) ÷ 2
        @test n == order(first(lgirsd[klab]))
        @test all(dlgir -> endswith(label(dlgir), 'ˢ'), dlgirs)
        @test all(dlgir -> klabel(dlgir) == klab, dlgirs)

        # translations are stored only if they give an αβγ-dependent phase; otherwise, the
        # phase is part of the matrices
        kabc = parts(position(dlg))[2]
        @test all(dlgirs) do dlgir
            τs = dlgir.translations
            all(iszero, τs) || any(τ -> norm(kabc'*τ) > 1e-10, τs)
        end

        # the barred half of the double group is represented by `D(Ēg) = -D(g)`
        @test all(dlgirs) do dlgir
            Ds = matrices(dlgir)
            all(i -> Ds[i+n] == -Ds[i], 1:n)
        end

        # the double-valued irreps are a complete set of the double group's irreps with
        # `D(Ē) = -𝟙`: their dimensions square-sum to `|G|`, not `2|G|`
        @test sum(dlgir -> irdim(dlgir)^2, dlgirs) == n

        # character orthogonality, summed over the double group
        χs = [characters(dlgir, αβγ) for dlgir in dlgirs]
        @test all(χ -> dot(χ, χ) ≈ 2n, χs)
        @test all(((a, b),) -> a == b || abs(dot(χs[a], χs[b])) < 1e-10,
                  Iterators.product(eachindex(χs), eachindex(χs)))

        # class characters at Γ (elsewhere, the irreps can be ray representations, whose
        # characters are not class functions); needs the classes of the double group
        if klab == "Γ"
            ct = classcharacters(dlgirs)
            ws, X = length.(classes(ct)), matrix(ct)
            @test sum(ws) == 2n
            @test X' * (ws .* X) ≈ 2n * I
            @test [subduction_count(a, b) for a in dlgirs, b in dlgirs] == I
        end

        # the irreps respect the double group's multiplication (computed in a primitive
        # basis, and including ray-representation phases at nonsymmorphic k-points)
        for dlgir in dlgirs, αβγ′ in (nothing, αβγ)
            @test all(check_multtable_vs_ir(dlgir, αβγ′))
        end

        # the ray-representation phases depend only on the spatial operations
        israyᵈ, αᵈ = israyrep(first(dlgirs))
        isray,  α  = israyrep(first(lgirsd[klab]))
        @test israyᵈ == isray
        @test αᵈ ≈ repeat(α, 2, 2)
    end
end

# Compatibility relations from a point onto a line, as tabulated by Bilbao's DCOMPREL; for
# each irrep at the point, its decomposition at the line. Single-valued irreps are included,
# via the spinless irreps
@testset "Compatibility relations (Bilbao's DCOMPREL)" begin
    function compatibility(sgnum, klabᴳ, klabᴴ)
        lgirsdˢ, lgirsdᵈ = lgirreps(sgnum, Val(3)), lgirreps(sgnum, Val(3), Val(true))
        d = Dict{String, Dict{String, Int}}()
        for lgirsd in (lgirsdˢ, lgirsdᵈ)
            lgirsᴳ, lgirsᴴ = lgirsd[klabᴳ], lgirsd[klabᴴ]
            # the αβγ that places the line at the point
            intersection = can_intersect(position(first(lgirsᴴ)), position(first(lgirsᴳ)))
            @test intersection.bool
            αβγ = intersection.αβγ
            for Dᴳ in lgirsᴳ
                ns = [label(Dᴴ) => subduction_count(Dᴳ, Dᴴ, αβγ) for Dᴴ in lgirsᴴ]
                d[label(Dᴳ)] = Dict(filter(p -> p[2] > 0, ns))
            end
        end
        return d
    end

    # Ia-3d (230), P → Λ
    @test compatibility(230, "P", "Λ") == Dict(
        "P₁"  => Dict("Λ₃" => 1),
        "P₂"  => Dict("Λ₃" => 1),
        "P₃"  => Dict("Λ₁" => 1, "Λ₂" => 1, "Λ₃" => 1),
        "P₄ˢ" => Dict("Λ₅ˢ" => 1),
        "P₅ˢ" => Dict("Λ₄ˢ" => 1),
        "P₆ˢ" => Dict("Λ₆ˢ" => 1),
        "P₇ˢ" => Dict("Λ₅ˢ" => 1, "Λ₆ˢ" => 1),
        "P₈ˢ" => Dict("Λ₄ˢ" => 1, "Λ₆ˢ" => 1))

    # Ia-3d (230), Γ → Σ
    @test compatibility(230, "Γ", "Σ") == Dict(
        "Γ₁⁺" => Dict("Σ₁" => 1),
        "Γ₁⁻" => Dict("Σ₄" => 1),
        "Γ₂⁺" => Dict("Σ₂" => 1),
        "Γ₂⁻" => Dict("Σ₃" => 1),
        "Γ₃⁺" => Dict("Σ₁" => 1, "Σ₂" => 1),
        "Γ₃⁻" => Dict("Σ₃" => 1, "Σ₄" => 1),
        "Γ₄⁺" => Dict("Σ₂" => 1, "Σ₃" => 1, "Σ₄" => 1),
        "Γ₄⁻" => Dict("Σ₁" => 1, "Σ₂" => 1, "Σ₃" => 1),
        "Γ₅⁺" => Dict("Σ₁" => 1, "Σ₃" => 1, "Σ₄" => 1),
        "Γ₅⁻" => Dict("Σ₁" => 1, "Σ₂" => 1, "Σ₄" => 1),
        "Γ₆ˢ" => Dict("Σ₅ˢ" => 1),
        "Γ₇ˢ" => Dict("Σ₅ˢ" => 1),
        "Γ₈ˢ" => Dict("Σ₅ˢ" => 1),
        "Γ₉ˢ" => Dict("Σ₅ˢ" => 1),
        "Γ₁₀ˢ" => Dict("Σ₅ˢ" => 2),
        "Γ₁₁ˢ" => Dict("Σ₅ˢ" => 2))

    # P-6m2 (187), H → P
    @test compatibility(187, "H", "P") == Dict(
        "H₁"  => Dict("P₁" => 1),  "H₂"  => Dict("P₁" => 1),
        "H₃"  => Dict("P₂" => 1),  "H₄"  => Dict("P₂" => 1),
        "H₅"  => Dict("P₃" => 1),  "H₆"  => Dict("P₃" => 1),
        "H₇ˢ" => Dict("P₄ˢ" => 1), "H₈ˢ" => Dict("P₄ˢ" => 1),
        "H₉ˢ" => Dict("P₅ˢ" => 1), "H₁₀ˢ" => Dict("P₅ˢ" => 1),
        "H₁₁ˢ" => Dict("P₆ˢ" => 1), "H₁₂ˢ" => Dict("P₆ˢ" => 1))
end

end # @testset "Double-valued little group irreps"
end # if isfile(datafile)

# Bilbao's single-valued irreps, written alongside the double-valued ones, must agree with
# the ISOTROPY irreps loaded by `lgirreps`, up to a change of basis (i.e., in characters)
datafile_bilbao = joinpath(dirname(datafile), "irreps_data_spinless_bilbao.jld2")
if isfile(datafile_bilbao)
@testset "Single-valued irreps: Bilbao vs. ISOTROPY" begin
    jldfile = Crystalline.JLD2.jldopen(datafile_bilbao, "r")
    try
        for sgnum in 1:MAX_SGNUM[3]
            lgirsd = lgirreps(sgnum, Val(3))
            lgirsd_bilbao = lgirreps(sgnum, Val(3), Val(false),
                                     Crystalline.LGS_JLDFILES[3][], jldfile)
            @test Set(keys(lgirsd_bilbao)) == Set(keys(lgirsd))
            for (klab, lgirs) in lgirsd
                lgirs_bilbao = lgirsd_bilbao[klab]
                @test sort(label.(lgirs)) == sort(label.(lgirs_bilbao))
                for lgir in lgirs
                    i = findfirst(ir -> label(ir) == label(lgir), lgirs_bilbao)
                    lgir_bilbao = lgirs_bilbao[i]
                    # the two datasets distribute the Bloch phase differently between
                    # matrices and translations, so compare characters instead: they agree
                    # for all αβγ iff they agree at αβγ = 0 and, wherever they are non-zero,
                    # the αβγ-dependent parts of their Bloch phases agree
                    χ, χ′ = characters(lgir), characters(lgir_bilbao)
                    kabc = parts(position(lgir))[2]
                    Δτs = lgir.translations .- lgir_bilbao.translations
                    @test χ ≈ χ′
                    @test all(zip(χ, Δτs)) do (c, Δτ)
                        abs(c) < 1e-10 || norm(kabc'*Δτ) < 1e-10
                    end
                end
            end
        end
    finally
        close(jldfile)
    end
end
end # if isfile(datafile_bilbao)
