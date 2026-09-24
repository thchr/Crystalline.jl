# Check Bilbao's single-valued little group irreps against the ISOTROPY irreps that
# `lgirreps` loads. The two must agree up to a change of basis, i.e. in their characters.
#
# Bilbao's single-valued irreps come along for free on the double space group pages, and
# `write_dsg_irreps.jl` stores them in `irreps_data_spinless_bilbao.jld2` for exactly this
# purpose: they are the part of the crawl whose answer is already known, so their agreement
# with ISOTROPY is evidence that the crawl, the parse, and the convention mapping of
# `parse_dsg_irreps.jl` are right — and hence that the double-valued irreps can be trusted.
#
# Crystalline never loads that file, so it is not shipped with the package: it is a lazy
# artifact (see `Artifacts.toml`), downloaded once per depot on the first access below.

using Crystalline, Test
using JLD2
using LinearAlgebra: norm
using Pkg.Artifacts: ensure_artifact_installed

"""
    validate_bilbao_vs_isotropy(path; sgnums = 1:230) --> Vector{String}

Describe every disagreement between the Bilbao irreps in `path` and the ISOTROPY irreps of
`lgirreps`; an empty return means the two datasets agree.
"""
function validate_bilbao_vs_isotropy(path::AbstractString;
                                     sgnums = 1:MAX_SGNUM[3], atol = 1e-10)
    problems = String[]
    JLD2.jldopen(path, "r") do jldfile
        for sgnum in sgnums
            lgirsd  = lgirreps(sgnum, Val(3))
            lgirsd′ = lgirreps(sgnum, Val(3), Crystalline.LGS_JLDFILES[3][], jldfile)
            if Set(keys(lgirsd′)) ≠ Set(keys(lgirsd))
                push!(problems, "sg $sgnum: k-label sets differ")
                continue
            end
            for (klab, lgirs) in lgirsd
                lgirs′ = lgirsd′[klab]
                if sort(label.(lgirs)) ≠ sort(label.(lgirs′))
                    push!(problems, "sg $sgnum, $klab: irrep labels differ")
                    continue
                end
                for lgir in lgirs
                    lgir′ = lgirs′[findfirst(ir -> label(ir) == label(lgir), lgirs′)]
                    # the two datasets distribute the Bloch phase differently between
                    # matrices and translations, so compare characters instead: they agree
                    # for all αβγ iff they agree at αβγ = 0 and, wherever they are non-zero,
                    # the αβγ-dependent parts of their Bloch phases agree
                    χ, χ′ = characters(lgir), characters(lgir′)
                    kabc = parts(position(lgir))[2]
                    Δτs = lgir.translations .- lgir′.translations
                    what = "sg $sgnum, $(label(lgir))"
                    χ ≈ χ′ || push!(problems, "$what: characters differ")
                    reality(lgir) == reality(lgir′) ||
                        push!(problems, "$what: realities differ")
                    all(((c, Δτ),) -> abs(c) < atol || norm(kabc'*Δτ) < atol,
                        zip(χ, Δτs)) ||
                        push!(problems, "$what: αβγ-dependent phases differ")
                end
            end
        end
    end
    return problems
end

# `@artifact_str` cannot be used here: it searches upwards from this file for an
# `Artifacts.toml` and stops at the first `Project.toml` it meets, which is `test/`'s own
datadir = try
    ensure_artifact_installed("bilbao_spinless_irreps",
                              joinpath(pkgdir(Crystalline), "Artifacts.toml"))
catch err
    @warn "could not obtain the `bilbao_spinless_irreps` artifact; skipping the comparison \
           of Bilbao's single-valued irreps against ISOTROPY's" err
    nothing
end

if !isnothing(datadir)

@testset "Bilbao's single-valued irreps vs. ISOTROPY's" begin
    problems = validate_bilbao_vs_isotropy(
                    joinpath(datadir, "irreps_data_spinless_bilbao.jld2"))
    @test isempty(problems)
    foreach(problem -> println("    ", problem), problems)
end

end # if !isnothing(datadir)
