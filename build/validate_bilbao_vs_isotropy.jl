# Check Bilbao's single-valued little group irreps against the ISOTROPY irreps that
# `lgirreps` loads. The two must agree up to a change of basis, i.e. in their characters.
#
# Bilbao's single-valued irreps come along for free on the double space group pages, and
# `write_dsg_irreps.jl` stores them in `irreps_data_spinless_bilbao.jld2` for exactly this
# purpose: they are the part of the crawl whose answer is already known, so their agreement
# with ISOTROPY is evidence that the crawl, the parse, and the convention mapping of
# `parse_dsg_irreps.jl` are right — and hence that the double-valued irreps can be trusted.
#
# That file is a build artifact and is not shipped with the package, so this check lives
# here rather than in the test suite.
#
#   julia --project=build build/validate_bilbao_vs_isotropy.jl [path]

using Crystalline
using JLD2
using LinearAlgebra: norm

const DEFAULT_PATH = joinpath(dirname(@__DIR__), "data", "irreps", "lgs", "3d",
                              "irreps_data_spinless_bilbao.jld2")

"""
    validate_bilbao_vs_isotropy(path = DEFAULT_PATH; sgnums = 1:230) --> Vector{String}

Describe every disagreement between the Bilbao irreps in `path` and the ISOTROPY irreps of
`lgirreps`; an empty return means the two datasets agree.
"""
function validate_bilbao_vs_isotropy(path::AbstractString = DEFAULT_PATH;
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

if abspath(PROGRAM_FILE) == @__FILE__
    problems = validate_bilbao_vs_isotropy(length(ARGS) ≥ 1 ? ARGS[1] : DEFAULT_PATH)
    if isempty(problems)
        println("Bilbao's single-valued irreps agree with ISOTROPY's, for all space groups")
    else
        println("$(length(problems)) disagreements:")
        foreach(problem -> println("  ", problem), problems)
        exit(1)
    end
end
