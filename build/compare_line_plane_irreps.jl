# Compare the single-valued little group irreps of the plane (2D) and line (1D) groups, as
# derived from their parent space groups (see `LinePlaneIrreps.jl`), with those of `lgirreps`.
# The irreps themselves must agree, up to their labels; the labels that differ are printed.
#
#   julia --project=. build/compare_line_plane_irreps.jl

include(joinpath(@__DIR__, "line_plane_irreps.jl"))
using .LinePlaneIrreps, Crystalline
using Crystalline: TEST_αβγs

for D in (2, 1), sgnum in 1:MAX_SGNUM[D]
    lgirsd = lgirreps(sgnum, Val(D))
    lgirsd′ = derive_lgirreps_from_spacegroup(sgnum, Val(D); spinful=Val(false))
    αβγ = TEST_αβγs[D]
    for (klab, lgirs′) in lgirsd′
        lgirs = lgirsd[klab]
        length(lgirs) == length(lgirs′) || error("⋕$sgnum ($(D)D), $klab: missing irreps")
        mismatches = String[]
        for lgir′ in lgirs′
            χ′ = characters(lgir′, αβγ)
            i = findfirst(lgir -> isapprox(characters(lgir, αβγ), χ′; atol=1e-10), lgirs)
            isnothing(i) && error("⋕$sgnum ($(D)D): derived irrep $(label(lgir′)) has no \
                                   counterpart in `lgirreps`")
            label(lgirs[i]) == label(lgir′) ||
                push!(mismatches, "$(label(lgirs[i])) ⇒ $(label(lgir′))")
        end
        isempty(mismatches) ||
            println("⋕$sgnum ($(D)D): `lgirreps` ⇒ derived labels: ", join(mismatches, ", "))
    end
end
