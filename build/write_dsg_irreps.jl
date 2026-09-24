# Convert the parsed Bilbao double space group pages into Crystalline's own JLD2 layout, i.e.
# the spinful counterpart of `data/irreps/lgs/3d/irreps_data.jld2`.
#
# Two files are written, both in `data/irreps/lgs/3d/`:
# - `irreps_data_spinful.jld2`: the double-valued irreps, loaded by
#   `lgirreps(…, Val(true))`;
# - `irreps_data_spinless_bilbao.jld2`: Bilbao's single-valued irreps, which are not loaded
#   by the package, but kept for comparison with the ISOTROPY irreps in `irreps_data.jld2`
#   (see `test/bilbao_vs_isotropy.jl`). This file is not committed: it is published as a
#   `data-v*` release asset and consumed as a lazy artifact (see `Artifacts.toml`).
#
# Mirrors `write_littlegroup_irreps.jl`: per space group, the file stores `matrices_list`,
# `translations_list`, `realities_list` and `cdml_list`, each indexed first by k-manifold and
# then by irrep. The **little groups themselves are not rewritten** — the operations are aligned
# onto those already in `littlegroups_data.jld2` (see `align_operations`), so the existing file
# is shared between the spinless and spinful irreps.
#
# ## What is stored
#
# For each irrep, `matrices` holds `P` in *Crystalline's* convention — the translation phase
# excluded, to be reapplied at evaluation time from `translations`:
#
#     P = D_bilbao(t=0) · exp(-2πi 𝐤⋅τ_bilbao)
#
# and `translations` holds **our** operation's `τ`. The two translations can differ by a
# lattice vector in centred lattices, where the phases then cancel exactly. `P` is evaluated
# at `αβγ = 0`: it is αβγ-independent by construction (verified on every page), and zero
# keeps the arithmetic cleanest.
#
# If the phase `exp(2πi 𝐤⋅τ)` does not depend on `αβγ` for any operation (i.e., if
# `kabcᵀτ = 0`; always so at special **k**-points), it is instead folded into the matrices,
# and `translations` is stored as `nothing`. So `translations` is non-`nothing` only if the
# irrep has a genuinely αβγ-dependent phase — as in ISOTROPY's data.
#
# ## Realities
#
# The realities are Bilbao's, as stated for the full space group representations on each
# page (see `parse_realities`); the Herring criterion (`calc_reality`) reproduces them.
#
# ## Usage
#
#   julia --project=build build/write_dsg_irreps.jl [outdir] [--sgnums 1:230]

include(joinpath(@__DIR__, "parse_dsg_irreps.jl"))
using JLD2
using LinearAlgebra: dot, norm

const DEFAULT_OUTDIR = joinpath(dirname(@__DIR__), "data", "irreps", "lgs", "3d")

# the crawler owns the cache layout; restate just the one path we need rather than including it
page_path(sgnum, klab) = joinpath(CRAWL_DIR, "out", "sg$(sgnum)-$(klab).html")

"""
    collect_sg(sgnum, lgs) --> (klabs, matrices_list, translations_list, realities_list, cdml_list)

Assemble one space group's spinful *and* single-valued irreps from its cached pages, with the
irrep matrices permuted into the operation order of `littlegroups(sgnum)`.
"""
function collect_sg(sgnum::Integer, lgs::AbstractDict)
    klabs = String[]
    matrices_list     = Vector{Vector{Vector{Matrix{ComplexF64}}}}()
    translations_list = Vector{Vector{Union{Nothing, Vector{Vector{Float64}}}}}()
    realities_list    = Vector{Vector{Int8}}()
    cdml_list         = Vector{Vector{String}}()

    for (klab, lg) in lgs
        # Bilbao spells five k-labels in ASCII; find the page for this manifold
        bilbao_klab = something(findfirst(==(klab), KLABEL_BILBAO2CDML), klab)
        path = page_path(sgnum, bilbao_klab)
        isfile(path) || error("no cached page for sg $sgnum, k-label $(klab) ($(bilbao_klab))")
        p = parse_page(path)

        r = align_operations(p, lg)
        r.ok || error("sg $sgnum, $(klab): cannot align operations — $(r.reason)")
        invp = invperm(r.perm)          # invp[j] = page index of our j-th operation
        n = length(lg)
        kv = position(lg)

        Ps = map(eachindex(p.irlabels)) do j
            [crystalline_matrix(p, invp[i], j, kv; αβγ = (0.0, 0.0, 0.0)) for i in 1:n]
        end
        τs = [Vector{Float64}(translation(op)) for op in lg]
        k₀, kabc = parts(kv)
        foldphase = all(τ -> norm(kabc' * τ) < 1e-10, τs) # phases independent of αβγ
        if foldphase
            for P in Ps, (i, τ) in enumerate(τs)
                P[i] *= cispi(2 * dot(k₀, τ))
            end
        end

        push!(klabs, klab)
        push!(matrices_list, Ps)
        push!(translations_list, [foldphase ? nothing : τs for _ in eachindex(p.irlabels)])
        realities = parse_realities(path)
        push!(realities_list, [realities[l * (d ? string(SPINFUL_MARK) : "")]
                               for (l, d) in zip(p.irlabels, p.isdouble)])
        push!(cdml_list, cdml_irlabels(p))
    end
    return klabs, matrices_list, translations_list, realities_list, cdml_list
end

function write_dsg_irreps(outdir::AbstractString = DEFAULT_OUTDIR; sgnums = 1:230)
    mkpath(outdir)
    path_double = joinpath(outdir, "irreps_data_spinful.jld2")
    path_single = joinpath(outdir, "irreps_data_spinless_bilbao.jld2")
    # Crystalline keeps the spinful data file open for reading, which blocks overwriting it
    isassigned(Crystalline.DLGIRREPS_JLDFILE) && close(Crystalline.DLGIRREPS_JLDFILE[])
    f_double = JLD2.jldopen(path_double, "w")
    f_single = JLD2.jldopen(path_single, "w")
    try
        for sgnum in sgnums
            lgs = littlegroups(sgnum, Val(3))
            klabs, ms, τs, rs, cs = collect_sg(sgnum, lgs)
            for (f, isdouble) in ((f_double, true), (f_single, false))
                # per k-manifold, the irreps of the requested kind
                idxs = [findall(l -> endswith(l, SPINFUL_MARK) == isdouble, c) for c in cs]
                f["$(sgnum)/klab_list"]         = klabs
                f["$(sgnum)/matrices_list"]     = [m[i] for (m, i) in zip(ms, idxs)]
                f["$(sgnum)/translations_list"] = [τ[i] for (τ, i) in zip(τs, idxs)]
                f["$(sgnum)/realities_list"]    = [r[i] for (r, i) in zip(rs, idxs)]
                f["$(sgnum)/cdml_list"]         = [c[i] for (c, i) in zip(cs, idxs)]
            end
            sgnum % 20 == 0 && println("  … sg $sgnum")
        end
    finally
        close(f_double)
        close(f_single)
    end
    for path in (path_double, path_single)
        println("wrote $(path)  ($(round(filesize(path)/2^20; digits=1)) MiB)")
    end
    return path_double, path_single
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) ≥ 1 && !startswith(ARGS[1], "--") ? ARGS[1] : DEFAULT_OUTDIR
    i = findfirst(==("--sgnums"), ARGS)
    sgnums = i === nothing ? (1:230) : eval(Meta.parse(ARGS[i+1]))
    write_dsg_irreps(outdir; sgnums)
end
