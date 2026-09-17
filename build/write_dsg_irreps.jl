# Convert the parsed Bilbao double space group pages into Crystalline's own JLD2 layout, i.e.
# the spinful counterpart of `data/irreps/lgs/3d/irreps_data.jld2`.
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
# and `translations` holds **our** operation's `τ`. Those differ by a lattice vector `Δ` for 181
# operations (all centred lattices), and the phases cancel exactly — see `parse_dsg_irreps.jl`.
# `P` is evaluated at `αβγ = 0`: it is αβγ-independent by construction (verified on every page),
# and zero keeps the arithmetic cleanest.
#
# ## Realities are NOT determined here
#
# Bilbao's little group tables do not state the reality type, and computing it for a spinful
# irrep needs the Herring criterion evaluated in the double group (the Z₂ cocycle σ). So every
# reality is written as `UNDEF`, to be filled in at Stage 5. Do not mistake this for a claim
# that the irreps are of undefined reality.
#
# ## Usage
#
#   julia --project=build build/write_dsg_irreps.jl [outpath] [--sgnums 1:230]

include(joinpath(@__DIR__, "parse_dsg_irreps.jl"))
using JLD2

const DEFAULT_OUT = joinpath(CRAWL_DIR, "irreps_data_spinful.jld2") # gitignored while in flux

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
        alltrivial = all(iszero, τs)

        push!(klabs, klab)
        push!(matrices_list, Ps)
        push!(translations_list, [alltrivial ? nothing : τs for _ in eachindex(p.irlabels)])
        push!(realities_list, fill(Int8(2), length(p.irlabels)))   # UNDEF; see header
        push!(cdml_list, cdml_irlabels(p))
    end
    return klabs, matrices_list, translations_list, realities_list, cdml_list
end

function write_dsg_irreps(outpath::AbstractString = DEFAULT_OUT; sgnums = 1:230)
    mkpath(dirname(outpath))
    JLD2.jldopen(outpath, "w") do f
        for sgnum in sgnums
            lgs = littlegroups(sgnum, Val(3))
            klabs, ms, τs, rs, cs = collect_sg(sgnum, lgs)
            f["$(sgnum)/klab_list"]         = klabs
            f["$(sgnum)/matrices_list"]     = ms
            f["$(sgnum)/translations_list"] = τs
            f["$(sgnum)/realities_list"]    = rs
            f["$(sgnum)/cdml_list"]         = cs
            sgnum % 20 == 0 && println("  … sg $sgnum")
        end
    end
    println("wrote $(outpath)  ($(round(filesize(outpath)/2^20; digits=1)) MiB)")
    return outpath
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = length(ARGS) ≥ 1 && !startswith(ARGS[1], "--") ? ARGS[1] : DEFAULT_OUT
    i = findfirst(==("--sgnums"), ARGS)
    sgnums = i === nothing ? (1:230) : eval(Meta.parse(ARGS[i+1]))
    write_dsg_irreps(out; sgnums)
end
