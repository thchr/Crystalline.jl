# Convert the crawled Bilbao double point group pages into Crystalline's JLD2 layout, i.e.
# the spinful counterpart of `data/irreps/pgs/3d/irreps_data.jld2`, written to
# `data/irreps/pgs/3d/irreps_data_spinful.jld2` in the same layout: per point group,
# `matrices`, `realities` and `cdmls`. Only the double-valued irreps are stored.
#
# ## What a cached page looks like
#
# Each `pg/<iuc>.html` (with '/' spelled '_') holds a character table, and then the section
# "Matrices of the representations of the group", a single table:
#
#   header row : [ "N" | "Matrix presentation" | "Seitz Symbol" | <irrep labels…> ]
#   data rows  : [ N | W 3×3 | SU(2) 2×2 | Seitz symbol | <irrep matrices…> ]
#
# As on the space group pages, the SU(2) column has no header, and the rows list the whole
# double group: the |G| operations, then the Ē-barred ones. Unlike the space group pages,
# each irrep label is followed by its reality, `(1)`, `(-1)` or `(0)`.
#
# ## Checks
#
# Each operation of `pointgroup(iuc)` is matched to the page row with the same rotation part
# *and* the same SU(2) element as `su2` assigns it, so no assumption is made about the order
# of rows. Bilbao's single-valued irreps must then reproduce `pgirreps(iuc)` exactly, since
# that data was crawled from the same Bilbao program.
#
#   julia --project=build build/write_dsg_pgirreps.jl [outpath]

include(joinpath(@__DIR__, "parse_dsg_irreps.jl"))
using JLD2

const DEFAULT_OUT = joinpath(dirname(@__DIR__), "data", "irreps", "pgs", "3d",
                             "irreps_data_spinful.jld2")
const NOVALS = Dict{Symbol,Real}()

pg_path(iuc) = joinpath(CRAWL_DIR, "pg", replace(iuc, '/' => "_") * ".html")

"""
    parse_pg_page(path) --> (page::DsgPage, realities::Vector{Int8})

Parse a double point group page into a `DsgPage` (with `sgnum = 0` and `klabel = "GM"`), so
that the space group page tooling applies, together with the stated irrep realities.
"""
function parse_pg_page(path::AbstractString)
    doc = parsehtml(read(path, String))
    rows = rows_of(littlegroup_table(doc))

    hdr = tag_children(first(rows), :td)
    irlabel_cells = hdr[4:end]
    isdouble = has_overline.(irlabel_cells)
    irlabels, realities = String[], Int8[]
    for c in irlabel_cells
        m = match(r"^\s*(\S+)\s*\((-?\d)\)\s*$", text_of(c))
        m === nothing && error("unexpected irrep label cell $(repr(text_of(c))) in $path")
        push!(irlabels, m[1])
        push!(realities, parse(Int8, m[2]))
    end

    datarows = [r for r in rows[2:end] if length(tag_children(r, :td)) == length(hdr) + 1]
    N, nir = length(datarows), length(irlabels)
    opmatrix = Vector{String}(undef, N)
    su2      = Vector{String}(undef, N)
    seitz    = Vector{String}(undef, N)
    entries  = Matrix{Vector{String}}(undef, N, nir)
    for (i, r) in enumerate(datarows)
        cs = tag_children(r, :td)
        opmatrix[i] = text_of(cs[2])
        su2[i]      = text_of(cs[3])
        seitz[i]    = has_overline(cs[4]) ? "‾" * text_of(cs[4]) : text_of(cs[4])
        for j in 1:nir
            entries[i, j] = entry_strings(cs[4+j])
        end
    end
    page = DsgPage(0, "GM", "0,0,0", seitz, opmatrix, su2, irlabels, isdouble, entries)
    return page, realities
end

# the rotation part of row `i`, from its 3×3 cell (row-major)
function rotation_of(p::DsgPage, i::Integer)
    v = parse.(Int, split(p.opmatrix[i]))
    length(v) == 9 || error("expected 9 entries in the W cell of row $i, got $(length(v))")
    return SMatrix{3,3,Int}(v[1], v[4], v[7], v[2], v[5], v[8], v[3], v[6], v[9])
end
function su2_of(p::DsgPage, i::Integer)
    return SU2(entry_value.(split(strip(p.su2[i])), Ref(NOVALS))[1:2]...)
end

"""
    collect_pg(iuc) --> (matrices, realities, cdmls)

The double-valued irreps of point group `iuc`, with matrices in the operation order of
`pointgroup(iuc)`.
"""
function collect_pg(iuc::String)
    p, realities = parse_pg_page(pg_path(iuc))
    check_barred_coset_numeric(p).ok || error("$iuc: D(Ē·g) ≠ ±D(g)")
    pg = pointgroup(iuc, Val(3))
    hexagonal = Crystalline._ishexagonal(pg)

    # the page row of each of our operations
    rows = map(operations(pg)) do op
        u = su2(op, hexagonal)
        W = round.(Int, rotation(op))
        i = findfirst(i -> rotation_of(p, i) == W && isapprox(su2_of(p, i), u),
                      eachindex(p.seitz))
        i === nothing && error("$iuc: no row matches $(seitz(op)) with SU(2) element $u")
        i
    end
    length(p.seitz) == 2length(pg) || error("$iuc: expected $(2length(pg)) rows")

    # Bilbao's single-valued irreps must be those of `pgirreps`
    pgirs = pgirreps(iuc, Val(3))
    for j in findall(!, p.isdouble)
        lab = cdml_irlabel(p.irlabels[j], false)
        k = findfirst(pgir -> label(pgir) == lab, pgirs)
        k === nothing && error("$iuc: no irrep $lab in `pgirreps`")
        all(i -> matrix_at(p, rows[i], j) ≈ pgirs[k].matrices[i], eachindex(rows)) ||
            error("$iuc: irrep $lab disagrees with `pgirreps`")
        realities[j] == Integer(reality(pgirs[k])) ||
            error("$iuc: reality of irrep $lab disagrees with `pgirreps`")
    end
    count(!, p.isdouble) == length(pgirs) || error("$iuc: single-valued irreps missing")

    js = findall(p.isdouble)
    matrices = [[matrix_at(p, i, j) for i in rows] for j in js]
    cdmls = [cdml_irlabel(p.irlabels[j], true) for j in js]
    return matrices, realities[js], cdmls
end

function write_dsg_pgirreps(outpath::AbstractString = DEFAULT_OUT)
    mkpath(dirname(outpath))
    # Crystalline keeps the spinful data file open for reading, which blocks overwriting it
    isassigned(Crystalline.DPGIRREPS_JLDFILE) && close(Crystalline.DPGIRREPS_JLDFILE[])
    JLD2.jldopen(outpath, "w") do f
        for iuc in Crystalline.PG_IUCs[3]
            matrices, realities, cdmls = collect_pg(iuc)
            key = Crystalline._unmangle_pgiuclab(iuc)
            f["$key/matrices"]  = matrices
            f["$key/realities"] = realities
            f["$key/cdmls"]     = cdmls
        end
    end
    println("wrote $(outpath)  ($(round(filesize(outpath)/2^10; digits=1)) KiB)")
    return outpath
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_dsg_pgirreps(length(ARGS) ≥ 1 ? ARGS[1] : DEFAULT_OUT)
end
