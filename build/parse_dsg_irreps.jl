# Parse the raw HTML cached by `crawl_dsg_irreps.jl` into little group irreps of the double
# space groups. Kept separate from the crawler because parsing is re-run often and the crawl
# only when the cache is incomplete.
#
# ## What a cached page looks like
#
# Each `out/sg<N>-<klab>.html` holds two sections; we want only the first, "Matrices of the
# representations of the *little group*" (the second gives the full space group
# representation, which we do not use). That section is a single top-level table:
#
#   header row : [ "Matrix presentation" | "Seitz Symbol" | <irrep labels…> ]
#   data rows  : [ (W|w) 3×4 | SU(2) 2×2 | Seitz symbol | <irrep matrices…> ]
#
# Note the SU(2) column carries no header, so a data row has one cell *more* than the header.
#
# ## The double group is given explicitly
#
# The data rows list the **whole** double group: first the |G| spatial operations, then the
# same operations again prefixed "d" in their Seitz symbol, meaning Ē·g. A double-valued
# irrep has `D(Ē·g) = -D(g)` and a single-valued one `D(Ē·g) = +D(g)`, which is what lets
# Crystalline store only the first |G| matrices and recover the rest from the sign.
# `check_barred_coset` checks this on a parsed page; `test/double/lgirreps.jl` checks it on
# the data that is ultimately written.
#
# ## Matrix entries
#
# Entries are `0`, `±1`, `±i`, small integers, or `e` with a superscript exponent such as
# `iπ(t₁+t₂)`, `i2πt₁`, `-iπ/4`, `iπ(1/2+w)`. Two kinds of symbol appear in those exponents:
#   - `t₁,t₂,t₃`: a *general* lattice translation. Bilbao writes every matrix as a function of
#     it, i.e. `D({R|w+t}) = exp(2πi k⋅t) D({R|w})`; the listed operation is `t = 0`.
#   - `u,v,w`: the free parameters of a non-special **k**, i.e. Crystalline's `αβγ`.
# Since `D({R|w})` already contains `exp(2πi k⋅w)`, Crystalline's stored `P` (which excludes
# the translation phase, applied instead at evaluation time from the `translations` field) is
#     P = D_bilbao(t=0) * exp(-2πi k⋅τ),   τ = w = translation(op).
# A useful consequence, asserted below: after dividing out that phase, `P` must be free of
# `u,v,w` — it is the `αβγ`-independent part.
#
# ## Bundled k-manifolds
#
# `crawl_dsg_irreps.jl` resolves Bilbao's bundled form entries into one page per manifold, so
# every `out/` page here is exactly one k-manifold. A page lacking a "Seitz Symbol" header is
# a stray selection page and is rejected rather than silently skipped.

using Gumbo
using Crystalline
using StaticArrays
using LinearAlgebra: dot
using Pkg.Artifacts: ensure_artifact_installed

const NOVALS = Dict{Symbol,Real}()  # for entries that carry no `t` or `αβγ` symbols

"""
    crawl_datadir() --> String

The directory holding the captured Bilbao pages, downloading them first if necessary: the
capture is not kept in the repository, but published as a release asset and declared as the
`dsg_crawl` artifact (see `Artifacts.toml` and `build/DATA-RELEASE.md`). It also contains
`crawl_dsg_irreps.jl`, the crawler that produced it.

`@artifact_str` is not usable here — it searches upwards from this file for an
`Artifacts.toml` and stops at `build/Project.toml` — so the declaration is named explicitly.
"""
crawl_datadir() = ensure_artifact_installed("dsg_crawl",
                                            joinpath(dirname(@__DIR__), "Artifacts.toml"))

# ---------------------------------------------------------------------------------------- #
# Gumbo helpers

children_elements(el::HTMLElement) = [c for c in el.children if c isa HTMLElement]
tag_children(el::HTMLElement, t::Symbol) = filter(c -> tag(c) == t, children_elements(el))

# HTML5 parsing inserts a <tbody> between <table> and its <tr>s, so a table's rows are not in
# general its direct children; Gumbo reproduces that faithfully.
function rows_of(tbl::HTMLElement)
    direct = tag_children(tbl, :tr)
    isempty(direct) || return direct
    return [r for b in tag_children(tbl, :tbody) for r in tag_children(b, :tr)]
end

function text_of(x)
    io = IOBuffer()
    _text_of(io, x)
    return strip(replace(String(take!(io)), r"\s+" => " "))
end
_text_of(io, t::HTMLText) = print(io, t.text)
_text_of(io, el::HTMLElement) = foreach(c -> _text_of(io, c), el.children)

"""
    has_overline(el) --> Bool

Whether `el` contains a `<font style="text-decoration:overline">`. Bilbao uses the overline
both for double-valued irrep labels (`X̄₃`) *and* for roto-inversions in Seitz symbols (`1̄`),
and plain text extraction loses it — `{1̄|0,0,0}` reads as `{1|0,0,0}`, i.e. inversion
masquerading as the identity. So the markup, never the text, is what must be inspected.
"""
function has_overline(el::HTMLElement)
    if tag(el) === :font
        sty = get(attrs(el), "style", "")
        occursin("overline", sty) && return true
    end
    return any(c -> c isa HTMLElement && has_overline(c), el.children)
end
has_overline(::HTMLText) = false

# ---------------------------------------------------------------------------------------- #
# Locating the little group table

"""
    littlegroup_table(doc) --> HTMLElement

The single top-level table of the "little group" section. Throws if the page has no such
table — which is how a stray selection page (see the module header) is caught.
"""
function littlegroup_table(doc::HTMLDocument)
    best = nothing
    walk(el) = begin
        if el isa HTMLElement
            if tag(el) === :table
                rows = rows_of(el)
                if !isempty(rows) && occursin("Seitz", text_of(first(rows))) && best === nothing
                    best = el
                end
            end
            foreach(walk, el.children)
        end
    end
    walk(doc.root)
    best === nothing && error("no little group table found: not an irrep page (a stray \
                               selection page?)")
    return best
end

# ---------------------------------------------------------------------------------------- #
# Row/column extraction

struct DsgPage
    sgnum     :: Int
    klabel    :: String
    kvec      :: String
    seitz     :: Vector{String}        # length 2N: the spatial ops, then the Ē-barred ones
    opmatrix  :: Vector{String}        # the (W|w) cell, verbatim
    su2       :: Vector{String}        # the SU(2) cell, verbatim
    irlabels  :: Vector{String}        # with 'ˢ' marking a double-valued (overlined) label
    isdouble  :: Vector{Bool}
    entries   :: Matrix{Vector{String}} # [op, irrep] -> flat list of entry strings
end

"""
    datarow_cells(rows, ncol, nir, off) --> (seitz, opmatrix, su2, entries)

The cells of the data rows of an irrep table, verbatim. Each such row holds the operation's
matrix cell, its SU(2) cell, its Seitz symbol, and then one cell per irrep; `off` is the
index of the matrix cell (1 on a space group page, 2 on a point group page, which prefixes
a numbering column). Rows with a cell count other than `ncol` are not data rows.
"""
function datarow_cells(rows, ncol::Integer, nir::Integer, off::Integer)
    datarows = [r for r in rows[2:end] if length(tag_children(r, :td)) == ncol]
    N = length(datarows)
    seitz    = Vector{String}(undef, N)
    opmatrix = Vector{String}(undef, N)
    su2      = Vector{String}(undef, N)
    entries  = Matrix{Vector{String}}(undef, N, nir)
    for (i, r) in enumerate(datarows)
        cs = tag_children(r, :td)
        opmatrix[i] = text_of(cs[off])
        su2[i]      = text_of(cs[off+1])
        # mark roto-inversions explicitly, since `text_of` drops the overline
        c = cs[off+2]
        seitz[i] = has_overline(c) ? "‾" * text_of(c) : text_of(c)
        for j in 1:nir
            entries[i, j] = entry_strings(cs[off+2+j])
        end
    end
    return seitz, opmatrix, su2, entries
end

function parse_page(path::AbstractString)
    src = read(path, String)
    cut = findfirst("Matrices of the representations of the group", src)
    src = cut === nothing ? src : src[1:first(cut)-1]
    doc = parsehtml(src)
    tbl = littlegroup_table(doc)
    rows = rows_of(tbl)

    hdr = tag_children(first(rows), :td)
    # header = [matrix presentation, Seitz symbol, irrep labels…]; data rows insert an
    # unheaded SU(2) column, hence the offset of 3 rather than 2 below.
    irlabel_cells = hdr[3:end]
    isdouble = has_overline.(irlabel_cells)
    irlabels = text_of.(irlabel_cells)

    seitz, opmatrix, su2, entries = datarow_cells(rows, length(hdr)+1, length(irlabels), 1)

    sgnum, klabel = let m = match(r"sg(\d+)-(.+)\.html$", basename(path))
        m === nothing && error("cannot read sgnum/klabel from filename $(basename(path))")
        parse(Int, m.captures[1]), String(m.captures[2])
    end
    kvec = let m = match(r"wave vector k\s*1?\s*=\s*\(([^)]*)\)", text_of(doc.root))
        m === nothing ? "" : String(m.captures[1])
    end
    return DsgPage(sgnum, klabel, kvec, seitz, opmatrix, su2, irlabels, isdouble, entries)
end

"""
    parse_realities(path) --> Dict{String, Int8}

The realities stated on the page, keyed by the Bilbao irrep label with `ˢ` marking a
double-valued irrep. They are given in the page's second section, the representations of
the full space group (the star of **k**), whose table header lists each irrep as e.g.
`*X3 (0)`: `(1)` real, `(-1)` pseudoreal, `(0)` complex. This is the reality that the
Herring criterion determines.
"""
function parse_realities(path::AbstractString)
    src = read(path, String)
    cut = findfirst("Matrices of the representations of the group", src)
    cut === nothing && error("no full space group representations on $(basename(path))")
    rows = rows_of(littlegroup_table(parsehtml(src[first(cut):end])))
    realities = Dict{String, Int8}()
    for c in tag_children(first(rows), :td)[3:end]
        m = match(r"^\s*\*?\s*(\S+)\s*\((-?\d)\)\s*$", text_of(c))
        m === nothing && error("unexpected irrep label cell $(repr(text_of(c))) on \
                                $(basename(path))")
        realities[m[1] * (has_overline(c) ? string(SPINFUL_MARK) : "")] = parse(Int8, m[2])
    end
    return realities
end

"""
    entry_strings(cell) --> Vector{String}

The matrix entries of one irrep cell, in row-major order, each kept verbatim (e.g. `"0"`,
`"-i"`, `"e^(iπ(t1+t2))"`). Each entry is an innermost `<td>` of the cell's nested table;
superscripts are rendered as `^(…)` so that the exponent survives as structure rather than
being run together with the base by text extraction.
"""
function entry_strings(cell::HTMLElement)
    tbls = innermost_tables(cell)
    # a 1×1 "matrix" is written as bare text, with no nested table at all
    isempty(tbls) && return [render_entry(cell)]
    return [render_entry(c) for t in tbls for r in rows_of(t) for c in tag_children(r, :td)]
end

# the tables below `el` that contain no further table, in document order
function innermost_tables(el, out = HTMLElement[])
    el isa HTMLElement || return out
    n = length(out)
    foreach(c -> innermost_tables(c, out), el.children)
    tag(el) === :table && length(out) == n && push!(out, el)
    return out
end

function render_entry(el)
    io = IOBuffer()
    _render(io, el)
    return replace(strip(String(take!(io))), r"\s+" => "")
end
_render(io, t::HTMLText) = print(io, t.text)
function _render(io, el::HTMLElement)
    if tag(el) === :sup
        print(io, "^(", text_of(el), ")")
    elseif tag(el) === :sub
        print(io, text_of(el))            # subscripts are indices: t<sub>1</sub> -> t1
    else
        foreach(c -> _render(io, c), el.children)
    end
end

# ---------------------------------------------------------------------------------------- #

pages(; dir = joinpath(crawl_datadir(), "out")) = sort(readdir(dir; join = true))

# ---------------------------------------------------------------------------------------- #
# Evaluating matrix entries
#
# An entry is `0`, a small integer, `±i`, a surd, a parenthesised complex scalar, or
# `e^(<exponent>)`. The exponents are **bilinear**, not linear — e.g. `i2πt3w`,
# `iπ(t1+t2+2t3w)` — which is exactly right: the phase is 2πi𝐤⋅𝐭 and 𝐤 carries the free
# parameters, so components of 𝐭 multiply αβγ.
#
# Symbols: `i` (imaginary unit), `π`, `t1,t2,t3` (a general lattice translation) and `u,v,w`
# (the free parameters of 𝐤, i.e. Crystalline's αβγ). Juxtaposition means multiplication.

const EXP_SYMBOL_CHARS = ('i', 'π', 'u', 'v', 'w')

"""
    expr_of(s) --> Expr

Turn a Bilbao scalar — a matrix entry *or* an exponent — into a Julia expression, making the
implicit multiplications explicit: `"iπ(t1+t2+2t3w)"` → `im*π*(t1+t2+2*t3*w)`,
`"e^(i3π/4)√2/2"` → `exp(im*3*π/4)*sqrt(2)/2`, `"(1-i)/2"` → `(1-im)/2`.

Deliberately one grammar for both, rather than an enumeration of the shapes seen so far: the
set of shapes grew repeatedly as the crawl advanced, breaking each such parser in turn.
Anything outside the grammar throws, so a genuinely new construct is a loud failure rather
than a silent misreading.
"""
function expr_of(s::AbstractString)
    toks = String[]
    i = firstindex(s)
    while i ≤ lastindex(s)
        c = s[i]
        if isdigit(c)
            j = i
            while j < lastindex(s) && isdigit(s[nextind(s, j)]); j = nextind(s, j); end
            push!(toks, s[i:j]); i = nextind(s, j)
        elseif c == 't'                       # t1, t2, t3
            j = nextind(s, i)
            (j ≤ lastindex(s) && isdigit(s[j])) ||
                error("expected a digit after 't' in exponent $(repr(s))")
            push!(toks, s[i:j]); i = nextind(s, j)
        elseif c == 'e' && i < lastindex(s) && s[nextind(s, i)] == '^'
            push!(toks, "exp"); i = nextind(s, i, 2)   # 'e^(' → exp(  (the '(' follows)
        elseif c == '√'
            m = match(r"^√(\d+)", SubString(s, i))
            m === nothing && error("expected digits after '√' in $(repr(s))")
            push!(toks, "sqrt(" * m.captures[1] * ")"); i += ncodeunits(m.match)
        elseif c ∈ EXP_SYMBOL_CHARS || c ∈ ('+', '-', '*', '/', '(', ')')
            push!(toks, string(c)); i = nextind(s, i)
        elseif isspace(c)
            i = nextind(s, i)
        else
            error("unexpected character $(repr(c)) in exponent $(repr(s))")
        end
    end
    io = IOBuffer()
    prev = ""
    for t in toks
        !isempty(prev) && _ends_value(prev) && _starts_value(t) && print(io, '*')
        print(io, t == "i" ? "im" : t)
        prev = t
    end
    return Meta.parse(String(take!(io)))
end
_ends_value(t)   = t == ")" || endswith(t, ")") || isdigit(first(t)) ||
                   (t ∉ ("exp",) && first(t) ∈ ('t', EXP_SYMBOL_CHARS...))
_starts_value(t) = t == "(" || t == "exp" || startswith(t, "sqrt") || isdigit(first(t)) ||
                   first(t) ∈ ('t', EXP_SYMBOL_CHARS...)

"""
    eval_expr(ex, vals) --> ComplexF64

Evaluate a parsed exponent against `vals`, a mapping of `:t1,:t2,:t3,:u,:v,:w` to numbers.
Deliberately a small interpreter rather than `eval`: the input is scraped from a web page, and
an interpreter that only knows `+ - * /` cannot be made to do anything else.
"""
function eval_expr(ex, vals::AbstractDict{Symbol,<:Real})
    ex isa Number && return ComplexF64(ex)
    if ex isa Symbol
        ex === :im && return im
        ex === :π  && return ComplexF64(π)
        haskey(vals, ex) && return ComplexF64(vals[ex])
        error("unknown symbol $(ex) in exponent")
    end
    (ex isa Expr && ex.head === :call) || error("unexpected expression $(ex) in exponent")
    op = ex.args[1]
    as = [eval_expr(a, vals) for a in ex.args[2:end]]
    op === :+ && return sum(as)
    op === :* && return prod(as)
    op === :/ && return as[1] / as[2]
    op === :- && return length(as) == 1 ? -as[1] : as[1] - as[2]
    op === :exp  && return exp(only(as))
    op === :sqrt && return sqrt(only(as))
    error("disallowed operator $(op) in expression")
end

"""
    entry_value(e, vals) --> ComplexF64

Numeric value of a single matrix entry. An entry is a product of optional factors — a leading
sign, `i`, an exponential `e^(…)`, a surd `√n`, an integer, and a divisor `/n` — e.g. `-1`,
`i`, `e^(iπ(1/2+w))`, `√2/2`, `i√2/2`, `e^(i7π/12)√2/2`. Unrecognised forms **throw**:
silently ignoring an unknown factor would corrupt a matrix.
"""
function entry_value(e::AbstractString, vals::AbstractDict{Symbol,<:Real})
    e == "0" && return zero(ComplexF64)
    return eval_expr(expr_of(e), vals)
end

# byte index of the ')' matching the '(' at byte index `open_at`
function _matching_paren(s::AbstractString, open_at::Int)
    depth = 0
    i = open_at
    while i ≤ lastindex(s)
        c = s[i]
        c == '(' && (depth += 1)
        c == ')' && (depth -= 1; depth == 0 && return i)
        i = nextind(s, i)
    end
    error("unbalanced parentheses in entry $(repr(s))")
end

"""
    sparse_entry(e) --> Union{Nothing, Tuple{Int,Int,String}}

Bilbao gives irreps of dimension > 4 **sparsely**, listing only the non-zero elements as
`(i;j):x` (the page says so in its own preamble). Returns `(i, j, x)` for such an entry, or
`nothing` for an ordinary dense one.
"""
function sparse_entry(e::AbstractString)
    m = match(r"^\((\d+);(\d+)\):(.*)$", e)
    m === nothing && return nothing
    return (parse(Int, m.captures[1]), parse(Int, m.captures[2]), String(m.captures[3]))
end

"""
    matrix_at(p, iop, iir; t = (0,0,0), αβγ = (0.123, 0.187, 0.243)) --> Matrix{ComplexF64}

The irrep matrix of operation `iop` in irrep `iir`, evaluated at lattice translation `t` and
free parameters `αβγ`. The default `t = 0` is what picks out the operation *as listed*; the
default `αβγ` is a generic point, deliberately not a special value that could make
inequivalent matrices coincide by accident.
"""
function matrix_at(p::DsgPage, iop::Integer, iir::Integer;
                   t = (0, 0, 0), αβγ = (0.123, 0.187, 0.243))
    vals = Dict{Symbol,Real}(:t1 => t[1], :t2 => t[2], :t3 => t[3],
                             :u => αβγ[1], :v => αβγ[2], :w => αβγ[3])
    es = p.entries[iop, iir]
    sparse = sparse_entry.(es)
    if any(!isnothing, sparse)
        all(!isnothing, sparse) ||
            error("irrep $(p.irlabels[iir]) at op $iop mixes sparse and dense entries")
        # the dimension cannot be read off the entry count here; take it from the largest
        # index, which is exact because every row and column of a unitary irrep matrix
        # carries at least one non-zero element
        d = maximum(x -> max(x[1], x[2]), sparse)
        M = zeros(ComplexF64, d, d)
        for (r, c, x) in sparse
            M[r, c] = entry_value(x, vals)
        end
        return M
    end
    d = isqrt(length(es))
    d*d == length(es) || error("irrep $(p.irlabels[iir]) at op $iop has $(length(es)) entries, \
                                not a square matrix")
    return ComplexF64[entry_value(es[(r-1)*d + c], vals) for r in 1:d, c in 1:d] # row-major
end

"""
    check_barred_coset(p; kws...) --> NamedTuple

Verify that the second half of the operation list is the Ē-barred coset of the first, by
evaluating both and checking `D(Ē·g) = ±D(g)`, with the sign set by whether the irrep is
double-valued. The comparison must be numeric: Bilbao routinely absorbs the sign into the
exponent as a shift by π, e.g. for sg 100 at `B̄₃`, where `D(g) = e^(-iπ(1/2-v))` and
`D(Ē·g) = e^(iπ(1/2+v))` differ by -1 although no string manipulation would show it.
"""
function check_barred_coset(p::DsgPage; atol = 1e-10, kws...)
    N = length(p.seitz)
    isodd(N) && return (ok = false, reason = "odd number of operations")
    n = N ÷ 2
    for i in 1:n, j in eachindex(p.irlabels)
        s = p.isdouble[j] ? -1 : +1
        A = matrix_at(p, i,   j; kws...)
        B = matrix_at(p, i+n, j; kws...)
        size(A) == size(B) || return (ok = false, reason = "shape mismatch at op $i, irrep $j")
        isapprox(B, s*A; atol) ||
            return (ok = false, reason = "D(Ē·g) ≠ $(s > 0 ? "+" : "-")D(g) at op $i, \
                                          irrep $(p.irlabels[j])")
    end
    return (ok = true, reason = "")
end

# ---------------------------------------------------------------------------------------- #
# Operations, and alignment with Crystalline's own little groups

"""
    operation_of(p, i) --> SymOperation{3}

The `i`th operation of the page, from its `(W|w)` cell: twelve whitespace-separated tokens in
row-major order. The identity row carries the symbolic general translation `t1,t2,t3`, which
is read as zero — that row *is* the identity, written with the general translation attached so
that the page can display the Bloch phase.
"""
function operation_of(p::DsgPage, i::Integer)
    toks = split(p.opmatrix[i])
    length(toks) == 12 ||
        error("expected 12 entries in the (W|w) cell of op $i, got $(length(toks))")
    v = map(toks) do s
        s ∈ ("t1", "t2", "t3") && return 0.0
        m = match(r"^(-?\d+)/(\d+)$", s)
        m === nothing ? parse(Float64, s) : parse(Float64, m[1]) / parse(Float64, m[2])
    end
    W = SMatrix{3,3,Float64}(v[1], v[5], v[9], v[2], v[6], v[10], v[3], v[7], v[11])
    w = SVector{3,Float64}(v[4], v[8], v[12])
    return SymOperation{3}(W, w)
end

"""
    su2_of(p, i) --> SU2

The SU(2) element of row `i`, from its 2×2 cell (`a` and `b` are its first two entries).
"""
function su2_of(p::DsgPage, i::Integer)
    a, b = entry_value.(split(strip(p.su2[i]))[1:2], Ref(NOVALS))
    return SU2(a, b)
end

"""
    KLABEL_BILBAO2CDML

Bilbao spells five CDML k-labels in ASCII. Verified 2026-09-17 that, with just these five, the
k-label sets of Bilbao and `lgirreps` agree for all 230 space groups in both directions.
"""
const KLABEL_BILBAO2CDML = Dict("GM" => "Γ", "LD" => "Λ", "DT" => "Δ", "SM" => "Σ", "GP" => "Ω")
cdml_klabel(klab::AbstractString) = get(KLABEL_BILBAO2CDML, klab, klab)

"""
    align_operations(p, lg) --> NamedTuple

Match the page's spatial operations (the first half; the second is the Ē-coset) onto the
operations of `lg`, returning the permutation `perm` with `operations(lg)[perm[i]]` ≡ the
page's `i`th operation.

Comparison is modulo the lattice translations of the **primitive** cell, i.e. it must be told
the `centering`: in a centred lattice the centring vector is a lattice translation, so two
operations differing by e.g. `(0,½,½)` in an F-centred group are the same operation.
Comparing only modulo *integer* translations fails on a handful of the centred groups.

Returns `ok = false` with a reason rather than throwing, so that a sweep can tabulate how the
settings differ instead of stopping at the first mismatch.
"""
function align_operations(p::DsgPage, lg, cntr::Char = centering(p.sgnum, 3))
    n = length(p.seitz) ÷ 2
    length(lg) == n || return (ok = false, perm = Int[],
                               reason = "order mismatch: page has $n, littlegroups has $(length(lg))")
    perm = Vector{Int}(undef, n)
    for i in 1:n
        opᵇ = operation_of(p, i)
        j = findfirst(opᶜ -> isapprox(opᵇ, opᶜ, cntr, true), operations(lg))
        j === nothing && return (ok = false, perm = Int[],
                                 reason = "no match for page op $i ($(p.seitz[i]))")
        perm[i] = j
    end
    length(unique(perm)) == n || return (ok = false, perm, reason = "match is not a bijection")
    return (ok = true, perm, reason = "")
end

"""
    crystalline_matrix(p, iop, iir, kv; αβγ) --> Matrix{ComplexF64}

The irrep matrix in **Crystalline's** convention, i.e. `LGIrrep.matrices`, which excludes the
translation phase (that is reapplied at evaluation time from the `translations` field).

Bilbao writes every matrix as a function of a general lattice translation 𝐭,
`D({R|w+t}) = exp(2πi𝐤⋅𝐭)·D({R|w})`, and the listed operation is `𝐭 = 0`; but `D({R|w})` still
contains `exp(2πi𝐤⋅w)`. So

    P = D_bilbao(t=0) · exp(-2πi 𝐤⋅τ),   τ = w = translation(op).

Since `P` is by construction the αβγ-independent part, evaluating it at two different αβγ must
give the same matrix — which is what `check_convention` asserts, and which is the empirical
test of this whole convention mapping.
"""
function crystalline_matrix(p::DsgPage, iop::Integer, iir::Integer, kv; αβγ)
    D = matrix_at(p, iop, iir; t = (0, 0, 0), αβγ)
    τ = translation(operation_of(p, iop))
    return D .* cispi(-2 * dot(kv(collect(αβγ)), τ))
end

"""
    check_convention(p, lg; αβγ₁, αβγ₂, atol) --> NamedTuple

Verify that `crystalline_matrix` is genuinely αβγ-independent, by evaluating it at two generic,
unrelated values of αβγ and comparing. A failure means the assumed relation between Bilbao's
matrices and Crystalline's `(matrices, translations)` split is wrong for that page.
"""
function check_convention(p::DsgPage, lg;
                          αβγ₁ = (0.123, 0.187, 0.243), αβγ₂ = (0.317, 0.059, 0.431),
                          atol = 1e-9)
    kv = position(lg)
    n = length(p.seitz) ÷ 2
    for i in 1:n, j in eachindex(p.irlabels)
        A = crystalline_matrix(p, i, j, kv; αβγ = αβγ₁)
        B = crystalline_matrix(p, i, j, kv; αβγ = αβγ₂)
        isapprox(A, B; atol) || return (ok = false,
            reason = "αβγ-dependence survives at op $i ($(p.seitz[i])), irrep $(p.irlabels[j])")
    end
    return (ok = true, reason = "")
end

# ---------------------------------------------------------------------------------------- #
# Label normalisation

"""
    SPINFUL_MARK

The modifier letter `ˢ` (`U+02E2`), appended to mark a double-valued irrep: `Γ₇⁺ˢ`, `WA₅ˢ`.

CDML and Bilbao instead overline the k-label symbol (`<font overline>WA</font><sub>5</sub>`),
but an overbar is awkward here: it is a combining mark that must follow the *complete* letter
run rather than each letter, it collides with the overbar of a roto-inversion, and it renders
unreliably. The `ˢ` spelling already exists in Crystalline for the spinful EBR labels of Bilbao's
tabulated band representations, which append it in the same position, so that `klabel` needs
no special casing.
"""
const SPINFUL_MARK = 'ˢ'

"""
    cdml_irlabel(lab, isdouble) --> String

Normalise a Bilbao irrep label to Crystalline's CDML spelling: `"GM1+"` → `"Γ₁⁺"`,
`"X3"` (overlined) → `"X₃ˢ"`, `"WA5"` (overlined) → `"WA₅ˢ"`.

Every label in the crawl matches `^[A-Z]+[0-9]+[+-]?\$` (checked over all 20088 of them), so the
parse is exhaustive rather than best-effort; anything else throws.
"""
function cdml_irlabel(lab::AbstractString, isdouble::Bool)
    m = match(r"^([A-Z]+)(\d+)([+-]?)$", lab)
    m === nothing && error("unexpected irrep label $(repr(lab))")
    klab_b, num, sgn = m.captures
    return string(cdml_klabel(klab_b),
                  Crystalline.subscriptify(num),
                  isempty(sgn) ? "" : Crystalline.supscriptify(sgn),
                  isdouble ? SPINFUL_MARK : "")
end

cdml_irlabels(p::DsgPage) = [cdml_irlabel(l, d) for (l, d) in zip(p.irlabels, p.isdouble)]
