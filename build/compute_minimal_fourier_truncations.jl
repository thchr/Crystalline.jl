# Compute the minimal Fourier truncations `idxmax` for `levelsetlattice` such that a
# *generically modulated* Fourier lattice realizes exactly the targeted space group —
# no accidental extra symmetry (supergroup) and no accidental extra periodicity.
#
# Background: `levelsetlattice(sgnum, Dᵛ, idxmax)` seeds orbits of reciprocal-lattice
# vectors G from the box |Gᵢ| ≤ idxmaxᵢ (orbits are *completed* beyond the box), and
# deletes orbits whose symmetry constraints force their coefficient to zero (systematic
# extinctions). For strongly nonsymmorphic/centered groups, so few symmetry-breaking
# orbits survive at small truncation that a generic modulation is invariant under a strict
# supergroup of the target group.
#
# Criterion (exact; no random sampling): with orbits Oᵢ and unit-modulus coefficient
# vectors vᵢ from `levelsetlattice`, the generic field is
#     f(r) = Re Σᵢ cᵢ Σₘ vᵢₘ exp(2πi Gᵢₘ·r),   cᵢ ∈ ℂ generic and independent.
# `idxmax` is *sufficient* iff the stabilizer of f — within the isometries of a lattice
# with generic parameters of the group's Bravais class — equals the target group exactly:
#
#  (T) no extra translations: {1|w} stabilizes f iff exp(2πi G·w) = 1 for every retained
#      G, so the ℤ-span L of the retained G's must equal
#      L_exp = {m ∈ ℤᴰ : m·t ∈ ℤ ∀ centering translations t}
#      (the reciprocal lattice of the *centered* direct lattice, in conventional coords).
#  (R) no extra point-group operations: for every W in the arithmetic holohedry of the
#      Bravais type (conventional basis) with W ∉ P(G), decide whether {W|w} stabilizes f
#      for SOME w — crucially, w is allowed to depend on the modulation draw, since the
#      physically relevant object is the stabilizer of one realization. Because f is
#      real (a₋G = conj(a_G)), the retained orbits organize into pairs {O, -O}; writing
#      a_G = z_p v_{i,G} on orbit O_i of pair p (z_p a generic complex combination of the
#      two independent modulations of O and -O), the structural requirement is that
#      σ = W⁻ᵀ maps every self-negating orbit onto itself and every pair {O,-O} onto
#      itself — either preserving or SWAPPING the two partner orbits. (Any other image
#      orbit forces |z_p| = |z_q| for independent generic moduli: impossible.) The phase
#      conditions are then, with n(m) the intra-orbit index of the relevant partner:
#        - orbit preserved:  (σGᵢₘ)·w ≡ θᵢₘ (mod 1),  exp(2πiθᵢₘ) = vᵢₘ/vᵢₙ₍ₘ₎
#          [matches the `conds` convention in `levelsetlattice`: vₙ = exp(-2πi(σGₘ)·w)vₘ]
#        - pair swapped (σOᵢ = -Oᵢ ≠ Oᵢ): (σGᵢₘ)·w - s_p ≡ θᵢₘ (mod 1), with
#          exp(2πiθᵢₘ) = vᵢₘ·vᵢₙ₍ₘ₎, n(m) = pos(-σGᵢₘ) in Oᵢ, and s_p = arg(z_p)/π a
#          GENERIC real (one per swapped pair, fixed by the draw).
#      {W|w} is a generic symmetry iff the system is solvable for generic (s_p): after
#      exact integer row reduction over columns (w | s | θ), there must be no reduced row
#      whose pivot lies in the s-block (such a row pins generic phases — measure zero)
#      and every zero-row must have integer θ. This swap channel is essential: e.g. for
#      sg 27 (Pcc2) at idxmax=(1,1,1) only one swapped pair survives the glide
#      extinctions and the generic lattice is exactly centrosymmetric about a
#      modulation-dependent center (realizing Pccm) — a symmetry invisible to any
#      criterion that demands a modulation-independent w. Given (T), candidates with
#      W ∈ P(G) contribute nothing new (composing with the group's own {W|w₀} reduces
#      them to pure translations, which have no draw-dependence).
#
# Monotonicity: enlarging the box only adds orbits, which can only shrink the generic
# stabilizer; hence sufficiency is upward-closed in the componentwise partial order and
# both the minimal isotropic N and the antichain of minimal anisotropic tuples are
# well-defined.
#
# Usage:  julia --project=<repo> [-t N] build/compute_minimal_fourier_truncations.jl \
#             [--emit-table=PATH] [--verify] [--outdir=PATH] \
#             [--dims=1,2,3] [--smoke] [--isocap=8] [--verify-extent=7]
#
# To regenerate the tabulated data shipped with the package (the canonical invocation):
#   julia --project=. build/compute_minimal_fourier_truncations.jl \
#       --verify --emit-table=src/tables/fourier_truncations.jl
#
# `--emit-table` writes the antichains of minimal symmetry-safe truncations as a Julia
# source file (refusing to do so if any group failed validation). `--verify` additionally
# checks that each antichain is complete over the grid {1..7}ᴰ. `--outdir` is optional and
# writes the *audit* artifacts — `{1,2,3}d.csv` (per-group diagnostics: minimal
# truncations, the accidental symmetry operations arising at the previous default
# `(2,…,2)`, orbit counts, validation status) and a human-readable `summary.md`. Those
# artifacts are not needed to build the package; they exist to explain *why* each tabulated
# entry is what it is. With neither flag, the run only validates and reports.

using Crystalline
using Crystalline: UnityFourierLattice, getorbits, getcoefs
using Crystalline.SmithNormalForm: smith
using Bravais: all_centeringtranslations, centering_volume_fraction
using LinearAlgebra, StaticArrays, Random, Printf

# ---------------------------------------------------------------------------------------
# Exact integer-lattice utilities (row-lattice reduction, membership, congruences)
# ---------------------------------------------------------------------------------------

_pivot(v::Vector{Int}) = findfirst(!=(0), v)

"""
Insert integer row `v` into `rows`, a generating set of a row lattice kept with distinct
pivot (first-nonzero) columns. The lattice generated by `rows` grows to include `v`.
"""
function _addrow!(rows::Vector{Vector{Int}}, v::Vector{Int})
    j = _pivot(v)
    while j !== nothing
        maximum(abs, v) < 10^14 || error("integer overflow risk in row reduction")
        k = findfirst(r -> _pivot(r) == j, rows)
        if k === nothing
            push!(rows, v)
            sort!(rows; by=_pivot)
            return rows
        end
        b = rows[k]
        g, x, y = gcdx(b[j], v[j])
        newb = x .* b .+ y .* v          # pivot j, value g (entries before j vanish)
        newv = (b[j] ÷ g) .* v .- (v[j] ÷ g) .* b  # pivot strictly beyond j
        rows[k] = newb
        v = newv
        j = _pivot(v)
    end
    return rows
end

"""
Exact membership test: is integer vector `b` in the row lattice generated by `rows`
(as maintained by `_addrow!`, i.e. distinct ascending pivots)?
"""
function _inlattice(rows::Vector{Vector{Int}}, b::Vector{Int})
    v = copy(b)
    for r in rows # rows sorted by ascending pivot
        j = _pivot(r)
        v[j] == 0 && continue
        v[j] % r[j] == 0 || return false
        v .-= (v[j] ÷ r[j]) .* r
    end
    return all(==(0), v)
end

"""
Decide solvability of the congruence system  a·w ≡ θₐ (mod 1)  over w ∈ ℝᴰ, given as
augmented integer rows [a₁,…,a_D, qθ] (common denominator `q`). Returns `nothing` if
unsolvable, else a particular rational solution `w::Vector{Rational{Int}}` (components
reduced mod 1; free components set to 0).

Method: integer row reduction of the augmented rows is a unimodular re-combination of the
constraints, hence solution-set preserving. In the reduced system the rows have distinct
ascending pivots; a row with zero G-part (pivot in the θ-column) is the constraint
0 ≡ t/q (mod 1), i.e. q | t — the only possible obstruction. The remaining rows are
solved exactly by back-substitution (constraints imposed with equality, hence mod 1).
"""
function _solve_congruences(augrows::Vector{Vector{Int}}, q::Int, D::Int)
    red = Vector{Int}[]
    for v in augrows
        _addrow!(red, copy(v))
    end
    for r in red
        if _pivot(r) == D + 1 # zero G-part: 0·w ≡ r[D+1]/q (mod 1)
            r[D+1] % q == 0 || return nothing
        end
    end
    w = zeros(Rational{Int}, D)
    for r in Iterators.reverse(red) # ascending pivots ⇒ reverse = back-substitution
        j = _pivot(r)
        j ≤ D || continue
        acc = r[D+1] // q
        for k in j+1:D
            acc -= r[k] * w[k]
        end
        w[j] = acc // r[j]
    end
    return mod.(w, 1)
end

"""
Float-rhs variant of `_solve_congruences` for a *specific* modulation draw: solve
a·w ≡ rhs (mod 1) with integer rows `As` and real right-hand sides. Returns `w` or
`nothing` if inconsistent (within `tol`). Used to cross-validate the exact generic
decision against an actual realization.
"""
function _solve_congruences_float(As::Vector{Vector{Int}}, rhs::Vector{Float64}, D::Int;
                                  tol::Real=1e-6)
    red = Tuple{Vector{Int},Float64}[]
    for (a0, r0) in zip(As, rhs)
        a, r = copy(a0), r0
        inserted = false
        j = _pivot(a)
        while j !== nothing
            k = findfirst(t -> _pivot(t[1]) == j, red)
            if k === nothing
                push!(red, (a, r))
                sort!(red; by=t -> _pivot(t[1]))
                inserted = true
                break
            end
            b, rb = red[k]
            g, x, y = gcdx(b[j], a[j])
            red[k] = (x .* b .+ y .* a, x * rb + y * r)
            a, r = (b[j] ÷ g) .* a .- (a[j] ÷ g) .* b, (b[j] ÷ g) * r - (a[j] ÷ g) * rb
            j = _pivot(a)
        end
        if !inserted # row reduced to zero: consistent only if rhs is (near-)integer
            abs(r - round(r)) < tol || return nothing
        end
    end
    w = zeros(D)
    for (a, r) in Iterators.reverse(red)
        j = _pivot(a)
        acc = r
        for k in j+1:D
            acc -= a[k] * w[k]
        end
        w[j] = acc / a[j]
    end
    return mod.(w, 1)
end

# ---------------------------------------------------------------------------------------
# Holohedry of the Bravais type (conventional setting) and expected reciprocal lattice
# ---------------------------------------------------------------------------------------

_metric(Rs) = (R = stack(Rs); R' * R)

"""
All rotation parts (integer matrices, conventional basis) of isometries of a lattice with
*generic* parameters of `sgnum`'s Bravais class: entries in {-1,0,1}, preserving three
independently sampled generic metrics of the crystal family and mapping the centered
translation lattice to itself. Deterministic per Bravais type (cached by caller).
"""
function _holohedry_rotations(sgnum::Integer, Dᵛ::Val{D}) where D
    cntr = centering(sgnum, D)
    Tc = all_centeringtranslations(cntr, Dᵛ)
    gs = [_metric(directbasis(sgnum, Dᵛ)) for _ in 1:3]
    Ws = SMatrix{D,D,Int,D*D}[]
    for t in Iterators.product(ntuple(_ -> -1:1, D * D)...)
        W = SMatrix{D,D,Int}(t)
        abs(round(Int, det(SMatrix{D,D,Float64}(W)))) == 1 || continue
        all(g -> norm(W' * g * W - g) < 1e-8 * norm(g), gs) || continue
        # W must map the centered lattice to itself: W·t ∈ ℤᴰ + {centering translations}
        ok = all(Tc) do tc
            wt = W * tc
            any(tc′ -> all(x -> abs(x - round(x)) < 1e-9, wt - tc′),
                Iterators.flatten(((zero(wt),), Tc)))
        end
        ok && push!(Ws, W)
    end
    return Ws
end

const HOLOHEDRY_CACHE = Dict{Tuple{Int,String},Vector}()
function holohedry_rotations(sgnum::Integer, Dᵛ::Val{D}) where D
    bt = bravaistype(sgnum, D)
    return get!(HOLOHEDRY_CACHE, (D, bt)) do
        Ws = _holohedry_rotations(sgnum, Dᵛ)
        length(Ws) == expected_holohedry_order(D, bt) ||
            error("holohedry order mismatch for ($D, $bt): got $(length(Ws))")
        Ws
    end :: Vector{SMatrix{D,D,Int,D*D}}
end

function expected_holohedry_order(D::Integer, bt::String)
    f = first(bt)
    if D == 1
        return 2
    elseif D == 2
        return f == 'm' ? 2 : f == 'o' ? 4 : f == 't' ? 8 : f == 'h' ? 12 :
               error("unknown 2D bravais type $bt")
    else
        return f == 'a' ? 2 : f == 'm' ? 4 : f == 'o' ? 8 : f == 't' ? 16 :
               f == 'h' ? (bt == "hR" ? 12 : 24) : f == 'c' ? 48 :
               error("unknown 3D bravais type $bt")
    end
end

"""
Integer basis (columns) of L_exp = {m ∈ ℤᴰ : m·t ∈ ℤ ∀ centering translations t}, the
reciprocal lattice of the centered direct lattice in conventional reciprocal coordinates.
"""
function Lexp_basis(cntr::Char, Dᵛ::Val{D}) where D
    Tc = all_centeringtranslations(cntr, Dᵛ)
    isempty(Tc) && return Matrix{Int}(I, D, D)
    q = 12
    C = zeros(Int, length(Tc), D)
    for (i, t) in enumerate(Tc)
        qt = q .* t
        all(x -> abs(x - round(x)) < 1e-9, qt) || error("unexpected centering translation $t")
        C[i, :] .= round.(Int, qt)
    end
    # m obeys C m ≡ 0 (mod q). With C = SΛT (Smith), y = Tm: λᵢyᵢ ≡ 0 (mod q), i.e.
    # yᵢ ∈ dᵢℤ with dᵢ = q ÷ gcd(λᵢ, q); basis of L_exp: columns dᵢ·T⁻¹[:,i].
    F = smith(C)
    d = [q ÷ gcd(i ≤ length(F.SNF) ? F.SNF[i] : 0, q) for i in 1:D]
    B = F.Tinv * Diagonal(d)
    abs(round(Int, det(Matrix{Float64}(B)))) == centering_volume_fraction(cntr, Dᵛ) ||
        error("L_exp index mismatch for centering $cntr")
    return B
end

# ---------------------------------------------------------------------------------------
# Master orbit data (computed once per (sgnum, cap); sub-boxes obtained by filtering)
# ---------------------------------------------------------------------------------------

struct Master{D}
    cap        :: Int
    orbits     :: Vector{Vector{SVector{D,Int}}}
    coefs      :: Vector{Vector{ComplexF64}}
    dicts      :: Vector{Dict{SVector{D,Int},Int}}   # member -> intra-orbit position
    whichorbit :: Dict{SVector{D,Int},Int}           # member -> orbit index
end

function Master(sgnum::Integer, Dᵛ::Val{D}, cap::Integer) where D
    flat = levelsetlattice(sgnum, Dᵛ, ntuple(_ -> Int(cap), Val(D)))
    orbits = getorbits(flat)
    coefs = getcoefs(flat)
    dicts = [Dict(G => i for (i, G) in enumerate(orb)) for orb in orbits]
    whichorbit = Dict(G => i for (i, orb) in enumerate(orbits) for G in orb)
    # realness closure: the negative partner orbit -O must be retained whenever O is
    for orb in orbits
        haskey(whichorbit, -first(orb)) ||
            error("orbit list not closed under negation (sg $sgnum, cap $cap)")
    end
    return Master{D}(cap, orbits, coefs, dicts, whichorbit)
end

# orbit indices retained for a sub-box `idxmax` (an orbit is seeded iff it meets the box)
function retained(mst::Master{D}, idxmax::NTuple{D,Int}) where D
    return [i for (i, orb) in enumerate(mst.orbits) if
            any(G -> all(abs.(G) .≤ idxmax), orb)]
end

# ---------------------------------------------------------------------------------------
# The sufficiency check
# ---------------------------------------------------------------------------------------

struct ExtraOp{D,L}
    W     :: SMatrix{D,D,Int,L}
    w     :: Union{Nothing,SVector{D,Float64}} # modulation-independent translation (nswap = 0)
    nswap :: Int                               # number of swapped orbit pairs
end

function oplabel(e::ExtraOp{D}) where D
    if e.w !== nothing
        op = SymOperation(SMatrix{D,D,Float64}(e.W), e.w)
        return try seitz(op) catch; xyzt(op) end
    else # translation part depends on the modulation draw
        op = SymOperation(SMatrix{D,D,Float64}(e.W), zero(SVector{D,Float64}))
        lab = try seitz(op) catch; xyzt(op) end
        return lab * "|w(c)"
    end
end

struct ExcessResult{D,L}
    sufficient   :: Bool
    spanok       :: Bool
    extra_transl :: Union{Nothing,SVector{D,Float64}} # representative, if any
    extra_ops    :: Vector{ExtraOp{D,L}}
    nretained    :: Int
end

"""
Collect the phase-congruence system for candidate rotation `Wit = W⁻ᵀ` on the retained
orbits `sel` of `mst`: rows  Gᵣ·w - s_{srow_r} ≡ θᵣ (mod 1), with `srow_r = 0` for
orbit-preserved rows and `srow_r = p > 0` for rows of swapped pair `p` (whose generic
phase s_p = arg(z_p)/π is set by the modulation). Returns `nothing` if `W` is
structurally excluded (some orbit image is neither the orbit itself nor its negative
partner, or not retained). `pairreps[p]` is the master orbit index defining z_p.
"""
function _candidate_system(mst::Master{D}, sel, Wit) where D
    selset = Set(sel)
    Grows = SVector{D,Int}[]
    srows = Int[]
    θs = Float64[]
    pairreps = Int[]
    for i in sel
        orb = mst.orbits[i]
        iszero(first(orb)) && continue # the G = 0 orbit is trivially symmetric
        v = mst.coefs[i]
        ι = mst.whichorbit[-first(orb)] # negative-partner orbit (= i if self-negating)
        j = get(mst.whichorbit, Wit * first(orb), 0)
        j in selset || return nothing
        if j == i # orbit preserved setwise: modulation-independent phase conditions
            dict = mst.dicts[i]
            for (m, G) in enumerate(orb)
                G′ = Wit * G
                n = get(dict, G′, 0)
                n == 0 && return nothing
                push!(Grows, G′); push!(srows, 0)
                push!(θs, angle(v[m] / v[n]) / (2π))
            end
        elseif j == ι && ι != i # swapped pair {O, -O}; conditions carry generic s_p
            i > ι && continue   # already handled from the partner (rows are conjugate-redundant)
            dictι = mst.dicts[ι]
            dicti = mst.dicts[i]
            push!(pairreps, i)
            p = length(pairreps)
            for (m, G) in enumerate(orb)
                G′ = Wit * G
                get(dictι, G′, 0) == 0 && return nothing # must land in the partner orbit
                n = get(dicti, -G′, 0)
                n == 0 && error("negation closure violated within orbit pair")
                push!(Grows, G′); push!(srows, p)
                push!(θs, angle(v[m] * v[n]) / (2π))
            end
        else
            return nothing # image is an unrelated orbit: |z_p| = |z_q| impossible generically
        end
    end
    return (; Grows, srows, θs, npairs=length(pairreps), pairreps)
end

"""
Exact decision: is the system from `_candidate_system` solvable in `w` for *generic*
swapped-pair phases `s`? If yes with no swapped pairs, also return the (modulation-
independent) representative translation.

Method: eliminate each s_p analytically. Rows of pair p are  G′ₘ·w - s_p ≡ θₘ; their
intra-pair differences (G′ₘ - G′₁)·w ≡ θₘ - θ₁ are s-free and join the preserved-orbit
rows in a standard congruence system (a). The remaining per-pair conditions
g_p·w ≡ θ₁ + s_p (g_p := G′₁ of pair p) are satisfiable for generic independent s_p iff
no integer combination Σnₚg_p with n ≠ 0 lies in the ℤ-span of the s-free rows' G-parts
(such a combination would pin Σnₚs_p — measure zero); since any rational dependence can
be scaled to such an integer combination, this holds iff the images of the g_p in
ℚᴰ/span_ℚ(free rows) are ℚ-linearly independent (b).
"""
function _decide_generic(sys, ::Val{D}) where D
    P = sys.npairs
    θrs = map(sys.θs) do θ
        θr = rationalize(Int, θ; tol=1e-7)
        abs(θ - θr) < 1e-7 || error("phase not rational: θ = $θ")
        denominator(θr) ≤ 48 || error("unexpectedly large phase denominator: $θr")
        θr
    end
    if P == 0
        q = reduce(lcm, denominator.(θrs); init=1)
        aug = [vcat(collect(G), Int(θ * q)) for (G, θ) in zip(sys.Grows, θrs)]
        w = _solve_congruences(aug, q, D)
        return w === nothing ? (false, nothing) : (true, SVector{D,Float64}(float.(w)))
    end
    # split: s-free rows (preserved orbits + intra-pair differences) and pair reps g_p
    firstrow = Dict{Int,Int}()
    freeG = Vector{Int}[]
    freeθ = Rational{Int}[]
    gs = Vector{SVector{D,Int}}(undef, P)
    for (k, p) in enumerate(sys.srows)
        if p == 0
            push!(freeG, collect(sys.Grows[k]))
            push!(freeθ, θrs[k])
        elseif !haskey(firstrow, p)
            firstrow[p] = k
            gs[p] = sys.Grows[k]
        else
            k0 = firstrow[p]
            push!(freeG, collect(sys.Grows[k] - sys.Grows[k0]))
            push!(freeθ, θrs[k] - θrs[k0])
        end
    end
    # (a) the modulation-independent subsystem must be solvable
    qf = reduce(lcm, denominator.(freeθ); init=1)
    aug = [vcat(a, Int(θ * qf)) for (a, θ) in zip(freeG, freeθ)]
    _solve_congruences(aug, qf, D) === nothing && return (false, nothing)
    # (b) generic phases must be absorbable by w
    BL = Vector{Int}[]
    for a in freeG
        any(!=(0), a) && _addrow!(BL, copy(a))
    end
    r = length(BL)
    r + P > D && return (false, nothing)
    M = Matrix{Rational{Int}}(undef, D, r + P)
    for (j, a) in enumerate(BL)
        M[:, j] .= a
    end
    for p in 1:P
        M[:, r+p] .= gs[p]
    end
    _qrank(M) == r + P || return (false, nothing)
    return (true, nothing)
end

function _qrank(M::AbstractMatrix{<:Rational})
    A = Matrix{Rational{BigInt}}(M)
    nr, nc = size(A)
    rk = 0
    for j in 1:nc
        i0 = findnext(i -> !iszero(A[i, j]), rk+1:nr, 1)
        i0 === nothing && continue
        i0 = (rk+1:nr)[i0]
        if i0 != rk + 1
            A[rk+1, :], A[i0, :] = A[i0, :], A[rk+1, :]
        end
        rk += 1
        for i in rk+1:nr
            iszero(A[i, j]) && continue
            A[i, :] .-= (A[i, j] // A[rk, j]) .* A[rk, :]
        end
    end
    return rk
end

function excess_symmetry(mst::Master{D}, hol, pgrots, Bexp, Tc,
                         idxmax::NTuple{D,Int}) where D
    sel = retained(mst, idxmax)

    # --- (T): ℤ-span of retained G's must equal L_exp ------------------------------------
    rowsT = Vector{Int}[]
    for i in sel, G in mst.orbits[i]
        iszero(G) && continue
        _addrow!(rowsT, collect(G))
    end
    spanok = length(rowsT) == D && all(j -> _inlattice(rowsT, Bexp[:, j]), 1:D)
    extra_transl = spanok ? nothing : _extra_translation(rowsT, Tc, Val(D))

    # --- (R): extra holohedry operations (draw-dependent translations allowed) -----------
    extra_ops = ExtraOp{D,D*D}[]
    for W in hol
        W in pgrots && continue
        Wi = round.(Int, inv(SMatrix{D,D,Float64}(W)))
        W * Wi == SMatrix{D,D,Int}(I) || error("non-unimodular holohedry candidate")
        Wit = SMatrix{D,D,Int}(transpose(Wi))
        sys = _candidate_system(mst, sel, Wit)
        sys === nothing && continue
        isextra, w = _decide_generic(sys, Val(D))
        isextra && push!(extra_ops, ExtraOp{D,D*D}(W, w, sys.npairs))
    end

    sufficient = spanok && isempty(extra_ops)
    return ExcessResult{D,D*D}(sufficient, spanok, extra_transl, extra_ops, length(sel))
end

# representative translation w with exp(2πi G·w) = 1 for all retained G but w ∉ ℤᴰ + Tc
function _extra_translation(rowsT::Vector{Vector{Int}}, Tc, ::Val{D}) where D
    pivots = Set(_pivot.(rowsT))
    for j in 1:D
        j ∉ pivots && return SVector{D,Float64}(ntuple(k -> k == j ? 0.437 : 0.0, D))
    end
    # full rank: dual basis vectors of span(L); at least one lies outside ℤᴰ + Tc
    B = Rational{Int}.(reduce(hcat, rowsT)') # rows = lattice basis
    for i in 1:D
        d = B \ [Rational{Int}(k == i) for k in 1:D]
        dm = mod.(float.(d), 1)
        inlat = any(tc′ -> all(x -> abs(x - round(x)) < 1e-9, dm - tc′),
                    Iterators.flatten((((@SVector zeros(D)),), Tc)))
        inlat || return SVector{D,Float64}(dm)
    end
    error("inconsistency: span deficient but all dual vectors are lattice vectors")
end

# ---------------------------------------------------------------------------------------
# Minimal truncation search
# ---------------------------------------------------------------------------------------

mutable struct GroupContext{D,L} # L = D*D (StaticArrays length parameter)
    sgnum   :: Int
    Dᵛ      :: Val{D}
    hol     :: Vector{SMatrix{D,D,Int,L}}
    pgrots  :: Set{SMatrix{D,D,Int,L}}
    Bexp    :: Matrix{Int}
    Tc      :: Vector{SVector{D,Float64}}
    mst     :: Master{D}
    memo    :: Dict{NTuple{D,Int},ExcessResult{D,L}}
end

function GroupContext(sgnum::Integer, Dᵛ::Val{D}; cap::Integer=4) where D
    sg = spacegroup(sgnum, Dᵛ)
    pgrots = Set(SMatrix{D,D,Int,D*D}(round.(Int, Matrix(rotation(op)))) for op in sg)
    hol = holohedry_rotations(sgnum, Dᵛ)
    cntr = centering(sgnum, D)
    Tc = all_centeringtranslations(cntr, Dᵛ)
    Bexp = Lexp_basis(cntr, Dᵛ)
    mst = Master(sgnum, Dᵛ, cap)
    return GroupContext{D,D*D}(sgnum, Dᵛ, hol, pgrots, Bexp, Tc, mst,
                               Dict{NTuple{D,Int},ExcessResult{D,D*D}}())
end

function ensure_cap!(ctx::GroupContext{D}, cap::Integer) where D
    if cap > ctx.mst.cap
        ctx.mst = Master(ctx.sgnum, ctx.Dᵛ, cap)
        # memoized results for boxes within the old cap remain valid (orbit filtering
        # yields identical retained orbit sets); no need to clear `memo`
    end
    return ctx
end

function check(ctx::GroupContext{D}, t::NTuple{D,Int}) where D
    maximum(t) ≤ ctx.mst.cap || error("box exceeds master cap")
    r = get!(ctx.memo, t) do
        excess_symmetry(ctx.mst, ctx.hol, ctx.pgrots, ctx.Bexp, ctx.Tc, t)
    end
    return r
end

function minimal_isotropic(ctx::GroupContext{D}; isocap::Integer=8) where D
    for N in 1:isocap
        ensure_cap!(ctx, max(N, ctx.mst.cap))
        check(ctx, ntuple(_ -> N, Val(D))).sufficient && return N
    end
    return 0 # not found within cap
end

function minimal_anisotropic(ctx::GroupContext{D}, Niso::Integer) where D
    Niso == 0 && return NTuple{D,Int}[]
    Ncap = Niso + 2
    while true
        ensure_cap!(ctx, Ncap)
        suff = Dict{NTuple{D,Int},Bool}()
        for t in sort!(vec(collect(Iterators.product(ntuple(_ -> 1:Ncap, Val(D))...)));
                       by=sum)
            # monotonic inference: any sufficient tuple ≤ t makes t sufficient
            if any(k -> t[k] > 1 && suff[ntuple(j -> j == k ? t[j] - 1 : t[j], Val(D))],
                   1:D)
                suff[t] = true
            else
                suff[t] = check(ctx, t).sufficient
            end
        end
        minimal = [t for (t, s) in suff if s &&
                   all(k -> t[k] == 1 ||
                            !suff[ntuple(j -> j == k ? t[j] - 1 : t[j], Val(D))], 1:D)]
        sort!(minimal)
        if any(t -> maximum(t) == Ncap, minimal) && Ncap < 12
            Ncap += 2 # a minimal element touches the grid boundary; enlarge and redo
            continue
        end
        return minimal
    end
end

"""
Verify that the antichain `aniso` is *complete* on the grid `{1..extent}ᴰ`, i.e. that every
sufficient tuple there dominates some element of `aniso`. Because sufficiency is
upward-closed, completeness of the antichain is exactly what makes the dominance test
`is_symmetry_safe` exact rather than conservative; `minimal_anisotropic` only searches a
grid of extent `Nmin+2`, so this checks well beyond it. Returns the list of missed
sufficient tuples (empty if complete).
"""
function verify_antichain_complete(ctx::GroupContext{D}, aniso::Vector{NTuple{D,Int}};
                                   extent::Integer=7) where D
    isempty(aniso) && return NTuple{D,Int}[]
    ensure_cap!(ctx, extent)
    missed = NTuple{D,Int}[]
    for t in Iterators.product(ntuple(_ -> 1:Int(extent), Val(D))...)
        any(p -> all(k -> p[k] ≤ t[k], ntuple(identity, Val(D))), aniso) && continue
        check(ctx, t).sufficient && push!(missed, t)
    end
    return missed
end

# ---------------------------------------------------------------------------------------
# Numerical cross-validation (algebra vs. direct evaluation of modulated lattices)
# ---------------------------------------------------------------------------------------

function _filtered_unity(mst::Master{D}, idxmax::NTuple{D,Int}) where D
    sel = retained(mst, idxmax)
    return UnityFourierLattice{D}(mst.orbits[sel], mst.coefs[sel])
end

function _op_residual(mflat, op::SymOperation{D}, rs) where D
    W = SMatrix{D,D,Float64}(Matrix(rotation(op)))
    w = SVector{D,Float64}(translation(op))
    return maximum(abs(mflat(W * r + w) - mflat(r)) for r in rs)
end

"""
Numerically validate the algebraic verdicts for `sgnum` at box `idxmax` against a
concrete random modulation:
- the group's own operations and any claimed extra translation must hold;
- every accepted extra rotation must hold with its translation part *solved from the
  draw* (via the swapped-pair phases s_p = arg(z_p)/π);
- every holohedry candidate that passed the structural stage but was rejected by the
  generic phase decision must admit NO translation part for this draw either (the
  float congruence system must be inconsistent).
Returns `(ok, msgs)`.
"""
function validate_numeric(ctx::GroupContext{D}, idxmax::NTuple{D,Int};
                          npts::Integer=120) where D
    res = check(ctx, idxmax)
    sel = retained(ctx.mst, idxmax)
    rng = Xoshiro(hash((ctx.sgnum, D, idxmax)))
    flatU = _filtered_unity(ctx.mst, idxmax)
    N = length(getcoefs(flatU))
    modulation = (0.5 .+ 0.5 .* rand(rng, N)) .* cis.(2π .* rand(rng, N))
    mflat = modulate(flatU, modulation)
    rs = [rand(rng, SVector{D,Float64}) for _ in 1:npts]
    scale = maximum(abs, mflat.(rs)) + 1e-3
    msgs = String[]

    # effective real-field coefficients a_G of the draw (Re folds in ±G conjugates)
    aG = Dict{SVector{D,Int},ComplexF64}()
    for (orb, cs) in zip(getorbits(mflat), getcoefs(mflat)), (G, c) in zip(orb, cs)
        aG[G] = get(aG, G, zero(ComplexF64)) + c / 2
        aG[-G] = get(aG, -G, zero(ComplexF64)) + conj(c) / 2
    end

    for op in spacegroup(ctx.sgnum, ctx.Dᵛ) # sanity: target group must be realized
        ρ = _op_residual(mflat, op, rs)
        ρ < 1e-8 * scale || push!(msgs, "group op $(seitz(op)) broken (ρ=$ρ)")
    end
    if res.extra_transl !== nothing # claimed extra translation must hold exactly
        ρ = maximum(abs(mflat(r + res.extra_transl) - mflat(r)) for r in rs)
        ρ < 1e-8 * scale || push!(msgs, "claimed extra translation fails (ρ=$ρ)")
    end

    accepted = Dict(e.W => e for e in res.extra_ops)
    for W in ctx.hol
        W in ctx.pgrots && continue
        Wi = round.(Int, inv(SMatrix{D,D,Float64}(W)))
        Wit = SMatrix{D,D,Int}(transpose(Wi))
        sys = _candidate_system(ctx.mst, sel, Wit)
        if sys === nothing # structurally excluded; support mismatch — spot check only
            haskey(accepted, W) &&
                push!(msgs, "structurally excluded op was accepted: $(oplabel(accepted[W]))")
            op = SymOperation(SMatrix{D,D,Float64}(W), rand(rng, SVector{D,Float64}))
            ρ = _op_residual(mflat, op, rs)
            ρ > 1e-6 * scale ||
                push!(msgs, "structurally excluded W numerically unbroken (ρ=$ρ)")
            continue
        end
        # draw-dependent phases: s_p = arg(z_p)/π with z_p = a_{G₁}/v_{i,1} at pair rep
        svals = map(sys.pairreps) do i
            G₁ = first(ctx.mst.orbits[i])
            angle(aG[G₁] / first(ctx.mst.coefs[i])) / π
        end
        rhs = [θ + (p == 0 ? 0.0 : svals[p]) for (θ, p) in zip(sys.θs, sys.srows)]
        wf = _solve_congruences_float([collect(G) for G in sys.Grows], rhs, D)
        if haskey(accepted, W)
            if wf === nothing
                push!(msgs, "accepted op $(oplabel(accepted[W])) has no per-draw translation")
            else
                op = SymOperation(SMatrix{D,D,Float64}(W), SVector{D,Float64}(wf))
                ρ = _op_residual(mflat, op, rs)
                ρ < 1e-8 * scale ||
                    push!(msgs, "accepted op $(oplabel(accepted[W])) fails numerically (ρ=$ρ)")
            end
        else # rejected by the generic phase decision: this draw must not admit it either
            if wf !== nothing
                op = SymOperation(SMatrix{D,D,Float64}(W), SVector{D,Float64}(wf))
                ρ = _op_residual(mflat, op, rs)
                ρ > 1e-6 * scale ||
                    push!(msgs, "rejected op IS a symmetry of the draw: $(seitz(op)) (ρ=$ρ)")
            end
        end
    end
    return isempty(msgs), msgs
end

# ---------------------------------------------------------------------------------------
# Self-tests (cheap; run at startup)
# ---------------------------------------------------------------------------------------

function selftests()
    # Smith normal form convention: A = S·Λ·T
    A = [2 4 4; -6 6 12; 10 4 16]
    F = smith(A)
    Λ = zeros(Int, size(A))
    for (i, λ) in enumerate(F.SNF)
        Λ[i, i] = λ
    end
    F.S * Λ * F.T == A || error("unexpected Smith normal form convention")

    # row-lattice utilities
    rows = Vector{Int}[]
    _addrow!(rows, [2, 0]); _addrow!(rows, [0, 3])
    (_inlattice(rows, [4, 3]) && !_inlattice(rows, [1, 0])) ||
        error("row-lattice membership self-test failed")

    # congruence solver: w ≡ θ₁/2 (mod 1/2) and w ≡ θ₂/3 (mod 1/3) solvable iff
    # 3θ₁ - 2θ₂ ∈ ℤ  (here θ = t/q with q = 12)
    sol = _solve_congruences([[2, 3], [3, 8]], 12, 1)   # θ₁ = 1/4, θ₂ = 2/3: 3/4-4/3 ∉ ℤ
    sol === nothing || error("congruence self-test (unsolvable case) failed")
    sol = _solve_congruences([[2, 3], [3, 14]], 12, 1)  # θ₁ = 1/4, θ₂ = 7/6: 3/4-7/3 ∉ ℤ
    sol === nothing || error("congruence self-test (unsolvable case 2) failed")
    sol = _solve_congruences([[2, 6], [3, 9]], 12, 1)   # θ₁ = 1/2, θ₂ = 3/4: 3/2-3/2 ∈ ℤ
    sol !== nothing || error("congruence self-test (solvable case) failed")
    (mod(2 * sol[1], 1) == mod(1//2, 1) && mod(3 * sol[1], 1) == mod(3//4, 1)) ||
        error("congruence self-test solution incorrect")

    # holohedry orders for a representative set (also exercises the cache)
    for (D, sgnums) in ((1, (1, 2)),
                        (2, (1, 3, 5, 10, 13)),
                        (3, (1, 3, 5, 16, 20, 22, 23, 38, 75, 79, 143, 146, 195, 196, 197)))
        for sgnum in sgnums
            holohedry_rotations(sgnum, Val(D)) # errors internally on order mismatch
        end
    end

    # sg 230 ground truth: exactly 3 orbits retained at the default (2,2,2)
    mst230 = Master(230, Val(3), 2)
    length(retained(mst230, (2, 2, 2))) == 3 ||
        error("sg 230 (2,2,2) orbit-count self-test failed")

    # --- swap-channel (draw-dependent center) anchors ---
    # 1D p1: a single retained ±G pair is a pure cosine, always centrosymmetric about a
    # modulation-dependent center; two harmonics generically are not → minimal N = 2
    minimal_isotropic(GroupContext(1, Val(1))) == 2 ||
        error("1D p1 swap-channel anchor failed (expected minimal N = 2)")
    # 3D P1 at (1,1,1): 13 swapped pairs with ℤ-relations among their G's pin the generic
    # phases → no accidental inversion → sufficient
    check(GroupContext(1, Val(3)), (1, 1, 1)).sufficient ||
        error("3D P1 anchor failed (expected (1,1,1) sufficient)")
    # sg 27 (Pcc2) at (1,1,1): glide extinctions leave a single swapped pair; inversion
    # with a modulation-dependent center is a generic symmetry → insufficient
    r27 = check(GroupContext(27, Val(3)), (1, 1, 1))
    (!r27.sufficient &&
     any(e -> e.nswap > 0 && e.W == -SMatrix{3,3,Int}(I), r27.extra_ops)) ||
        error("sg 27 swap-channel anchor failed (expected inversion with w(c))")

    return nothing
end

# ---------------------------------------------------------------------------------------
# Per-group driver, CSV/summary output, main
# ---------------------------------------------------------------------------------------

function process_group(sgnum::Integer, Dᵛ::Val{D}; isocap::Integer=8,
                       validate::Bool=false, verify::Bool=false,
                       verify_extent::Integer=7) where D
    ctx = GroupContext(sgnum, Dᵛ)
    default = ntuple(_ -> 2, Val(D))
    dres = check(ctx, default)
    Niso = minimal_isotropic(ctx; isocap)
    aniso = minimal_anisotropic(ctx, Niso)

    completeness_msg = ""
    if verify
        missed = verify_antichain_complete(ctx, aniso; extent=verify_extent)
        isempty(missed) ||
            (completeness_msg = "antichain incomplete: missed $(missed)")
    end

    dovalidate = validate || !dres.sufficient
    valstatus = :skipped
    valmsgs = String[]
    if dovalidate
        okd, msgs1 = validate_numeric(ctx, default)
        append!(valmsgs, msgs1)
        okm = true
        if Niso != 0
            okm, msgs2 = validate_numeric(ctx, ntuple(_ -> Niso, Val(D)))
            append!(valmsgs, msgs2)
        end
        # monotonicity spot-check at the isotropic threshold
        if Niso != 0 && Niso + 1 ≤ 12
            ensure_cap!(ctx, Niso + 1)
            check(ctx, ntuple(_ -> Niso + 1, Val(D))).sufficient ||
                push!(valmsgs, "monotonicity violated at N = $(Niso + 1)")
        end
        valstatus = (okd && okm && isempty(valmsgs)) ? :ok : :failed
    end
    if !isempty(completeness_msg)
        push!(valmsgs, completeness_msg)
        valstatus = :failed
    end

    return (; D, sgnum,
              label = iuc(sgnum, D),
              bt = bravaistype(sgnum, D),
              norb_default = dres.nretained,
              default_ok = dres.sufficient,
              Niso, aniso,
              extra_ops_default = [oplabel(e) for e in dres.extra_ops],
              extra_transl_default = dres.extra_transl === nothing ? "" :
                                     string(round.(dres.extra_transl; digits=4)),
              valstatus,
              valmsg = join(valmsgs, " / "),
              err = "")
end

_csvq(x) = (s = string(x); '"' * replace(s, '"' => "\"\"") * '"')

function write_csv(path::String, rows)
    open(path, "w") do io
        println(io, join(["D", "sgnum", "label", "bravaistype", "norbits_default",
                          "default_sufficient", "min_isotropic",
                          "min_anisotropic", "extra_ops_at_default",
                          "extra_translation_at_default", "validation", "notes"], ','))
        for r in rows
            println(io, join([r.D, r.sgnum, _csvq(r.label), r.bt, r.norb_default,
                              r.default_ok, r.Niso,
                              _csvq(join(("(" * join(t, ',') * ")" for t in r.aniso), "; ")),
                              _csvq(join(r.extra_ops_default, "; ")),
                              _csvq(r.extra_transl_default), r.valstatus,
                              _csvq(r.valmsg * (isempty(r.err) ? "" : " ERROR: " * r.err))],
                             ','))
        end
    end
end

function write_summary(path::String, allrows)
    open(path, "w") do io
        println(io, "# Minimal Fourier truncations for `levelsetlattice`\n")
        println(io, """
        A truncation `idxmax` is *sufficient* when a generically modulated
        `levelsetlattice(sgnum, Val(D), idxmax)` has exactly the symmetry of the target
        group — no extra point-group operations and no extra periodicity. Extra operations
        whose translation part depends on the modulation draw (e.g. an inversion center
        that shifts with the random coefficients) count as accidental symmetry and are
        labeled `…|w(c)`. Since orbits are completed beyond the seed box, `idxmax` here
        means "all orbits meeting the box |Gᵢ| ≤ idxmaxᵢ".
        `min_isotropic` is the smallest sufficient N = idxmax₁ = … ;
        `min_anisotropic` lists the minimal sufficient tuples (an antichain — no component
        of a listed tuple can be lowered without losing sufficiency).\n""")
        for D in sort!(unique(r.D for r in allrows))
            rows = sort!([r for r in allrows if r.D == D]; by=r -> r.sgnum)
            nbad = count(!r.default_ok for r in rows)
            println(io, "## $(D)D  ($(length(rows)) groups; ",
                    "$nbad insufficient at default `(", join(fill(2, D), ","), ")`)\n")
            hist = Dict{Int,Int}()
            for r in rows
                hist[r.Niso] = get(hist, r.Niso, 0) + 1
            end
            println(io, "Distribution of minimal isotropic N: ",
                    join(("N=$k: $(hist[k])" for k in sort!(collect(keys(hist)))), ", "),
                    " (N=0 ⇒ not found ≤ cap)\n")
            if nbad > 0
                println(io, "| sg | label | Bravais | min. iso. N | min. aniso. | extra ops at default | extra transl. | valid. |")
                println(io, "|---|---|---|---|---|---|---|---|")
                for r in rows
                    r.default_ok && continue
                    println(io, "| $(r.sgnum) | $(r.label) | $(r.bt) | $(r.Niso) | ",
                            join(("(" * join(t, ',') * ")" for t in r.aniso), " "), " | ",
                            join(r.extra_ops_default, ", "), " | ",
                            r.extra_transl_default, " | ", r.valstatus, " |")
                end
                println(io)
            end
        end
        nfail = count(r -> r.valstatus == :failed || !isempty(r.err), allrows)
        println(io, "Validation failures/errors: $nfail\n")
    end
end

"""
Write the `src/tables/fourier_truncations.jl` data file: for each dimension and space-group
number, the antichain of minimal symmetry-safe `idxmax` tuples. Only these antichains are
tabulated; the default (isotropic) truncation is derived from them at run time.
"""
function emit_table(path::AbstractString, allrows)
    # all of 1D, 2D, & 3D must be present: the emitted accessors reference every
    # `SAFE_IDXMAX_ANTICHAINS_$(D)D` constant, so a partial table would not even load
    sort!(unique(r.D for r in allrows)) == [1, 2, 3] ||
        error("cannot emit table: all of `--dims=1,2,3` must be computed")
    open(path, "w") do io
        println(io, """
        # Antichains of minimal symmetry-safe Fourier truncations `idxmax` for
        # `levelsetlattice`, i.e. the minimal elements of the (upward-closed) set of `idxmax`
        # for which a generically modulated Fourier lattice has exactly the symmetry of the
        # associated group — and no more. A given `idxmax` is symmetry-safe if and only if it
        # dominates (componentwise) at least one element of the associated antichain; see
        # `is_symmetry_safe` and `default_idxmax`.
        #
        # Generated (incl. exactness proof sketch, numerical cross-validation, and a
        # completeness check of each antichain over the grid {1..7}ᴰ) via:
        #   julia --project=. build/compute_minimal_fourier_truncations.jl \\
        #       --verify --emit-table=src/tables/fourier_truncations.jl
        """)
        for D in sort!(unique(r.D for r in allrows))
            rows = sort!([r for r in allrows if r.D == D]; by=r -> r.sgnum)
            length(rows) == MAX_SGNUM[D] ||
                error("cannot emit table: $(D)D is incomplete ($(length(rows)) groups)")
            println(io, "const SAFE_IDXMAX_ANTICHAINS_$(D)D = Vector{NTuple{$D,Int}}[")
            for r in rows
                isempty(r.aniso) && error("empty antichain for ($(D)D, sg $(r.sgnum))")
                isempty(r.err) || error("errored group in table emission: sg $(r.sgnum)")
                # NB: a 1-tuple must be emitted as `(n,)`; `(n)` would parse as an `Int`
                entries = join(("(" * join(t, ", ") * (D == 1 ? ",)" : ")")
                                for t in r.aniso), ", ")
                println(io, "    [", entries, "], # ", r.sgnum, ": ", iuc(r.sgnum, D))
            end
            println(io, "]\n")
        end
        println(io, """
        safe_idxmax_antichain(sgnum::Integer, ::Val{1}) = SAFE_IDXMAX_ANTICHAINS_1D[sgnum]
        safe_idxmax_antichain(sgnum::Integer, ::Val{2}) = SAFE_IDXMAX_ANTICHAINS_2D[sgnum]
        safe_idxmax_antichain(sgnum::Integer, ::Val{3}) = SAFE_IDXMAX_ANTICHAINS_3D[sgnum]""")
    end
    @info "wrote table to $path"
end

function main(args=ARGS)
    dims = [1, 2, 3]
    outdir = nothing # opt-in: only write the audit artifacts if `--outdir` is given
    smoke = false
    isocap = 8
    verify = false
    verify_extent = 7
    tablepath = nothing
    for a in args
        if startswith(a, "--dims=")
            dims = parse.(Int, split(chopprefix(a, "--dims="), ','))
        elseif startswith(a, "--outdir=")
            outdir = chopprefix(a, "--outdir=")
        elseif a == "--smoke"
            smoke = true
        elseif startswith(a, "--isocap=")
            isocap = parse(Int, chopprefix(a, "--isocap="))
        elseif a == "--verify"
            verify = true
        elseif startswith(a, "--verify-extent=")
            verify_extent = parse(Int, chopprefix(a, "--verify-extent="))
        elseif startswith(a, "--emit-table=")
            tablepath = chopprefix(a, "--emit-table=")
        else
            error("unrecognized argument: $a")
        end
    end
    outdir === nothing || mkpath(outdir)
    outdir === nothing && tablepath === nothing &&
        @info "no `--outdir` or `--emit-table` given: validating only, writing nothing"

    @info "running self-tests"
    selftests()

    allrows = Any[]
    for D in dims
        sgmax = MAX_SGNUM[D]
        sgnums = smoke && D == 3 ?
            [1, 2, 5, 24, 43, 61, 70, 76, 92, 110, 142, 146, 161, 199, 205, 218, 220, 222, 227, 228, 230] :
            collect(1:sgmax)
        @info "processing $(length(sgnums)) groups in $(D)D"
        foreach(n -> holohedry_rotations(n, Val(D)), sgnums) # warm cache serially
        rows = Vector{Any}(undef, length(sgnums))
        Threads.@threads :dynamic for i in eachindex(sgnums)
            sgnum = sgnums[i]
            rows[i] = try
                process_group(sgnum, Val(D); isocap, validate=(sgnum % 11 == 0),
                              verify, verify_extent)
            catch e
                (; D, sgnum, label=iuc(sgnum, D), bt=bravaistype(sgnum, D),
                   norb_default=-1, default_ok=false, Niso=-1, aniso=NTuple{D,Int}[],
                   extra_ops_default=String[], extra_transl_default="",
                   valstatus=:failed, valmsg="", err=sprint(showerror, e))
            end
        end
        for r in rows
            isempty(r.err) || @warn "error for ($(D)D, sg $(r.sgnum)): $(r.err)"
            r.valstatus == :failed && @warn "validation failed for ($(D)D, sg $(r.sgnum)): $(r.valmsg)"
        end
        outdir === nothing || write_csv(joinpath(outdir, "$(D)d.csv"), rows)
        append!(allrows, rows)
        nbad = count(r -> !r.default_ok, rows)
        @info "$(D)D done: $nbad/$(length(rows)) groups insufficient at default truncation"
    end
    if outdir !== nothing
        write_summary(joinpath(outdir, "summary.md"), allrows)
        @info "wrote audit artifacts to $outdir"
    end

    if tablepath !== nothing
        nfail = count(r -> r.valstatus == :failed || !isempty(r.err), allrows)
        nfail == 0 || error("refusing to emit table: $nfail group(s) failed validation")
        verify || @warn "emitting table without `--verify`: antichain completeness unchecked"
        smoke && error("refusing to emit table from a `--smoke` run")
        emit_table(tablepath, allrows)
    end
    return allrows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
