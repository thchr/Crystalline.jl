## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION NOTES / THEORY

#=
# Finding an explicitly real form of irrep matrices (spinless case)

An explicitly real, or physically real, form of a set of irrep matrices is one where the
associated matrices D(g) have the property

    D(g) = D*(g)   ∀ g ∈ G,

for all operations g of the considered group G.

The standard listings of irreps are not explicitly real. However, if an irrep is either
intrinsically real - or has been made into a corep in the complex or pseudoreal case - it is
always equivalent to an intrinsically real form. That is, there exists a unitary transform Γ
such that

    Γ D*(g) Γ⁻¹ = Γ D*(g) Γ† = D(g).

Γ is the unitary part of time reversal in the basis of D (see the spinful section below).
Note that its direction matters: the transform mapping D onto D*, rather than the other way,
is Γ⁻¹ = Γ*, and building W below from it gives a different - though equally real - form.

Suppose we can find Γ by some means. What we are interested in is finding a related
transform W, defining an explicitly real form Dʳ of the irrep,

    Dʳ(g) = W⁻¹ D(g) W = W† D(g) W,

where W is some other unitary transformation and where Dʳ(g) is an intrinsically real form
of D(g), i.e., where

    Dʳ(g) = Dʳ*(g)   ∀ g ∈ G.

Our aim is to find W, assuming we know Γ. For brevity, we often write Dg in place of D(g).
First, note that Γ is not merely a unitary matrix: rather, since D(g) is by assumption
equivalent to a real matrix, Γ is also a *symmetric* unitary matrix, i.e., Γ = Γᵀ and
Γ⁻¹ = Γ† (implying, jointly, Γ* = Γ† = Γ⁻¹); this is derived in e.g. Inui p. 74 (bottom) to
75 (top). Accordingly, Γ is also normal, i.e., Γ Γ* = Γ* Γ.

This property in turn implies that we can express Γ as the square of another symmetric
unitary matrix, say W, in the sense that Γ = W². This follows from the manipulations below
(Inui, p. 75 bottom), involving the eigendecomposition Γ = V Λ V⁻¹, where Λ is a diagonal
matrix with unit-modulus values and V is a set of real eigenvectors (real because Γ is
symmetric unitary), with V⁻¹ = V† = Vᵀ (since Γ is normal):

    Γ = V Λ V⁻¹ = V Λ^½ Λ^½ Vᵀ = (V Λ^½ Vᵀ)(V Λ^½ Vᵀ),

so we can pick W = V Λ^½ Vᵀ (the square root of Λ exists since it is diagonal). Hence
W* = V (Λ^½)* Vᵀ = V Λ^-½ Vᵀ = W⁻¹, and Wᵀ = (V Λ^½ Vᵀ)ᵀ = V Λ^½ Vᵀ = W, i.e., W is also
unitary, symmetric, and normal.

Now, let us rewrite Γ D*(g) Γ⁻¹ = D(g) in terms of W:

    W W Dg* W⁻¹ W⁻¹ = Dg.

Multiplying from the left by W⁻¹ and from the right by W:

    W Dg* W⁻¹ = W⁻¹ Dg W.

Identifying Dʳ(g) = W⁻¹ D(g) W, we obtain the desired invariance under complex conjugation,
since

    Dʳ*(g) = (W⁻¹ Dg W)* = (W⁻¹)* Dg* W* = W Dg* W⁻¹ = W⁻¹ Dg W = Dʳ(g),

where we used W* = W⁻¹ (and hence (W⁻¹)* = W) in the third step.

Equivalently, and in the language of the spinful case below: the basis change D → W†DW sends
Γ → W†ΓW* = WW* = 𝟙, i.e., it brings the unitary part of time reversal to its canonical
spinless value, just as the spinful procedure brings it to J.

# Spinful (double-valued) irreps

For spinful irreps, time reversal squares to -1 and it is not always possible to make a
basis choice that makes the matrices real. The canonical form is instead fixed through the
unitary part Γ of time reversal T = ΓK, with K denoting complex conjugation, by requiring

    Γ D*(g) Γ† = D(g)   ∀ g ∈ G,

with Γ = J ≡ iσʸ ⊗ 𝟙ₙ for an irrep of dimension 2n. The spinless case is the same statement
with Γ = 𝟙, which reduces to D(g) = D*(g); `timereversal_unitary` returns the applicable Γ.
The outer index of J is the Kramers index, matching how `realify` glues an irrep to its
partner. Under a basis change ψ → Wψ, the matrices transform as D → W D W†, but Γ transforms
as Γ → W Γ Wᵀ - with a transpose, rather than an adjoint, since K W† = Wᵀ K. A global phase
W = exp(iα)𝟙 therefore leaves D alone while sending Γ → exp(2iα)Γ, so the phase of Γ is free
and can always be fixed to make Γ = J.

Which Γ occurs follows from Schur's lemma: as above, ΓΓ* = c𝟙 with c = ±1, and c = +1
exactly when Γ is symmetric (the real case), c = -1 exactly when it is antisymmetric (the
pseudoreal case). An invertible antisymmetric matrix has even dimension, so pseudoreal
irreps are even-dimensional; there, T² = ΓΓ* = -1 is Kramers' theorem.

Accordingly, there are three cases, each reached differently:

- REAL (glued by `realify` into D₀ ⊕ D₀): the *undoubled* D₀ is equivalent to a real irrep
  and is made explicitly real by the spinless procedure above; the doubled diag(D₀, D₀) then
  commutes with J, as required. Note that the procedure cannot be applied to the doubled
  matrices directly: their symmetric intertwiner is not unique (the doubling introduces
  additional intertwiners).
- COMPLEX (glued into D₀ ⊕ D₀′ with D₀′ ≃ D₀*): mapping D₀′ onto D₀* by an intertwiner and
  then transforming by S = 2^-½ [𝟙 𝟙; -i𝟙 i𝟙] takes diag(D₀, D₀*) to [A -B; B A], writing
  D₀ = A + iB with real A and B. This is real, and matrices of this form commute with J.
- PSEUDOREAL (left alone by `realify`): the irrep has no real form, and its unique
  intertwiner Γ is antisymmetric. A basis of Kramers pairs {vₖ, -Tvₖ} brings Γ to J, since
  v ⟂ Tv and the orthogonal complement of a pair is T-invariant. The matrices stay complex.

In every case, Γ D*(g) Γ† = D(g) with Γ = J is equivalent to the block form

    D(g) = [A B; -B* A*],

the n-block analogue of the SU(2) form [a b; -b* a*].
=#

## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION

"""
    timereversal_unitary(ir) --> AbstractMatrix{<:Real}

Return the unitary part `Γ` of time reversal `T = ΓK`, with `K` denoting complex
conjugation, in the basis established by [`physical_realify`](@ref), for a point group or
site symmetry irrep `ir`.

For spinless irreps ([`PGIrrep`](@ref), [`SiteIrrep`](@ref)), `Γ` is the identity matrix
(returned as `I(n)`), so that time reversal is plain complex conjugation. For spinful irreps
([`DPGIrrep`](@ref), [`DSiteIrrep`](@ref)) of dimension `2n`, `Γ` is
`J = iσʸ ⊗ 𝟙ₙ = [0 𝟙ₙ; -𝟙ₙ 0]`, so that `T² = -𝟙`; the outer index of this block structure
is the Kramers index, matching the way [`realify`](@ref) glues an irrep together with its
time-reversal partner. `Γ` is real in either case, and is returned as a real-valued matrix.

The matrices `Ds = physical_realify(ir).matrices` obey `Γ*conj(D)*Γ' ≈ D` for every `D` in
`Ds`, which for spinless irreps is just `D ≈ conj(D)`.

`Γ` is a choice of convention: its phase is free, and, for spinful irreps, so is a
subsequent symplectic rotation. Any use of the irrep matrices alongside `Γ` — e.g., in
time-reversal symmetric tight-binding models — must adopt the same convention.
"""
function timereversal_unitary(ir::Union{AbstractPGIrrep, AbstractSiteIrrep})
    N = irdim(ir)
    if !isspinful(ir)
        return I(N) # 𝟙
    else
        if isodd(N)
            error("a spinful irrep of odd dimension ($N) has no associated unitary: it \
                   must first be glued to its time-reversal partner by `realify`")
        end
        # `J = iσʸ ⊗ 𝟙ₙ = [0 𝟙ₙ; -𝟙ₙ 0]`; equivalent to `kron([0 1; -1 0], I(n))`, but
        # written out, since the Kronecker product only places these ±1s
        n = N ÷ 2
        J = zeros(Float64, N, N)
        @inbounds for i in OneTo(n)
            J[i, n+i] = 1
            J[n+i, i] = -1
        end
        return J
    end
end

"""
    physical_realify(ir::Union{AbstractPGIrrep, AbstractSiteIrrep})

Return a canonical form of an input irrep `ir` under time reversal.

For spinless irreps ([`PGIrrep`](@ref), [`SiteIrrep`](@ref)), this is a manifestly real
form, also called a physically real irrep. For spinful irreps ([`DPGIrrep`](@ref),
[`DSiteIrrep`](@ref)), time reversal squares to `-1` and a real form does not generally
exist; the returned matrices instead obey `J*conj(D)*J' ≈ D` with `J = iσʸ ⊗ 𝟙ₙ`, and are
additionally real whenever that is possible (i.e., unless `ir` is pseudoreal). See
[`timereversal_unitary`](@ref), which returns the relevant matrix in either case. This
canonical form is unique only up to conjugation by a real orthogonal matrix for spinless
irreps, or by a unitary symplectic matrix for spinful irreps.

The input irrep must be one that [`realify`](@ref) leaves unchanged: i.e., either an irrep
that has already been glued together with its time-reversal partner (`iscorep(ir) = true`)
or one whose [`Reality`](@ref) type requires no gluing (`REAL` for spinless irreps,
`PSEUDOREAL` for spinful ones).

See also [`physical_realify(::Collection{<:Union{<:PGIrrep, <:SiteIrrep}})`](@ref)
for application to a collection of irreps.

## Examples
```jldoctest; filter = [r"-" => "", r" +" => " "]
julia> pgir = pgirreps(19, 3)[end]
Γ₃┌     1: ⎡ 1  0 ⎤
  │        ⎣ 0  1 ⎦
  ├ 3₀₀₁⁺: ⎡ -0.5+0.866im             0 ⎤
  │        ⎣            0  -0.5-0.866im ⎦
  ├ 3₀₀₁⁻: ⎡ -0.5-0.866im             0 ⎤
  │        ⎣            0  -0.5+0.866im ⎦
  ├  m₁₁₀: ⎡ 0  1 ⎤
  │        ⎣ 1  0 ⎦
  ├  m₁₀₀: ⎡            0  -0.5-0.866im ⎤
  │        ⎣ -0.5+0.866im             0 ⎦
  ├  m₀₁₀: ⎡            0  -0.5+0.866im ⎤
  └        ⎣ -0.5-0.866im             0 ⎦

julia> reality(pgir)
REAL::Reality = 1

julia> physical_realify(pgir)
Γ₃┌     1: ⎡ 1  0 ⎤
  │        ⎣ 0  1 ⎦
  ├ 3₀₀₁⁺: ⎡  -0.5  -0.866 ⎤
  │        ⎣ 0.866    -0.5 ⎦
  ├ 3₀₀₁⁻: ⎡   -0.5  0.866 ⎤
  │        ⎣ -0.866   -0.5 ⎦
  ├  m₁₁₀: ⎡ 0  1 ⎤
  │        ⎣ 1  0 ⎦
  ├  m₁₀₀: ⎡ 0.866    -0.5 ⎤
  │        ⎣  -0.5  -0.866 ⎦
  ├  m₀₁₀: ⎡ -0.866   -0.5 ⎤
  └        ⎣   -0.5  0.866 ⎦
```
"""
function physical_realify(ir::Union{AbstractPGIrrep, AbstractSiteIrrep})
    spinful = isspinful(ir)
    if !iscorep(ir)
        # the input must be an irrep that `realify` leaves alone; that is the REAL ones for
        # spinless irreps and the PSEUDOREAL ones for spinful irreps (`realify` doubles or
        # glues to a partner in every other case)
        r_untouched = spinful ? PSEUDOREAL : REAL
        if reality(ir) ≠ r_untouched
            error("cannot build a canonical time-reversal form of a $(reality(ir)) irrep \
                   that is not already a corep: pass the input through `realify` first")
        end
    end

    # a real form exists in every case but the spinful pseudoreal one
    hasrealform = !spinful || reality(ir) ≠ PSEUDOREAL

    Xs = ir.matrices
    if _iscanonical(Xs, spinful, hasrealform)
        # `ir` is already in the form that the transformations below would produce, so there
        # is nothing to do; return a copy, to be safe against any subsequent mutation
        return _irrep_from_matrices(ir, [copy(X) for X in Xs])
    end

    Zs = spinful ? _spinful_canonical_matrices(ir) : _real_matrices(Xs)
    _iscanonical(Zs, spinful, hasrealform) ||
        error("obtained matrices are not in canonical time-reversal form")

    return _irrep_from_matrices(ir, Zs)
end

# whether `Zs` is already in the canonical form: explicitly real, if a real form exists, and
# obeying `Γ*conj(Z)*Γ' = Z` for the `Γ` of `timereversal_unitary`. For spinless irreps the
# latter is implied by the former (`Γ = 𝟙`); for spinful ones it is the block form checked
# by `_iskramersblocked`. Both checks avoid forming matrix products, since they run on every
# call, whether or not there is anything to do
function _iscanonical(Zs, spinful::Bool, hasrealform::Bool)
    hasrealform && !all(_isapproxreal, Zs) && return false
    return !spinful || all(_iskramersblocked, Zs)
end

_isapproxreal(Z) = all(z -> abs(imag(z)) ≤ DEFAULT_ATOL, Z)

# whether `Z` has the block form `[A B; -B* A*]`, i.e. whether `J*conj(Z)*J' = Z` for
# `J = iσʸ ⊗ 𝟙ₙ`; written out blockwise, since `J`'s products merely permute and negate
function _iskramersblocked(Z)
    n = LinearAlgebra.checksquare(Z) ÷ 2
    for j in OneTo(n), i in OneTo(n)
        abs(Z[i+n, j+n] - conj(Z[i, j]))   ≤ DEFAULT_ATOL || return false # A* block
        abs(Z[i+n, j]   + conj(Z[i, j+n])) ≤ DEFAULT_ATOL || return false # -B* block
    end
    return true
end

function _irrep_from_matrices(ir::AbstractPGIrrep, Zs::Vector{Matrix{ComplexF64}})
    return typeof(ir)(ir.cdml, group(ir), Zs, reality(ir), ir.iscorep)
end
function _irrep_from_matrices(ir::AbstractSiteIrrep, Zs::Vector{Matrix{ComplexF64}})
    return typeof(ir)(ir.cdml, group(ir), Zs, reality(ir), ir.iscorep, ir.pglabel)
end

"""
    physical_realify(irs::Collection{T})
            where T <: Union{AbstractPGIrrep, AbstractSiteIrrep}

Return a canonical form of `irs` under time reversal (see
[`physical_realify(::Union{<:PGIrrep, <:SiteIrrep})`](@ref) for the individual irrep case
and for the meaning of "canonical"), where `irs` is a [`Collection`](@ref) of point group or
site symmetry irreps.

The input irreps may or may not have already been passed through [`realify`](@ref) (and thus
already glued together with any pseudoreal or complex partners); if they have not, the input
is first passed through `realify`.

## Examples
```jldoctest; filter = [r"-" => "", r" +" => " "]
julia> pgirs = pgirreps(9,2);

julia> physical_realify(pgirs)
4-element Collection{PGIrrep{2}} for ⋕9 (6):
Γ₁┌  1: 1
  ├ 3⁺: 1
  ├ 3⁻: 1
  ├  2: 1
  ├ 6⁻: 1
  └ 6⁺: 1

Γ₂┌  1: 1
  ├ 3⁺: 1
  ├ 3⁻: 1
  ├  2: -1
  ├ 6⁻: -1
  └ 6⁺: -1

Γ₃Γ₅┌  1: ⎡ 1  0 ⎤
    │     ⎣ 0  1 ⎦
    ├ 3⁺: ⎡  -0.5  -0.866 ⎤
    │     ⎣ 0.866    -0.5 ⎦
    ├ 3⁻: ⎡   -0.5  0.866 ⎤
    │     ⎣ -0.866   -0.5 ⎦
    ├  2: ⎡ 1  0 ⎤
    │     ⎣ 0  1 ⎦
    ├ 6⁻: ⎡  -0.5  -0.866 ⎤
    │     ⎣ 0.866    -0.5 ⎦
    ├ 6⁺: ⎡   -0.5  0.866 ⎤
    └     ⎣ -0.866   -0.5 ⎦

Γ₄Γ₆┌  1: ⎡ 1  0 ⎤
    │     ⎣ 0  1 ⎦
    ├ 3⁺: ⎡  -0.5  -0.866 ⎤
    │     ⎣ 0.866    -0.5 ⎦
    ├ 3⁻: ⎡   -0.5  0.866 ⎤
    │     ⎣ -0.866   -0.5 ⎦
    ├  2: ⎡ -1   0 ⎤
    │     ⎣  0  -1 ⎦
    ├ 6⁻: ⎡    0.5  0.866 ⎤
    │     ⎣ -0.866    0.5 ⎦
    ├ 6⁺: ⎡   0.5  -0.866 ⎤
    └     ⎣ 0.866     0.5 ⎦
```
"""
function physical_realify(
            irs::Collection{T}
) where T<:Union{AbstractPGIrrep, AbstractSiteIrrep}
    if any(ir -> ir.iscorep, irs)
        return Collection{T}(map(physical_realify, irs))
    else
        return Collection{T}(map(physical_realify, realify(irs)))
    end
end

## --------------------------------------------------------------------------------------- #
# UTILITIES FOR FINDING THE UNITARY TRANSFORM BETWEEN TWO SETS OF IRREP MATRICES

# Find the "intertwiner" `U` that maps the matrices `froms` onto the matrices `tos`, i.e.,
# the `U` with `U*froms[i] = tos[i]*U ∀i` - or, equivalently for invertible `U`, with
# `U*froms[i]*U⁻¹ = tos[i] ∀i`. The direction matters, and the arguments are named for it:
# the `U` mapping the other way is `U⁻¹ = U†`, differing from `U` unless `U` is symmetric.
# `transpose_symmetry` optionally constrains `U` under transposition: `+1` requires
# `U = Uᵀ`, `-1` requires `U = -Uᵀ`, and `0` (the default) imposes no constraint.
# The returned `U` need not be unitary; `mapto_canonical_unitary` maps it to a related
# unitary matrix. Errors if there is more than one choice of intertwiner, up to a scalar.
# Follows https://mathoverflow.net/a/391741; approach borrowed & adapted from KdotP.jl.
function find_intertwiner(froms, tos; transpose_symmetry::Int=0)
    N = LinearAlgebra.checksquare(first(froms))
    M = length(froms)
    if !all(X->LinearAlgebra.checksquare(X)==N, froms) ||
       !all(Y->LinearAlgebra.checksquare(Y)==N, tos) || length(tos) ≠ M
        error("`froms` and `tos` matrices are not square, or not of equal size & number")
    end

    N² = N^2
    T = eltype(first(froms))
    Q = Matrix{T}(undef, M*N², N²)
    for i in 1:M
        Q[(i-1)*N² .+ (1:N²), 1:N²] .= (kron(I(N), transpose(froms[i]))  .-  
                                        kron(tos[i], I(N)))
    end

    # remove zero-rows from Q for performance and numerical stability
    vs = Vector{Vector{T}}()
    for q in eachrow(Q)
        norm(q) < DEFAULT_ATOL || push!(vs, q)
    end

    # optionally also enforce that U is symmetric or antisymmetric, i.e., that U ∓ Uᵀ = 0
    if !iszero(transpose_symmetry)
        S = zeros(T, N², N²)
        for i in 1:N²
            j = Base._sub2ind((N, N), reverse(Base._ind2sub((N, N), i))...)
            S[i, i] += 1                   # equiv. to "U[n,m]" of the `Uv` vector
            S[i, j] -= transpose_symmetry  # equiv. to the transposed "U[m,n]" position
        end
        for s in eachrow(S)
            iszero(s) || push!(vs, s)
        end
    end
    if isempty(vs) # if there were no constraints; e.g., already real or trivial irrep
        return Matrix{T}(I(N))
    end
    C = stack(vs; dims=1)

    # finally, to find a matrix U which obeys the above constraints, we solve for the
    # nullspace of `C` (note that `U` is not necessarily unitary at this point)
    Uv = nullspace(C, atol=1e-10)
    if size(Uv, 2) ≠ 1
        error("failed to determine a unique intertwiner: nullspace dimension different \
               from 1")
    end
    return permutedims(reshape(Uv[:,1], N, N))
end

function mapto_canonical_unitary(U)
    # we assume that `U` is at most "a scalar away" from being unitary, and we want to
    # mutate `U` to become this related unitary matrix; additionally, we want to map it to
    # a canonical form. We have two steps: 

    # 1. "unitarize" `U` (a unitary operator must have UU† = I; let's ensure this - assuming
    # that `U` is indeed a "scalar away" from being unitary). Trick is from Wigner p. 78-79:
    K = inv(sqrt(U*U')) # a multiple of the identity matrix
    U = K*U

    # 2. "canonicalize" `U` (currently has an arbitrary complex phase of norm 1) by
    # requiring that it is a special unitary matrix (i.e., in SU(N)), so that `det(U) = 1`.
    # Already, `U` is in U(N), i.e. `norm(det(U)) = 1`; below, we find the phase needed to
    # "rotate" it into SU(N), using that det(cA) = cᴺdet(A)
    N = LinearAlgebra.checksquare(U)
    c = det(U)^(1/N)
    U ./= c

    return U
end

## --------------------------------------------------------------------------------------- #
# CANONICAL FORMS (see the theory notes above)

# an explicitly real form of `Xs`, which must be equivalent to a set of real matrices.
# The early return tests realness with the same tolerance as `_iscanonical`, so a caller
# that found `Xs` non-canonical cannot be handed `Xs` straight back
function _real_matrices(Xs)
    all(_isapproxreal, Xs) && return [copy(X) for X in Xs]

    # the unitary part of time reversal, in the same direction as for spinful irreps; here
    # it is symmetric, which is what makes a real form possible at all (see the notes above)
    Γ = mapto_canonical_unitary(find_intertwiner(conj.(Xs), Xs; transpose_symmetry=+1))
    if !(Γ*conj(Γ) ≈ I)
        error("Γ is not a symmetric unitary matrix; input cannot be mapped to a real form")
    end
    if !(all(X -> Γ*conj(X)*Γ' ≈ X, Xs))
        error("Γ does not map the conjugated matrices onto the originals: input might not \
               be equivalent to a real form")
    end

    # with Γ = W² for a symmetric unitary W, the basis change `D -> W'*D*W` sends Γ to 𝟙
    λ :: Vector{ComplexF64}, V = eigen(Γ)
    rV = real(V)
    V ≈ rV || error(lazy"expected real V, got $V")
    V = rV        # eigvecs of a symmetric unitary matrix can always be chosen real
    λ ./= abs.(λ) # enforce unitarity explicitly (cf. floating point errors)
    W = V * sqrt(Diagonal(λ)) * V'

    Zs = map(X->W'*X*W, Xs)
    all(_isapproxreal, Zs) || error("obtained real matrices are not real")

    return Zs
end

# the canonical time-reversal form of a spinful irrep, i.e., matrices obeying
# `J*conj(D)*J' = D`; real, except in the pseudoreal case (see the notes above)
function _spinful_canonical_matrices(ir::Union{AbstractPGIrrep, AbstractSiteIrrep})
    Xs = ir.matrices
    N = irdim(ir)
    n = N ÷ 2

    if !iscorep(ir)
        # PSEUDOREAL (and spinful): `realify` leaves these unglued, and there is no real
        # form; bring `Γ` to `J` in a basis of Kramers pairs, leaving the matrices complex
        Γ = mapto_canonical_unitary(find_intertwiner(conj.(Xs), Xs; transpose_symmetry=-1))
        Γ ≈ -transpose(Γ) || error("obtained intertwiner is not antisymmetric")
        V = _kramers_basis(Γ)
        return map(X -> V'*X*V, Xs)
    end

    D₀ = [X[1:n, 1:n] for X in Xs] # the irrep that `realify` doubled or glued to a partner
    if reality(ir) == REAL
        # REAL, i.e. `D₀ ⊕ D₀`: realify `D₀` alone, then double it again
        all(((X, A),) -> X ≈ _blockdiag2x2(A), zip(Xs, D₀)) ||
            error("real corep is not in the block form `diag(D₀, D₀)` built by `realify`")
        return [_blockdiag2x2(Z) for Z in _real_matrices(D₀)]
    else
        # COMPLEX, i.e. `D₀ ⊕ D₀′` with `D₀′ ≅ D₀*`: map `D₀′` onto `D₀*`, then rotate to a
        # real form
        D₀′ = [X[n+1:N, n+1:N] for X in Xs]
        V = mapto_canonical_unitary(find_intertwiner(D₀′, conj.(D₀)))
        S = kron(ComplexF64[1 1; -im im], I(n)) ./ sqrt(2)
        W = S * _blockdiag2x2(Matrix{ComplexF64}(I, n, n), V)
        return map(X -> W*X*W', Xs)
    end
end

# a unitary `V` whose columns are the Kramers pairs `(v₁, …, vₙ, -Tv₁, …, -Tvₙ)` of the
# antiunitary `T = ΓK`, so that `V'*Γ*conj(V) = J`. Assumes `Γ` unitary and antisymmetric,
# whence `T² = ΓΓ* = -𝟙` and `v ⟂ Tv`; the orthogonal complement of a Kramers pair is itself
# `T`-invariant, so the pairs can be built up one at a time
function _kramers_basis(Γ)
    N = LinearAlgebra.checksquare(Γ)
    n = N ÷ 2
    V = Matrix{ComplexF64}(undef, N, N)
    pairs = Matrix{ComplexF64}(undef, N, 0) # the pairs found so far, as orthonormal columns
    for k in OneTo(n)
        # start the next pair from any unit vector orthogonal to the pairs already found;
        # `nullspace` returns an orthonormal basis of exactly that complement, so take its
        # first column
        v = nullspace(pairs')[:, 1]
        # `w = -Tv` completes the pair: it is a unit vector, since `Γ` is unitary, and is
        # orthogonal to `v`, since `T² = -𝟙`. Together with `Tw = v`, this sign is what
        # makes `V'*Γ*conj(V)` come out as `J` rather than as `-J`
        w = -Γ * conj(v)
        # the `v`s fill the first half of `V` and the `w`s the second, so that the Kramers
        # index ends up as the outer index, as `J = iσʸ ⊗ 𝟙ₙ` requires
        V[:, k]   = v
        V[:, k+n] = w
        pairs = hcat(pairs, v, w)
    end
    return V
end
