## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION NOTES / THEORY

#=
# Finding an explicitly real form of irrep matrices

An explicit real, or physically real, form of a set of irrep matrices is one where the
associated matrices $D(g)$ have the following property:

$$
D(g) = D^*(g),
$$

for all operations $g$ in the considered group $G$.

The standard listings of irreps is not explicitly real. However, if an irrep is either
intrinsically real - or has been made into a corep in the complex or pseudoreal case - it is
always equivalent to an intrinsically real form. That is, there exists a unitary transform
$S$ such that:

$$
S D(g) S^{-1} = S D(g) S^\dagger = D^*(g).
$$

Suppose we can find this unitary transformation $S$ by some means. What we are interested
in, is finding a related transform $W$, defining an explicitly real form of the irrep:

$$
\tilde{D}(g) = W D(g) W^{-1} = W D(g) W^\dagger,
$$

where $W$ is some other unitary transformation and where $\tilde{D}(g)$ is an intrinsically
real form of $D(g)$, i.e., where

$$
\tilde{D}(g) = \tilde{D}^*(g) \quad \forall g\in G.
$$

Our aim is to find $W$, assuming we know $S$. For brevity, we will often write $D_g$ in
place of $D(g)$.
First, note that $S$ is not merely a unitary matrix: rather, since, by assumption $D(g)$
is equivalent to a real matrix, what we really mean is that $S$ is also a _symmetric_
unitary matrix, i.e., $S = S^{\mathrm{T}}$ and $S^{-1} = S^\dagger$ (implifying, jointly,
$S^* = S^\dagger = S^{-1}$); this is e.g., derived in Inui p. 74 (bottom) to 75 (top).
Accordingly $S$ is also normal, i.e., $S S^* = S^* S$.

This property in turn implies that we can express $S$ as the square of another symmetric
unitary matrix, say $W$, in the sense that $S = W^2$. This follows from the following
manipulations (Inui, p. 75 bottom), involving the eigendecomposition $S = V Λ V^{-1}$ where
$\Lambda$ is a diagonal matrix with unit-modulus values and $V$ are a set of real
eigenvectors (real because $S$ is symmetric unitary) and
$V^{-1} = V^\dagger = V^{\mathrm{T}}$ (since $S$ is normal).

$$
S = VΛV^{-1} = VΛ^{1/2}Λ^{1/2}V^{\mathrm{T}}
= (VΛ^{1/2}V^{\mathrm{T}})(VΛ^{1/2}V^{\mathrm{T}})
$$,

so we can pick $W = VΛ^{1/2}V^{\mathrm{T}}$ (note also that the square root of $\Lambda$
exist since it is diagonal).
Hence $W^* = V(Λ^{1/2})^*V^{\mathrm{T}} = VΛ^{-1/2}V^{\mathrm{T}} = W^{-1}$ and
$W^{\mathrm{T}} = (VΛ^{1/2}V^{\mathrm{T}})^{\mathrm{T}} =
(V^{\mathrm{T}})^{\mathrm{T}}(Λ^{1/2})^{\mathrm{T}} V^{\mathrm{T}}
= VΛ^{1/2}V^{\mathrm{T}} = W$. I.e., $W$ is also unitary symmetric and normal.

Now, let us rewrite $S D(g) S^{-1} = D^*(g)$ in terms of $W$:

$$
WW D_g W^{-1}W^{-1} = D_g^* \\
$$

Multiply from LHS by $W^*$ and from RHS by $W$:

$$
W^*WW D_g W^{-1}W^{-1} W = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^* D_g^* (W^{-1})^*
$$
where we have used the properties of $W$ to reduce the expressions. 

Identifying $\tilde{D}(g) = W D(g) W^{-1}$ we clearly obtain the desired invariance under
complex conjugation since 
$\tilde{D}^*(g) = (W D_g W^{-1})^* = W^* D_g^* (W^{-1})^* = W D_g W^{-1} = \tilde{D}(g)$.
=#

## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION

"""
    physical_realify(ir::Union{<:PGIrrep, <:SiteIrrep})

Return a manifestly real form of an input irrep `ir` (also called a physically real irrep).

The input irrep must be either a [`PGIrrep`](@ref) or a [`SiteIrrep`](@ref) and must be
equivalent to a real irrep: i.e., the irrep has [`Reality`](@ref) type `REAL` or is a
`PSEUDOREAL` or `COMPLEX` that has already been passed through `realify` and glued together
with its partner (i.e., `iscorep(ir) = true`).

See also [`physical_realify(::Collection{<:Union{<:PGIrrep, <:SiteIrrep}})`](@ref)
for application to a collection of irreps.

## Implementation
A symmetric, unitary transformation is found that maps the irrep matrices to a manifestly
real form. The resulting transformed irrep matrices `Zs` have the property that
`all(Zᵢ -> Zᵢ ≈ conj(Zᵢ), Zs)` is true.

## Examples
```jldoctest
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
  ├ 3₀₀₁⁺: ⎡   -0.5  0.866 ⎤
  │        ⎣ -0.866   -0.5 ⎦
  ├ 3₀₀₁⁻: ⎡  -0.5  -0.866 ⎤
  │        ⎣ 0.866    -0.5 ⎦
  ├  m₁₁₀: ⎡ 0  1 ⎤
  │        ⎣ 1  0 ⎦
  ├  m₁₀₀: ⎡ -0.866   -0.5 ⎤
  │        ⎣   -0.5  0.866 ⎦
  ├  m₀₁₀: ⎡ 0.866    -0.5 ⎤
  └        ⎣  -0.5  -0.866 ⎦
```
"""
function physical_realify(ir::Union{<:PGIrrep, <:SiteIrrep})
    if !ir.iscorep && reality(ir) != REAL
        error("cannot build physically real irreps for input that are complex/pseudoreal \
               and not already converted to a corep")
    end
    Xs = ir.matrices

    if all(X -> real(X) ≈ X, Xs)
        # already in real form: return a copy, to be safe against any subsequent mutation
        return _pgirrep_or_siteirrep_from_matrices(ir, [copy(X) for X in Xs])
    end

    _S = find_symmetric_intertwiner(Xs)
    S = mapto_canonical_unitary(_S)
    if !(S*conj(S) ≈ I)
        error("S is not a symmetric unitary matrix; input cannot be mapped to a real form")
    end
    if !(all(X -> S*X*S' ≈ conj(X), Xs))
        error("S is not a mapping between complex conjugates: input might not be \
               equivalent to a real form")
    end

    λ :: Vector{ComplexF64}, V = eigen(S)
    rV = real(V)
    V ≈ rV || error(lazy"expected real V, got $V")
    V = rV        # eigvecs of symmetry unitary matrix can always be chosen real
    λ ./= abs.(λ) # enforce unitarity explicitly (cf. floating point errors)
    W = V * sqrt(Diagonal(λ)) * V'

    Zs = map(X->W*X*W', Xs)
    all(Z ≈ real(Z) for Z in Zs) || error("obtained real matrices are not real")

    return _pgirrep_or_siteirrep_from_matrices(ir, Zs)
end

function _pgirrep_or_siteirrep_from_matrices(
        ir::PGIrrep{D}, Zs::Vector{Matrix{ComplexF64}}
) where D
    return PGIrrep{D}(ir.cdml, group(ir), Zs, reality(ir), ir.iscorep)
end
function _pgirrep_or_siteirrep_from_matrices(
        ir::SiteIrrep{D}, Zs::Vector{Matrix{ComplexF64}}
) where D
    return SiteIrrep{D}(ir.cdml, group(ir), Zs, reality(ir), ir.iscorep, ir.pglabel)
end

"""
    physical_realify(irs::Collection{T}) where T <: Union{<:PGIrrep, <:SiteIrrep}

Return a manifestly real form of `irs` (also known as physically real irreps),
where `irs` is a [`Collection`](@ref) of either [`PGIrrep`](@ref)s or [`SiteIrrep`](@ref)s.

The input irreps may or may not have already been passed through [`realify`](@ref) (and thus
already glued together with any pseudoreal or complex partners); if they have not, the input
is first passed through `realify`.

See also [`physical_realify(::Union{<:PGIrrep, <:SiteIrrep})`](@ref) for application to
individual irreps.

## Examples
```jldoctest
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
    ├ 3⁺: ⎡   -0.5  0.866 ⎤
    │     ⎣ -0.866   -0.5 ⎦
    ├ 3⁻: ⎡  -0.5  -0.866 ⎤
    │     ⎣ 0.866    -0.5 ⎦
    ├  2: ⎡ 1  0 ⎤
    │     ⎣ 0  1 ⎦
    ├ 6⁻: ⎡   -0.5  0.866 ⎤
    │     ⎣ -0.866   -0.5 ⎦
    ├ 6⁺: ⎡  -0.5  -0.866 ⎤
    └     ⎣ 0.866    -0.5 ⎦

Γ₄Γ₆┌  1: ⎡ 1  0 ⎤
    │     ⎣ 0  1 ⎦
    ├ 3⁺: ⎡   -0.5  0.866 ⎤
    │     ⎣ -0.866   -0.5 ⎦
    ├ 3⁻: ⎡  -0.5  -0.866 ⎤
    │     ⎣ 0.866    -0.5 ⎦
    ├  2: ⎡ -1   0 ⎤
    │     ⎣  0  -1 ⎦
    ├ 6⁻: ⎡   0.5  -0.866 ⎤
    │     ⎣ 0.866     0.5 ⎦
    ├ 6⁺: ⎡    0.5  0.866 ⎤
    └     ⎣ -0.866    0.5 ⎦
```
"""
function physical_realify(irs::Collection{T}) where T<:Union{<:PGIrrep, <:SiteIrrep}
    if any(ir -> ir.iscorep, irs)
        return Collection{T}(map(physical_realify, irs))
    else
        return Collection{T}(map(physical_realify, realify(irs)))
    end
end

## --------------------------------------------------------------------------------------- #
# UTILITIES FOR FINDING UNITARY TRANSFORM BETWEEN IRREP AND ITS COMPLEX CONJUGATE MATRICES

# find an unitary symmetric transform that maps Xᵢ and to Xᵢ* by UXᵢU⁻¹ = Xᵢ* ∀i
# for some set {i}, following https://mathoverflow.net/a/391741. Or more precisely - and
# this is not quite the same - finds the "intertwiner" U such that UXᵢ = Xᵢ*U ∀i.
# We explicitly ensure that `U` is symmetric - but it may not initially be in a form that
# is unitary. Errors if there is more than one possible choice of intertwiner, up a scalar.
# The approach is borrowed and adapted from a similar approach in KdotP.jl.
function find_symmetric_intertwiner(Xs)
    N = LinearAlgebra.checksquare(first(Xs))
    M = length(Xs)
    if !all(X->LinearAlgebra.checksquare(X)==N, Xs)
        error("`Xs` matrices are not square or not of identical size")
    end
    
    N² = N^2
    T = eltype(first(Xs))
    Q = Matrix{T}(undef, M*N², N²)
    for i in 1:M
        Q[(i-1)*N² .+ (1:N²), 1:N²] .= (kron(I(N), transpose(Xs[i]))  .-  
                                        kron(conj.(Xs[i]), I(N)))
    end

    # also enforce that U is symmetric, i.e., that U - Uᵀ = 0
    S = zeros(T, N², N²)
    for i in 1:N²
        S[i, i] += 1  # equiv. to "U[n,m]" of `Uv` vector
                      # ↓ equiv. to transposed `Uv`; effectively "U[m,n]" position in `Uv`
        S[i, Base._sub2ind((N, N), reverse(Base._ind2sub((N, N), i))...)] += -1 # ↑  
    end

    # remove zero-rows from Q and S for performance and numerical stability
    vs = Vector{Vector{T}}()
    for q in eachrow(Q)
        norm(q) < DEFAULT_ATOL || push!(vs, q)
    end
    for s in eachrow(S)
        iszero(s) || push!(vs, s)
    end
    if isempty(vs) # if there were no constraints; e.g., already real or trivial irrep
        return Matrix{T}(I(N))
    end
    C = stack(vs; dims=1)

    # finally, to find a matrix U which obeys UXᵢ = Xᵢ*U ∀i and is also symmetric, we solve
    # for the nullspace of `C` (note that `U` is not necessarily unitary at this point)
    Uv = nullspace(C, atol=1e-10)
    if size(Uv, 2) ≠ 1
        error("failed to determine a unique symmetric intertwiner: nullspace dimension \
               different from 1")
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
