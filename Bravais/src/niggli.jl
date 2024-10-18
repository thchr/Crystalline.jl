"""
    nigglibasis(Rs; rtol=1e-5, max_iterations=200)

Given a set of primitive basis vectors `Rs`, return a basis `Rs′` for the corresponding
Niggli-reduced unit cell, as well as a transformation matrix `P`, such that
`Rs′ = transform(Rs, P)` (see [`transform`](@ref)).

## Definition

A Niggli-reduced basis ``(\\mathbf{a}, \\mathbf{b}, \\mathbf{c})`` represents a unique
choice of basis for any given lattice (or, more precisely, a unique choice of the basis
vector lengths ``|\\mathbf{a}|, |\\mathbf{b}|, |\\mathbf{c}|``, and mutual angles between 
``\\mathbf{a}, \\mathbf{b}, \\mathbf{c}``). This uniqueness is one of the main motivations
for computing the Niggli reduction procedure, as it enables easy comparison of lattices.
Additionally, the associated Niggli-reduced basis vectors ``(\\mathbf{a}, \\mathbf{b},
\\mathbf{c})``, fulfil several conditions [^3]:

1. **"Main" conditions:**
    - The basis vectors are sorted by increasing length:
      ``|\\mathbf{a}| ≤ |\\mathbf{b}| ≤ |\\mathbf{c}|``.
    - The angles between basis vectors are either all acute or all non-acute.
2. **"Special" conditions:**
    - Several special conditions, applied in "special" cases, such as
      ``|\\mathbf{a}| = |\\mathbf{b}|`` or 
      `\\mathbf{b}\\cdot\\mathbf{c} = \\tfrac{1}{2}|\\mathbf{b}|^2`. See Ref. [^3] for
      details.

Equivalently, the Niggli-reduced basis fulfils the following geometric conditions (Section
9.3.1 of Ref. [^3]):
- The basis vectors are sorted by increasing length.
- The basis vectors have least possible total length, i.e., ``|\\mathbf{a}| + |\\mathbf{b}|
  + |\\mathbf{c}|`` is minimum. I.e., the associated Niggli cell is a Buerger cell.
- The associated Buerger cell has maximum deviation among all other Buerger cells, i.e., the
  basis vector angles ``α, β, γ`` maximize ``|90° - α| + |90° - β| + |90° - γ|``.
  
## Keyword arguments

- `rtol :: Real`: relative tolerance used in the Grosse-Kunstleve approach for floating point
  comparisons (default: `1e-5`).
- `max_iterations :: Int`: maximum number of iterations in which to cycle the Krivy-Gruber
  steps (default: `200`).

## Implementation

Implementation follows the algorithm originally described by Krivy & Gruber [^1], with the
stability modificiations proposed by Grosse-Kunstleve et al. [^2] (without which the 
algorithm proposed in [^1] simply does not work on floating point hardware).

[^1] I. Krivy & B. Gruber. A unified algorithm for determinign the reduced (Niggli) cell,
    [Acta Crystallogr. A **32**, 297 (1976)](https://doi.org/10.1107/S0567739476000636).
[^2] R.W. Grosse-Kunstleve, N.K. Sauter, & P.D. Adams, Numerically stable algorithms for the
    computation of reduced unit cells,
    [Acta Crystallogr. A **60**, 1 (2004)](https://doi.org/10.1107/S010876730302186X)
[^3] Sections 9.2 & 9.3, International Tables of Crystallography, Volume A, 5th ed. (2005).
"""
function nigglibasis(
            Rs :: DirectBasis{3};
            rtol :: Real = 1e-5, # default relative tolereance, following [2]
            max_iterations :: Int = 200
    )

    # check input
    if length(Rs) ≠ 3
        error("the Niggli reduction implementation only supports 3D settings: at least 3 \
               basis vectors must be given")
    end
    if any(R -> length(R) ≠ 3, Rs)
        error("the Niggli reduction implementation only supports 3D settings: a basis \
               vector was supplied that does not have 3 components")
    end
    rtol < 0 && throw(DomainError(rtol, "relative tolerance `rtol` must be non-negative"))
    max_iterations ≤ 0 && throw(DomainError(max_iterations, "`max_iterations` must be positive"))
    
    # tolerance
    D = 3 # algorithm assumes 3D setting (TODO: extend to 2D?)
    ϵ = rtol * abs(volume(Rs))^(1/D)

    # NB: it would make sense in principle to start by doing a conventional lattice
    #     reduction step here, e.g., `lll` from LLLplus.jl, to make sure we have a good
    #     "starting point" for the algorithm; this would mean we'd need to do fewer of the
    #     iterations below, especially for very skewed input. But LLLplus.jl is not very
    #     fast - usually slower than `nigglibasis` itself, so don't do it for now.
	# Rs_lll = lll(stack(Rs))[1] # LLL-reduced basis
	# Rs = DirectBasis(eachcol(Rs_lll))
    
    # initialization
    A, B, C, ξ, η, ζ = niggli_parameters(Rs)
    P = @SMatrix [1 0 0; 0 1 0; 0 0 1]

    # performing steps A1-A8 until no condition is met
    iteration = 0
    while iteration < max_iterations
        iteration += 1

        # At each step below, we update the transformation matrix `P` by post-multiplying it
        # with `P′` in the sense P ← P * P′. We also update the associated Niggli parameters
        # (A, B, C, ξ, η, ζ) by applying the transformation `P` to the original basis 
        # vectors `Rs`. If we reach the end of the set of steps, without being returned to 
        # step A1, we are done. Executing steps A2, A5, A6, A7, & A8 subsequently returns us
        # to step A1.

        # step A1                                      A > B || (A == B && abs(ξ) > abs(η))
        if A > B + ϵ || (abs(A-B) < ϵ && abs(ξ) > abs(η) + ϵ)
            P′ = @SMatrix [0 -1 0; -1 0 0; 0 0 -1] # swap (A,ξ) ↔ (B,η)
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
        end

        # step A2                                      B > C || (B == C && abs(η) > abs(ζ))
        if B > C + ϵ || (abs(B-C) < ϵ && abs(η) > abs(ζ) + ϵ)
            P′ = @SMatrix [-1 0 0; 0 0 -1; 0 -1 0] # swap (B,η) ↔ (C,ζ)
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
            continue # restart from step A1
        end

        # step A3                                                                 ξ*η*ζ > 0
        if tolsign(ξ, ϵ) * tolsign(η, ϵ) * tolsign(ζ, ϵ) == 1
            i, j, k = ifelse(ξ < -ϵ, -1, 1), ifelse(η < -ϵ, -1, 1), ifelse(ζ < -ϵ, -1, 1)
            P′ = @SMatrix [i 0 0; 0 j 0; 0 0 k] # update (ξ,η,ζ) to (|ξ|,|η|,|ζ|) 
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
        end

        # step A4                                                                 ξ*η*ζ ≤ 0
        l, m, n = tolsign(ξ, ϵ), tolsign(η, ϵ), tolsign(ζ, ϵ)
        if !(l == m == n == -1) && (l*m*n == -1 || l*m*n == 0)
            i, j, k = _stepA4_ijk(l, m, n)
            P′ = @SMatrix [i 0 0; 0 j 0; 0 0 k] # update (ξ,η,ζ) to (-|ξ|,-|η|,-|ζ|)
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
        end

        # step A5                      abs(ξ) > B || (ξ == B && 2η < ζ) || (ξ == -B, ζ < 0)
        if abs(ξ) > B + ϵ || (abs(B-ξ) < ϵ && 2η < ζ - ϵ) || (abs(ξ + B) < ϵ && ζ < -ϵ)
            P′ = @SMatrix [1 0 0; 0 1 ifelse(ξ > 0, -1, 1); 0 0 1]
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
            continue # restart from step A1
        end

        # step A6                    abs(η) > A || (η == A && 2ξ < ζ) || (η == -A && ζ < 0)
        if abs(η) > A + ϵ || (abs(η-A) < ϵ && 2ξ < ζ - ϵ) || (abs(η+A) < η && ζ < -ϵ)
            P′ = @SMatrix [1 0 ifelse(η > 0, -1, 1); 0 1 0; 0 0 1]
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
            continue # restart from step A1
        end

        # step A7                    abs(ζ) > A || (ζ == A && 2ξ < η) || (ζ == -A && η < 0)
        if abs(ζ) > A + ϵ || (abs(ζ-A) < ϵ && 2ξ < η - ϵ) || (abs(ζ+A) < ϵ && η < -ϵ)
            P′ = @SMatrix [1 ifelse(ζ > 0, -1, 1) 0; 0 1 0; 0 0 1]
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
            continue # restart from step A1
        end

        # step A8     ξ + η + ζ + A + B < 0 || (ξ + η + ζ + A + B == 0 && 2(A + η) + ζ > 0)
        if ξ + η + ζ + A + B < -ϵ || (abs(ξ + η + ζ + A + B) < ϵ && 2(A + η) + ζ > ϵ)
            P′ = @SMatrix [1 0 1; 0 1 1; 0 0 1]
            P *= P′
            A, B, C, ξ, η, ζ = niggli_parameters(transform(Rs, P))
            continue # restart from step A1
        end

        # end of steps, without being returned to step A1 → we are done
        Rs′ = transform(Rs, P)
        return Rs′, P
    end

    error("Niggli reduction did not converge after $max_iterations iterations")
end

function niggli_parameters(Rs)
    # The Gram matrix G (or metric tensor) is related to the Niggli parameters (A, B, C, ξ,
    # η, ζ) via
    #       G = [A ζ/2 η/2; ζ/2 B ξ/2; η ξ/2 C/2]
    # and `G = stack(Rs)'*stack(Rs)` = VᵀV.
    a, b, c = Rs[1], Rs[2], Rs[3]
    A = dot(a, a)
    B = dot(b, b)
    C = dot(c, c)
    ξ = 2dot(b, c)
    η = 2dot(c, a)
    ζ = 2dot(a, b)
    return A, B, C, ξ, η, ζ
end

intsign(x::Real) = signbit(x) ? -1 : 1
tolsign(x::Real, ϵ::Real) = abs(x) > ϵ ? intsign(x) : 0 # -1 if x<-ϵ, +1 if x>ϵ, else 0

function _stepA4_ijk(l::Int, m::Int, n::Int)
    i = j = k = 1
    r = 0 # reference to variables i, j, k: r = 0 → undef, r = 1 → i, r = 2 → j, r = 3 → k
    if l == 1 # ξ > ϵ
        i = -1
    elseif l == 0 # ξ ≈ 0
        r = 1 # → `i`
    end
    if m == 1 # η > ϵ
        j = -1
    elseif m == 0 # η ≈ 0
        r = 2 # → `j`
    end
    if n == 1 # ζ > ϵ
        k = -1
    elseif n == 0 # ζ ≈ 0
        r = 3 # → `k`
    end
    if i*j*k == -1
        if r == 1
            i = -1
        elseif r == 2
            j = -1
        elseif r == 3
            k = -1
        else # r == 0 (unset, but needed)
            error("unhandled error in step A4 of Niggli reduction; \
                   failed to set reference `r` to variables `i, j, k`, but \
                   expected reference to be set")
        end
    end
    return i, j, k
end

# ---------------------------------------------------------------------------------------- #
# Implementations in 1D & 2D

niggle_reduce(Rs :: DirectBasis{1}; kws...) = Rs # 1D (trivial)
function nigglibasis( # 2D by 3D-piggybacking
            Rs :: DirectBasis{2}; 
            rtol :: Real = 1e-5, kws...) 
    max_norm = maximum(norm, Rs)
    # create trivially 3D-extended basis by appending a basis vector along z (& make this
    # interim vector longer than any other basis vector, so it is necessarily last in the
    # Niggli basis so it can be readily removed)
    Rs³ᴰ = DirectBasis{3}(SVector(Rs[1]..., 0.0),
                          SVector(Rs[2]..., 0.0), 
                          SVector(0.0, 0.0, 2max_norm)) 
    Rs³ᴰ′, P³ᴰ = nigglibasis(Rs³ᴰ; rtol=rtol, kws...)
    # extract associated 2D basis and transformation
    Rs′ = DirectBasis{2}(Rs³ᴰ′[1][1:2], Rs³ᴰ′[2][1:2])
    P   = P³ᴰ[SOneTo(2), SOneTo(2)]
    # the 2D basis is now Niggli-reduced, but we want to maintain its original handedness: 
    # even though P³ᴰ is guaranteed to have det(P³ᴰ) = 1, the same is not necessarily true
    # for P = P³ᴰ[1:2,1:2], if P³ᴰ[3,3] == -1; in that case, we have to fix it - doing so is
    # not trickier than one might naively think: one cannot simply rotate the basis (this
    # could correspond to rotating the lattice, which might not preserve the lattice), nor 
    # simply swap the signs of `Rs′[1]` or `Rs′[2]`, since that change their mutual angles -
    # nor even generically swap `Rs[1]` to `Rs′[2]`; instead, the appropriate change depends
    # on the relative lengths and angles of `Rs′[1]` & `Rs′[2]`
    if P³ᴰ[3,3] == -1
        ϵ = rtol * abs(volume(Rs′))^(1/2)
        local P_flip :: SMatrix{2,2,Int,4}
        if abs(dot(Rs′[1], Rs′[1]) - dot(Rs′[2], Rs′[2])) < ϵ # Rs′[1] ≈ Rs′[2]
            # then we can just swap 𝐚 and 𝐛 to get a right-handed basis, without worrying
            # about breaking the Niggli-rule |𝐚|≤|𝐛|
            P_flip = @SMatrix [0 1; 1 0]
        else # norm(Rs′[1]) > norm(Rs′[2])
            # we either swap the sign of 𝐛 or check if 𝐛=𝐚-𝐛 gives a bigger ∠(𝐚, 𝐛) while
            # having the same length as |𝐛| (cf. Sec. 9.3.1 ITA5)
            candidate1 = -Rs′[2]
            candidate2 = Rs′[1] - Rs′[2]
            if abs(dot(candidate1, candidate1) - dot(candidate2, candidate2)) < ϵ
                # pick the candidate that maximizes the term |π/2 - ∠(𝐚,𝐛)| (Buerger 
                # condition (iv))
                γ1 = signed_angle²ᴰ(Rs′[1], candidate1) # ∠(𝐚, candidate1)
                γ2 = signed_angle²ᴰ(Rs′[1], candidate2) # ∠(𝐚, candidate2)
                if abs(π/2 - γ1) < abs(π/2 - γ2)  # pick `candidate2` for 𝐛
                    P_flip = @SMatrix [1 0; 1 -1]
                else                              # pick `candidate1` for 𝐛
                    P_flip = @SMatrix [1 0; 0 -1]
                end
            else # use `candidate1` for 𝐛
                P_flip = @SMatrix [1 0; 0 -1]
            end
        end           
        P *= P_flip
        Rs′ = transform(Rs′, P_flip)      
    end
    if det(P) < 0 # basis change did not preserve handedness, contrary to intent
        error("2D Niggli reduction failed to produce a preserve handedness")
    end
    if P³ᴰ[1,3] ≠ 0 || P³ᴰ[2,3] ≠ 0 || P³ᴰ[3,1] ≠ 0 || P³ᴰ[3,2] ≠ 0
        error(lazy"interim 3D transformation has unexpected nonzero elements: P³ᴰ=$P³ᴰ")
    end
    return Rs′, P
end
# computes the angle between 2D vectors a & b in [-π,π)
signed_angle²ᴰ(a::StaticVector{2}, b::StaticVector{2}) = atan(a[1]*b[2]-a[2]*b[1], dot(a,b))

# ---------------------------------------------------------------------------------------- #
# Reciprocal lattice vector by Niggli reduction of the direct lattice first
function nigglibasis(Gs :: ReciprocalBasis{D}; kws...) where D
    Rs = DirectBasis(reciprocalbasis(Gs).vs) # TODO: replace by `dualbasis` when implemented
    Rs′, P = nigglibasis(Rs; kws...)
    return reciprocalbasis(Rs′), P
end
