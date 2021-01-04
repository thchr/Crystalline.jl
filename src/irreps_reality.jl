"""
    realify(lgirs::AbstractVector{<:LGIrrep}, verbose::Bool=false)
                                                        --> AbstractVector{<:LGIrrep}

From `lgirs`, a vector of `LGIrrep`s, determine the associated (gray) co-representations,
i.e. the "real", or "physical" irreps that are relevant in scenarios with time-reversal
symmetry.

For `LGIrrep` that are `REAL`, or that characterize a k-point 𝐤 which is not
equivalent to -𝐤 (i.e. its star does not include both 𝐤 and -𝐤; equivalently, the little
group includes time-reversal symmetry), the associated co-representations are just the 
original irreps themselves. 
For `PSEUDOREAL` and `COMPLEX` `LGIrrep`s where ±𝐤 are equivalent, the
associated co-representations are built from pairs of irreps that "stick" together. This
method computes this pairing and sets the `LGIrrep` field `iscorep` to true, to indicate
that the resulting "paired irrep" (i.e. the co-representation) should be doubled with 
itself (`PSEUDOREAL` reality) or its complex conjugate (`COMPLEX` reality).

### Background
For background, see p. 650-652 (and 622-626 for point groups) in Bradley & Cracknell's book.
Their discussion is for magnetic groups (the "realified" irreps are really correspond to
co-representations of "gray" magnetic groups). 
Cornwell's book also does a good job of explicating this, as does Inui (p. 296-299).

### Keyword arguments
- `verbose::Bool`: if set to `true`, prints details about mapping from small irrep to small
corep for each `LGIrrep` (default: `false`).
"""
function realify(lgirs::AbstractVector{LGIrrep{D}}, verbose::Bool=false) where D
    Nirr = length(lgirs)
    lg = group(first(lgirs))
    kv = kvec(lg) # must be the same for all irreps in list
    αβγ    = D == length(TEST_αβγ) ? TEST_αβγ : TEST_αβγ[OneTo(D)]
    kv_αβγ = kv(αβγ)
    sgnum = num(lg)
    lgops = operations(lg)
    Nops = order(lg) # order of little group (= number of operations)

    cntr = centering(sgnum, D)
    sgops = operations(spacegroup(sgnum, D))

    verbose && print(klabel(lg), " │ ")

    # Check if -𝐤 is in the star of 𝐤, or if 𝐤 is equivalent to -𝐤: 
    # if so, TR is an element of the little group; if not, it isn't 
    # ║ 𝐑𝐞𝐚𝐬𝐨𝐧: if there is an element g of the (unitary) 𝑠𝑝𝑎𝑐𝑒 group G   
    # ║   that takes 𝐤 to -𝐤 mod 𝐆, then (denoting the TR element by Θ, 
    # ║   acting as θ𝐤 = -𝐤) the antiunitary element θg will take 𝐤 to  
    # ║   𝐤 mod 𝐆, i.e. θg will be an element of the little group of 𝐤
    # ║   M(k) associated with the **gray** space group M ≡ G + θG.
    # ║   Conversely, if no such element g exists, there can be no anti-
    # ║   unitary elements in the little group derived from M; as a result, 
    # ║   TR is not part of the little group and so does not modify its 
    # ║   small irreps (called "co-reps" for magnetic groups).
    # ║   There can then only be type 'x' degeneracy (between 𝐤 and -𝐤)
    # ║   but TR will not change the degeneracy at 𝐤 itself. Cornwall
    # ║   refers to this as "Case (1)" on p. 151.
    if !isapproxin(-kv, kstar(sgops, kv, cntr), cntr; atol=DEFAULT_ATOL)
        corep_idxs = [[i] for i in OneTo(Nirr)] # TR ∉ M(k) ⇒ smalls irrep (... small co-reps) not modified by TR
        verbose && println(klabel(lg), "ᵢ ∀i (type x) ⇒  no additional degeneracy (star{k} ∌ -k)")

    else
        # Test if 𝐤 is equivalent to -𝐤, i.e. if 𝐤 = -𝐤 + 𝐆
        k_equiv_kv₋ = isapprox(-kv, kv, cntr; atol=DEFAULT_ATOL)

        # Find an element in G that takes 𝐤 → -𝐤 (if 𝐤 is equivalent to -𝐤, 
        # then this is just the unit-element I (if `sgops` is sorted conven-
        # tionally, with I first, this is indeed what the `findfirst(...)`  
        # bits below will find)
        if !k_equiv_kv₋
            g₋ = sgops[findfirst(g-> isapprox(g∘kv, -kv, cntr; atol=DEFAULT_ATOL), sgops)]
        else
            # This is a bit silly: if k_equiv_kv₋ = true, we will never use g₋; but I'm not sure if 
            # the compiler will figure that out, or if it will needlessly guard against missing g₋?
            g₋ = SymOperation{D}(hcat(I, zeros(D))) # ... the unit element I
        end

        # -𝐤 is part of star{𝐤}; we infer reality of irrep from ISOTROPY's data (could also 
        # be done using `calc_reality(...)`) ⇒ deduce new small irreps (... small co-reps)
        corep_idxs = Vector{Vector{Int64}}()
        skiplist = Vector{Int64}()
        for (i, lgir) in enumerate(lgirs)
            if i ∈ skiplist; continue; end # already matched to this irrep previously; i.e. already included now
            iscorep(lgir) && throw(DomainError(iscorep(lgir), "should not be called with LGIrreps that have iscorep=true"))
            verbose && i ≠ 1 && print("  │ ")

            if reality(lgir) == REAL
                push!(corep_idxs, [i])
                if verbose
                    println(formatirreplabel(label(lgir)), 
                            " (real) ⇒  no additional degeneracy")
                end

            elseif reality(lgir) == PSEUDOREAL
                # doubles irrep on its own
                push!(corep_idxs, [i, i])
                if verbose
                    println(formatirreplabel(label(lgir)^2), 
                            " (pseudo-real) ⇒  doubles degeneracy"); 
                end

            elseif reality(lgir) == COMPLEX
                # In this case, there must exist a "partner" irrep (say, Dⱼ) which is
                # equivalent to the complex conjugate of the current irrep (say, Dᵢ), i.e.
                # an equivalence Dⱼ ∼ Dᵢ*; we next search for this equivalence.
                # When we check for equivalence between irreps Dᵢ* and Dⱼ we must account 
                # for the possibility of a 𝐤-dependence in the matrix-form of the irreps; 
                # specifically, for an element g, its small irrep is
                #     Dᵢ[g] = exp(2πik⋅τᵢ[g])Pᵢ[g],
                # where, crucially, for symmetry lines, planes, and general points 𝐤 depends
                # on (one, two, and three) free parameters (α,β,γ).
                # Thus, for equivalence of irreps Dᵢ* and Dⱼ we require that
                #     Dᵢ*[g] ~ Dⱼ[g]       ∀g ∈ G(k)
                #  ⇔ exp(-2πik⋅τᵢ[g])Pᵢ*[g] ~ exp(2πik⋅τⱼ[g])Pⱼ[g]
                # It seems rather tedious to prove that this is the case for all 𝐤s along a
                # line/plane (α,β,γ). Rather than attempt this, we simply test against an
                # arbitrary value of (α,β,γ) [superfluous entries are ignored] that is
                # non-special (i.e. ∉ {0,0.5,1}); this is `αβγ`.

                # Characters of the conjugate of Dᵢ, i.e. tr(Dᵢ*) = tr(Dᵢ)*
                θχᵢ = conj.(characters(lgir, αβγ))
                
                # Find matching complex partner
                partner = 0
                for j = i+1:Nirr
                    # only check if jth irrep has not previously matched and is complex
                    if j ∉ skiplist && reality(lgirs[j]) == COMPLEX

                        # Note that we require only equivalence of Dᵢ* and Dⱼ; not equality.
                        # Cornwell describes (p. 152-153 & 188) a neat trick for checking this
                        # efficiently: specifically, Dᵢ* and Dⱼ are equivalent irreps if
                        #     χⁱ(g)* = χʲ(g₋⁻¹gg₋) ∀g ∈ G(k)
                        # with g₋ an element of G that takes 𝐤 to -𝐤, and where χⁱ (χʲ) denotes
                        # the characters of the respective irreps.
                        χⱼ = characters(lgirs[j], αβγ)
                        match = true
                        for n in OneTo(Nops)
                            if k_equiv_kv₋ # 𝐤 = -𝐤 + 𝐆 ⇒ g₋ = I (the unit element), s.t. g₋⁻¹gg₋ = I⁻¹gI = g    (Cornwall's case (3))
                                χⱼ_g₋⁻¹gg₋ = χⱼ[n]
                            else           # 𝐤 not equivalent to -𝐤, i.e. 𝐤 ≠ -𝐤 + 𝐆, but -𝐤 is in the star of 𝐤 (Cornwall's case (2))
                                g₋⁻¹gg₋ = compose(compose(inv(g₋), lgops[n], false), g₋, false)
                                n′, Δw = findequiv(g₋⁻¹gg₋, lgops, cntr)
                                χⱼ_g₋⁻¹gg₋ = cis(2π*dot(kv_αβγ, Δw)) .* χⱼ[n′] # cis(x) = exp(ix)
                            end
                            
                            match = isapprox(θχᵢ[n], χⱼ_g₋⁻¹gg₋; atol=DEFAULT_ATOL)
                            if !match # ⇒ not a match
                                break
                            end
                        end

                        if match # ⇒ a match
                            partner = j
                            if verbose; 
                                println(formatirreplabel(label(lgir)*label(lgirs[j])), " (complex) ⇒  doubles degeneracy")
                            end
                        end
                    end
                end
                partner === 0 && throw(ErrorException("Didn't find a matching complex partner for $(label(lgir))"))
                push!(skiplist, partner)

                push!(corep_idxs, [i, partner])
                
            else
                throw(ArgumentError("Invalid real/pseudo-real/complex reality = $(reality(lgir))"))
            end
        end
    end

    Ncoreps = length(corep_idxs)

    # New small co-rep labels (composite)
    newlabs = Tuple(join(label(lgirs[i]) for i in corep_idxs[i′]) for i′ in OneTo(Ncoreps))

    # Build a vector of "new" small irreps (small co-reps), following B&C p. 616 & Inui p.
    # 298-299. For pseudo-real and complex co-reps, we set a flag `iscorep = true`, to
    # indicate to "evaluation" methods, such as `irreps(::LGIrrep)`, that a diagonal
    # "doubling" is required (see below).
    lgirs′ = Vector{LGIrrep{D}}(undef, Ncoreps)
    for i′ in OneTo(Ncoreps)
        idxs = corep_idxs[i′]
        if length(idxs) == 1      # ⇒ real or type x (unchanged irreps)
            lgirs′[i′] = lgirs[idxs[1]] # has iscorep = false flag set already

        elseif idxs[1] == idxs[2] # ⇒ pseudoreal     ("self"-doubles irreps)
            # The resulting co-rep of a pseudo-real irrep Dᵢ is
            #   D = diag(Dᵢ, Dᵢ)
            # See other details under complex case.
            lgir = lgirs[idxs[1]]
            blockmatrices = _blockdiag2x2.(lgir.matrices)
            lgirs′[i′] = LGIrrep{D}(newlabs[i′], lg, blockmatrices, lgir.translations,
                                    PSEUDOREAL, true)
            
        else                      # ⇒ complex        (doubles irreps w/ complex conjugate)
            # The co-rep of a complex irreps Dᵢ and Dⱼ is 
            #   D = diag(Dᵢ, Dⱼ)
            # where we know that Dⱼ ∼ Dᵢ*. Note that this is _not_ generally the same as
            # diag(Dⱼ, Dⱼ*), since we have only established that Dⱼ ∼ Dᵢ*, not Dⱼ = Dᵢ*.
            # Note also that we require the Dᵢ and Dⱼ irreps to have identical free phases,
            # i.e. `translations` fields, so that the overall irrep "moves" with a single
            # phase factor in k-space - we check for that explicitly for now, to be safe.
            τsᵢ = lgirs[idxs[1]].translations; τsⱼ = lgirs[idxs[2]].translations
            @assert τsᵢ == τsⱼ
            blockmatrices = _blockdiag2x2.(lgirs[idxs[1]].matrices, lgirs[idxs[2]].matrices)
            
            lgirs′[i′] = LGIrrep{D}(newlabs[i′], lg, blockmatrices, τsᵢ, COMPLEX, true)
        end
    end
    
    return lgirs′
end

# returns the block diagonal matrix `diag(A1, A2)` (and assumes identically sized and
# square `A1` and `A2`).
function _blockdiag2x2(A1::AbstractMatrix{T}, A2::AbstractMatrix{T}) where T
    n = LinearAlgebra.checksquare(A1)
    LinearAlgebra.checksquare(A2) == n || throw(DimensionMismatch())

    B = zeros(T, 2n, 2n)
    for i in OneTo(n)
        i′ = i+n
        @inbounds for j in OneTo(n)
            j′ = j+n
            B[i,j]   = A1[i,j]
            B[i′,j′] = A2[i,j]
        end
    end
    return B
end
# returns the block diagonal matrix `diag(A, A)` (and assumes square `A`)
function _blockdiag2x2(A::AbstractMatrix{T}) where T
    n = LinearAlgebra.checksquare(A)

    B = zeros(T, 2n, 2n)
    for i in OneTo(n)
        i′ = i+n
        @inbounds for j in OneTo(n)
            j′ = j+n
            aᵢⱼ = A[i,j]
            B[i,j]   = aᵢⱼ
            B[i′,j′] = aᵢⱼ
        end
    end
    return B
end

@doc raw"""
    calc_reality(lgir::LGIrrep, 
                 sgops::AbstractVector{SymOperation{D}},
                 αβγ::Union{Vector{<:Real},Nothing}=nothing) --> ::(Enum Reality)

Compute and return the reality of a `lgir::LGIrrep` using the Herring criterion.

The computed value is one of three integers in ``{1,-1,0}``.
In practice, this value is returned via a member of the Enum `Reality`, which has instances
`REAL = 1`, `PSEUDOREAL = -1`, and `COMPLEX = 0`.

## Optional arguments
As a sanity check, a value of `αβγ` can be provided to check for invariance along a symmetry
symmetry line/plane/general point in k-space. The reality must be invariant to this choice.

## Note 
The provided space group operations `sgops` **must** be the set reduced by primitive
translation vectors; i.e. using `spacegroup(...)` directly is **not** allowable in general
(since the irreps we reference only include these "reduced" operations). This reduced set
of operations can be obtained e.g. from the Γ point irreps of ISOTROPY's dataset, or
alternatively, from `reduce_ops(spacegroup(...), true)`.

## Implementation
The Herring criterion evaluates the following sum

``[∑ χ({β|b}²)]/[g_0/M(k)]``

over symmetry operations ``{β|b}`` that take ``k → -k``. Here ``g_0`` is the order of the
point group of the space group and ``M(k)`` is the order of star(``k``) [both in a primitive
basis].

See e.g. Cornwell, p. 150-152 & 187-188 (which we mainly followed), Inui Eq. (13.48), 
Dresselhaus, p. 618, or [Herring's original paper](https://doi.org/10.1103/PhysRev.52.361).
"""
function calc_reality(lgir::LGIrrep, 
                      sgops::AbstractVector{SymOperation{D}}, 
                      αβγ::Union{Vector{<:Real},Nothing}=nothing) where D
    iscorep(lgir) && throw(DomainError(iscorep(lgir), "method should not be called with LGIrreps where iscorep=true"))
    lgops = operations(lgir)
    kv = kvec(lgir)
    kv₋ = -kv
    cntr = centering(num(lgir), D)
    Ds = irreps(lgir, αβγ) # irrep matrices
    kv_αβγ = kv(αβγ)

    s = zero(ComplexF64)
    for op in sgops
        if isapprox(op∘kv, kv₋, cntr, atol=DEFAULT_ATOL) # check if op∘k == -k; if so, include in sum
            op² = compose(op, op, false) # this is op∘op, _including_ trivial lattice translation parts
            # find the equivalent of `op²` in `lgops`; this may differ by a number of 
            # primitive lattice vectors `w_op²`; the difference must be included when 
            # we calculate the trace of the irrep 𝐃: the irrep matrix 𝐃 is ∝exp(2πi𝐤⋅𝐭)
            idx_of_op²_in_lgops, Δw_op² = findequiv(op², lgops, cntr)
            ϕ_op² = cis(2π*dot(kv_αβγ, Δw_op²)) # phase accumulated by "trivial" lattice translation parts [cis(x) = exp(ix)]
            χ_op² = ϕ_op²*tr(Ds[idx_of_op²_in_lgops]) # χ(op²)

            s += χ_op²
        end
    end

    pgops = pointgroup(sgops) # point group assoc. w/ space group
    g₀ = length(pgops) # order of pgops (denoted h, or macroscopic order, in Bradley & Cracknell)
    Mk = length(kstar(pgops, kv, cntr)) # order of star of k (denoted qₖ in Bradley & Cracknell)
    normalization = convert(Int, g₀/Mk) # order of G₀ᵏ; the point group derived from the little group Gᵏ (denoted b in Bradley & Cracknell; [𝐤] in Inui)
    
    # s = ∑ χ({β|b}²) and normalization = g₀/M(k) in Cornwell's Eq. (7.18) notation
    type_float = real(s)/normalization
    type       = round(Int8, type_float)
    # check that output is a valid: real integer in (0,1,-1)
    isapprox(imag(s),    0.0,  atol=DEFAULT_ATOL) || _throw_reality_not_real(s)
    isapprox(type_float, type, atol=DEFAULT_ATOL) || _throw_reality_not_integer(real(s))
    
    return Reality(type) # return [∑ χ({β|b}²)]/[g₀/M(k)]
end

# Frobenius-Schur criterion for point group irreps (Inui p. 74-76):
#   |g|⁻¹∑ χ(g²) = {1 (≡ real), -1 (≡ pseudoreal), 0 (≡ complex)}
function calc_reality(pgir::PGIrrep)
    χs = characters(pgir)
    pg = group(pgir)

    s = zero(eltype(χs))
    for op in pg
        op² = op∘op
        idx = findfirst(≈(op²), pg)
        idx == nothing && error("unexpectedly did not find group element for op²")

        s += χs[idx]
    end

    type_float = real(s)/order(pg)
    type      = round(Int8, type_float)
    isapprox(imag(s),    0.0,  atol=DEFAULT_ATOL) || _throw_reality_not_real(s)
    isapprox(type_float, type, atol=DEFAULT_ATOL) || _throw_reality_not_integer(real(s))

    return Reality(type) # return |g|⁻¹∑ χ(g²)
end

@noinline _throw_reality_not_integer(x) = error("Criterion must produce an integer; obtained non-integer value = $(x)")
@noinline _throw_reality_not_real(x)    = error("Criterion must yield a real value; obtained complex value = $(x)")