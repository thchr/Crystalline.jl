@doc raw"""
    subduction_count(Dᴳᵢ, Dᴴⱼ[, αβγᴴⱼ]) --> Int

For two groups ``G`` and ``H``, where ``H`` is a subgroup of ``G``, i.e. ``H<G``, with
associated irreducible representations `Dᴳᵢ` = ``D^G_i(g)`` and `Dᴴⱼ` = ``D^H_j(g)`` over
operations ``g∈G`` and ``h∈H<G``, compute the compatibility relation between the two irreps
from the subduction reduction formula (or "magic" formula/Schur orthogonality relation),
returning how many times ``n^{GH}_{ij}`` the subduced representation ``D^G_i↓H`` contains 
the irrep ``D^H_j``; in other words, this gives the compatibility between the two irreps.

Optionally, a vector `αβγᴴⱼ` may be provided, to evaluate the characters/irreps 
of `Dᴳᵢ` at a concrete value of ``(α,β,γ)``. This is e.g. meaningful for `LGIrrep`s at
non-special **k**-vectors. Defaults to `nothing`.

The result is computed using the reduction formula [see e.g. Eq. (15) of
[arXiv:1706.09272v2](https://arxiv.org/abs/1706.09272)]:

``n^{GH}_{ij} = |H|^{-1} \sum_h \chi^G_i(h)\chi^H_j(h)^*``

## Example
Consider the two compatible **k**-vectors Γ (a point) and Σ (a line) in space group 207:
```jl
lgirsd  = lgirreps(207, Val(3));
Γ_lgirs = lgirsd["Γ"]; # at Γ ≡ [0.0, 0.0, 0.0]
Σ_lgirs = lgirsd["Σ"]; # at Σ ≡ [α, α, 0.0]
```
We can test their compatibility via:
```jl
[[subduction_count(Γi, Σj) for Γi in Γ_lgirs] for Σj in Σ_lgirs]
> # Γ₁ Γ₂ Γ₃ Γ₄ Γ₅
>  [ 1, 0, 1, 1, 2] # Σ₁
>  [ 0, 1, 1, 2, 1] # Σ₂
```
With following enterpretatation for compatibility relations between irreps at Γ and Σ:

| Compatibility relation | Degeneracies |
|------------------------|--------------|
| Γ₁ → Σ₁                | 1 → 1        |
| Γ₂ → Σ₂                | 1 → 1        |
| Γ₃ → Σ₁ + Σ₂           | 2 → 1 + 1    |
| Γ₄ → Σ₁ + 2Σ₂          | 3 → 1 + 2    |
| Γ₅ → 2Σ₁ + Σ₂          | 3 → 2 + 1    |

where, in this case, all the small irreps are one-dimensional.
"""
function subduction_count(Dᴳᵢ::T, Dᴴⱼ::T, 
                          αβγᴴⱼ::Union{<:AbstractVector{<:Real},Nothing}=nothing
                          ) where T<:AbstractIrrep
    # find matching operations between H & G and verify that H<G 
    boolsubgroup, idxsᴳ²ᴴ = _findsubgroup(operations(Dᴳᵢ), operations(Dᴴⱼ))
    !boolsubgroup && throw(DomainError("Provided irreps are not H<G subgroups"))

    # compute characters 
    # TODO: Care should be taken that the irreps 
    # actually can refer to identical k-points; that should be a check 
    # too, and then we should make sure that the characters are actually
    # evaluated at that KVec
    χᴳᵢ = characters(Dᴳᵢ)
    χᴴⱼ = characters(Dᴴⱼ, αβγᴴⱼ)

    # compute number of times that Dᴴⱼ occurs in the reducible 
    # subduced irrep Dᴳᵢ↓H
    s = zero(ComplexF64)
    @inbounds for (idxᴴ, χᴴⱼ′) in enumerate(χᴴⱼ)
        s += χᴳᵢ[idxsᴳ²ᴴ[idxᴴ]]*conj(χᴴⱼ′)
    end
    (abs(imag(s)) > DEFAULT_ATOL) && error("unexpected finite imaginary part")
    nᴳᴴᵢⱼ_float = real(s)/order(Dᴴⱼ)
    # account for the fact that the orthogonality relations are changed for coreps; seems
    # only the subgroup is important to include here - not exactly sure why
    nᴳᴴᵢⱼ_float /= corep_orthogonality_factor(Dᴴⱼ)
    nᴳᴴᵢⱼ = round(Int, nᴳᴴᵢⱼ_float)
    abs(nᴳᴴᵢⱼ - nᴳᴴᵢⱼ_float) > DEFAULT_ATOL && error("unexpected non-integral compatibility count")
    
    return nᴳᴴᵢⱼ
end

"""
$(TYPEDSIGNATURES)
"""
function find_compatible(kv::KVec{D}, kvs′::AbstractVector{KVec{D}}) where D
    isspecial(kv) || throw(DomainError(kv, "input kv must be a special k-point"))

    compat_idxs = Int[]
    @inbounds for (idx′, kv′) in enumerate(kvs′)
        isspecial(kv′) && continue # must be a line/plane/general point to match a special point kv
        can_intersect(kv, kv′).bool && push!(compat_idxs, idx′) 
    end

    return compat_idxs
end

"""
$(TYPEDSIGNATURES)

Check whether two `AbstractVec`s `v` and `v′` can intersect, i.e., whether there exist free
parameters such that they are equivalent modulo an integer lattice vector.

Returns a `NamedTuple` `(; bool, αβγ, αβγ′, L)`. If `bool = true`, `kv` and `kv′` are
compatible in the sense that `v(αβγ) == v′(αβγ′) + L` if `bool == true` where `L` is
an integer-valued lattice vector.
If `bool = false`, they are incompatible (and zero-valued vectors are returned for `αβγ`,
`αβγ′`, and `L`).

## Extended help

- The keyword argument `atol` (default, $DEFAULT_ATOL) specifies the absolute tolerance for
  the comparison.
- If both `v` and `v′` are not special, i.e., both have free parameters, the intersection
  point may not be unique (e.g., for co-linear `v` and `v′`).
- The implementation currently only checks the immediately adjacent lattice vectors for
  equivalence; if there is equivalence, but the the required elements of `L` would have
  `|Lᵢ| > 1`, the currently implementation will not identify the equivalence.
- This operation is usually only meaningful if the bases of `kv` and `kv′` agree and are
  primitive.
"""
function can_intersect(v::T, v′::T; atol::Real=DEFAULT_ATOL) where T<:AbstractVec{D} where D
    # check if solution exists to [A] v′ = v(αβγ) or [B] v′(αβγ′) = v(αβγ) by solving
    # a least squares problem and then checking if it is a strict solution. Details:
    #   Let v(αβγ) = v₀ + V*αβγ and v′(αβγ′) = v₀′ + V′*αβγ′
    #   [A] v₀′ = v₀ + V*αβγ             ⇔  V*αβγ = v₀′-v₀
    #   [B] v₀′ + V′*αβγ′ = v₀ + V*αβγ   ⇔  V*αβγ - V′*αβγ′ = v₀′-v₀  
    #                                    ⇔  hcat(V,-V′)*vcat(αβγ,αβγ′) = v₀′-v₀
    # these equations can always be solved in the least squares sense using the
    # pseudoinverse; we can then subsequently check if the residual of that solution is in
    # fact zero, in which can the least squares solution is a "proper" solution, signaling
    # that `v` and `v′` can intersect (at the found values of `αβγ` and `αβγ′`)
    Δcnst = constant(v′) - constant(v)
    if isspecial(v′)
        Δfree = free(v)                                     # D×D matrix
        return _can_intersect_equivalence_check(Δcnst, Δfree, atol, false)
    elseif isspecial(v)
        Δfree = -free(v′)                                   # D×D matrix
        return _can_intersect_equivalence_check(Δcnst, Δfree, atol, true)
    else                     # neither `v′` nor `v` are special
        Δfree = hcat(free(v), -free(v′))                    # D×2D matrix
        return _can_intersect_equivalence_check(Δcnst, Δfree, atol, false)
    end
    # NB: the above seemingly trivial splitting of return statements is intentional & to
    #     avoid type-instability (because the type of `Δfree` differs in the brances)
end

function _can_intersect_equivalence_check(Δcnst::StaticVector{D}, Δfree::StaticMatrix{D},
                                          atol::Real, inverted_order::Bool=false) where D
    # to be safe, we have to check for equivalence between `v` and `v′` while accounting
    # for the fact that they could differ by a lattice vector; in practice, for the wyckoff
    # listings that we have have in 3D, this seems to only make a difference in a single 
    # case (SG 130, wyckoff position 8f) - but there the distinction is actually needed
    Δfree⁻¹ = pinv(Δfree)
    for _L in Iterators.product(ntuple(_->(0, -1, 1), Val(D))...) # loop over adjacent lattice vectors
        L = SVector{D,Int}(_L)
        Δcnst_plus_L = Δcnst + L
        _αβγ = Δfree⁻¹*Δcnst_plus_L   # either `D`-dim `αβγ` or `2D`-dim `vcat(αβγ, αβγ′)`
        Δ = Δcnst_plus_L - Δfree*_αβγ # residual of least squares solve
        if norm(Δ) < atol
            if length(_αβγ) == D
                αβγ  = _αβγ
                αβγ′ = zero(αβγ)
                if inverted_order
                    αβγ, αβγ′ = αβγ′, αβγ
                end
            else # size(_αβγ, 2) == 2D
                αβγ  = _αβγ[SOneTo{D}()]
                αβγ′ = _αβγ[StaticArrays.SUnitRange{D+1,D}()]
            end
            return (; bool=true, αβγ=αβγ, αβγ′=αβγ′, L=L)
        end
    end
    sentinel = zero(SVector{D, Int})
    return (; bool=false, αβγ=sentinel, αβγ′=sentinel, L=sentinel)
end

"""
    remap_to_kstar(
        lgirs::AbstractVector{LGIrrep{D}},
        kv′::KVec{D},
        coset_representatives::AbstractVector{SymOperation{D}}
        ) --> Collection{LGIrrep{D}}

Given an set of `LGIrrep`s `lgirs` defined at a **k**-vector `kv`, remap the irrep data to
a different **k**-vector `kv′` in the star of `kv`.

The remapping is done by identifying an operation `g` s.t. `kv′ = g * kv` with `g` in the
space group of `lgirs` (more precisely, from among the coset representatives of the little
group in the space group). The original irreps ``D(h)`` with ``h`` in the little group of
`kv` are then transformed according to ``D′(h′) = D(h) = D(g⁻¹h′g)`` with ``h′`` from the
little group of `kv′`.

The coset representatives can be specified as an optional argument, to avoid repeated
recomputation and simplify the associated computation of `g`. The coset representatives
generate the star of `kv`.
"""
function remap_to_kstar(
            lgirs::AbstractVector{LGIrrep{D}},
            kv′::KVec{D},
            coset_representatives::AbstractVector{SymOperation{D}} = 
                        cosets(reduce_ops(spacegroup(num(first(lgirs)), Val{D}()),
                                          centering(num(first(lgirs)))),
                               group(first(lgirs)))
            ) where D
    
    kv = position(first(lgirs))
    kv′ == kv && return Collection(lgirs)

    # feasibility checks
    if freeparams(kv) != freeparams(kv′)
        error(lazy"kv=$kv and kv′=$kv′ do not have the same free parameters")
    end
    special_bool = isspecial(kv)

    # check if `kv′` is in the star of `kv`
    kv_star = map(coset_representatives) do g # compute {star(k)}
        g * kv
    end
    idx = begin
        idx′ = findfirst(≈(kv′), kv_star)
        if !isnothing(idx′)
            # as first priority, we return an exact match if it exists
            idx′
        else
            # otherwise, we look for any compatible match
            findfirst(kv_star) do kv′′
                if special_bool
                    can_intersect(kv′′, kv′).bool
                else
                    # nonspecial pts: check if parallel & possibly separated by a reciprocal
                    # vector; this is slightly more annoying because we might then have to deal
                    # later with a nonzero reciprocal vector
                    (kv′′.free == kv′.free || kv′′.free == -kv′.free) && 
                    all(isinteger, kv′′.cnst - kv′.cnst)
                end
            end
        end
    end
    isnothing(idx) && error(lazy"kv′=$kv′ is not compatible with any element in star(k)=$kv_star")
    g = coset_representatives[something(idx)] # g ∘ kv = kv′

    # remap irrep operations according to D′(h′) = D(h) = D(g⁻¹h′g),  (w/ D referencing
    # `kv′`, and D′ referencing `kv`). I.e., we have h = g⁻¹h′g s.t. h′ = ghg⁻¹
    lg = group(first(lgirs))
    ops′ = similar(operations(lg));
    for (idx, h) in enumerate(lg)
        h′ = compose(g, compose(h, inv(g), #=modτ=#false), #=modτ=#false)
        ops′[idx] = h′
    end
    lg′ = LittleGroup{D}(num(lg), kv′, klabel(lg), ops′)

    # build provisional `LGIrrep`s with above operator-sorting
    lgirs′ = map(lgirs) do lgir
        matrices = [copy(m) for m in lgir.matrices]
        translations = [rotation(g) * τ for τ in lgir.translations] # see (⋆)
        # (⋆) note the rotation of the translation vector: this is necessary since the 
        #     translation part of the original irrep has a form exp(ik⋅τ) - and in the new
        #     setting, it needs to be in the form exp(ik′⋅τ′) but to agree with the original
        #     phase for every free parameter, i.e., we need exp(ik⋅τ)=exp(ik′⋅τ′). Since
        #     k′(G) = g(R)⁻¹ᵀk(G), this translates to the requirement that 
        #     τ′(R) = rotation(g)(R)τ(R).
        LGIrrep{D}(lgir.cdml, lg′, matrices, translations, lgir.reality, lgir.iscorep)
    end

    # if the equivalence between kv and kv′ involves a nonzero G-vector, _AND_ if the 
    # any of `lgirs` depends on k explicitly, i.e., if any τ∈`translations.(lgirs)` are
    # nonzero, and this dependence could impart a dependence on the nonzero G-vector,
    # i.e., if exp(2πik(αβγ)⋅τ) could actually vary with αβγ, we need to revise the irrep
    # data accordingly, to account for the new "starting point"; for now, we don't do this,
    # but just throw an error to be conservative. To figure out if we're in this case, 
    # recall that k(αβγ) = constant(k) + free(k)⋅αβγ, so that the αβγ-dependent part of the
    # phase factor depends on a term αβγ ⋅ free(k)ᵀτ - so if free(k)ᵀτ is zero, we are safe
    # TODO: Try to actually do this without failing; should be possible if we decompose G
    #       into parts that are ∥/⟂ to free(k)ᵀτ (only parallel parts matter)
    ΔG = kv_star[something(idx)].cnst - kv′.cnst
    if (!special_bool && 
        norm(ΔG) > DEFAULT_ATOL &&
        any(lgir -> 
            any(τ -> norm(transpose(free(kv)) * τ) > DEFAULT_ATOL, lgir.translations), 
            lgirs)
        )
        error("nonzero reciprocal vector between nonspecial kv and kv′ requires explicit handling: not yet implemented")
    end

    # the remapped little group `lg′` might have translations that are not in a primitive 
    # setting; that's not great for comparison with tables later, so we need to change the
    # irreps accordingly now. This may affect the irreps of nonsymmorphic groups. We correct
    # by effectivly "multiplying" - via the translations term - with a suitable Bloch phase
    lg′_reduced = reduce_ops(lg′, centering(num(lg′), D))
    for i in eachindex(lg′)
        Δt = translation(lg′_reduced[i]) - translation(lg′[i])
        norm(Δt) > Crystalline.DEFAULT_ATOL || continue
        for lgir′ in lgirs′
            lgir′.translations[i] += Δt
        end
        lg′.operations[i] = lg′_reduced[i]
    end

    return Collection(lgirs′)
end