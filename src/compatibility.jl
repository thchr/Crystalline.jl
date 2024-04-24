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

    compat_idxs = Vector{Int}()
    @inbounds for (idx′, kv′) in enumerate(kvs′)
        isspecial(kv′) && continue # must be a line/plane/general point to match a special point kv
        is_compatible(kv, kv′).bool && push!(compat_idxs, idx′) 
    end

    return compat_idxs
end

#=
function compatibility_matrix(BRS::BandRepSet)
    lgirs_in, lgirs_out = matching_lgirreps(BRS::BandRepSet)
    for (iᴳ, Dᴳᵢ) in enumerate(lgirs_in)         # super groups
        for (jᴴ, Dᴴⱼ) in enumerate(lgirs_out)    # sub groups
            # we ought to only check this on a per-kvec basis instead of 
            # on a per-lgir basis to avoid redunant checks, but can't be asked...
            compat_bool, αβγ, αβγ′, G  = is_compatible(position(Dᴳᵢ), position(Dᴴⱼ))
            # TODO: Use associated (αβγ, αβγ′, G) "matching" values that makes kvⱼ and kvᵢ 
            #       and compatible; use to get correct lgirs at their "intersection".
            if compat_bool
                nᴳᴴᵢⱼ = subduction_count(Dᴳᵢ, Dᴴⱼ, αβγ)
                if !iszero(nᴳᴴᵢⱼ)
                    # TODO: more complicated than I thought: have to match across different
                            special lgirreps.
                end 
            end
        end
    end
end
=#