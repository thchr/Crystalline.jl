"""
    primitivize(irs::Collection{<:Union{LGIrrep, SiteIrrep}}, [cntr::Char])
                                            -> Collection{<:Union{LGIrrep, SiteIrrep}}

Given a collection of irreps, whose underlying group is specified in a conventional basis
(as e.g., returned by [`lgirreps`](@ref) or [`siteirreps`](@ref)), return a new collection
of irreps, referenced relative to a group specified in a primitive basis.

The centering symbol `cntr` is optional and will be inferred from the `irs` if not
explicitly provided.

The returned "primitivized" little group irreps may share data with the input `irs`:
if subsequent mutation is desired, use `deepcopy` on the returned collection.
"""
function primitivize(
    irs::Collection{T},
    cntr::Char = centering(num(irs), D)
) where T <: Union{LGIrrep{D}, SiteIrrep{D}} where D
    if ((D == 3 && cntr == 'P') || (D ≠ 3 && cntr == 'p'))
        return irs # already primitive; return as-is
    end

    # not already primitive; primitivize underlying group elements & reconstruct
    g = group(irs) # little group or site group
    # NB: for both `LGIrrep` & `SiteIrrep` input, `g` will not contain centering-copies
    #     centering copies of the same operations; so `primitivize(g)` and `g` contain
    #     the same operations (in the same order) - just in different bases. Cf. the fact
    #     that `lgirreps(...)[klab]` and `siteirreps(...)` only return irreps sampled at
    #     non-centering-repeated group elements
    g′ = if T === LGIrrep{D}
        primitivize(g::LittleGroup{D}, #=modw: do not reduce translations=# false)
    elseif T === SiteIrrep{D}
        primitivize(g::SiteGroup{D})
    else
        error("unreachable")
    end

    irs′ = Vector{T}(undef, length(irs))
    for (i, ir) in enumerate(irs)
        lgir′ = _rebuild_irrep_with_modified_group(ir, g′)
        irs′[i] = lgir′
    end
    return Collection(irs′)
end
function primitivize(lgirsd::Dict{String, Collection{LGIrrep{D}}}) where D
    cntr = centering(num(first(values(lgirsd))), D)
    return Dict(klab => primitivize(lgirs, cntr) for (klab, lgirs) in lgirsd)
end

function _rebuild_irrep_with_modified_group(ir::LGIrrep{D}, g′::LittleGroup{D}) where D
    # we have to also update the τᵢ = `ir.translations[i]` field, since if `g′` now refers to a
    # a k-point in a new basis, say, `k′`, while the original `g` referred to `k`, we must
    # ensure that k′⋅τᵢ′ = k⋅τᵢ, so the k-τ products are invariant (→ invariant phase factors
    # in the representation matrices)
    k′ = position(g′)
    k = position(ir)
    τs = ir.translations
    τs′ = if k == k′ || all(iszero, τs) # unchanged momentum or zero-translates: keep `τs`
        τs
    else # changed momentum: convert τᵢ as well
        # use invariance of dot-product under transformation:
        #    k′ = Pᵀk ⇒ k′⋅τ′ = k⋅τ ⇒ τ′ = P⁻¹τ
        P = primitivebasismatrix(centering(num(ir)), Val(D))
        [P\τ for τ in τs]
    end :: typeof(τs)


    # TODO: This is wrong if unless `ir.translations` is all zeros. Remember: eventually,
    #       the translations (τ) come into play as phase factors `cispi(2k⋅τ)` (nevermind
    #       sign now): so if we have updated `k` to a primitive counterpart `k′`, then we
    #       need to update the translations to `τ′` such that `cispi(2k⋅τ) == cispi(2k′⋅τ′)`
    return LGIrrep{D}(ir.cdml, g′, ir.matrices, τs′, ir.reality, ir.iscorep)
end
function _rebuild_irrep_with_modified_group(ir::SiteIrrep{D}, g′::SiteGroup{D}) where D
    return SiteIrrep{D}(ir.cdml, g′, ir.matrices, ir.reality, ir.iscorep, ir.pglabel)
end

"""
    primitivize(brs::Collection{<:NewBandRep}, [cntr::Char]) -> Collection{<:NewBandRep}

Analogous to `primitivize(::Collection{<:Union{LGIrrep, SiteIrrep}}, ::Char)` but for
band representations.

Primitivizes the groups associated with both the underlying little group irreps and the
site irreps.
"""
function primitivize(
    brs::Collection{NewBandRep{D}},
    cntr::Char = centering(num(brs), D)
) where D
    # --- early termination; don't need to do anything if already primitive ---
    ((D == 3 && cntr == 'P') || (D ≠ 3 && cntr == 'p')) && return brs

    # --- primitivize little group irreps ---
    # NB: all elements of `brs` point to the same set of irreps, by assumption
    lgirsv = irreps(brs)
    lgirsv′ = Vector{Collection{LGIrrep{D}}}(undef, length(lgirsv))
    for (i, lgirs) in enumerate(lgirsv)
        lgirsv′[i] = primitivize(lgirs, cntr)
    end

    # --- primitivize siteirreps & update each band rep ---
    vs′ = Vector{NewBandRep{D}}(undef, length(brs))
    for (i, br) in enumerate(brs)
        siteg = group(br)
        siteg′ = primitivize(siteg)
        siteir = br.siteir
        siteir′ = SiteIrrep{D}(siteir.cdml, siteg′, siteir.matrices, siteir.reality,
                               siteir.iscorep, siteir.pglabel)
        n = br.n
        n′ = SymmetryVector{D}(lgirsv′, multiplicities(n), occupation(n))
        br′ = NewBandRep{D}(siteir′, n′, br.timereversal, br.spinful)
        vs′[i] = br′
    end
    brs′ = Collection(vs′)

    return brs′
end

isspecial(c::Collection{<:AbstractIrrep}) = isspecial(first(c))