using LinearAlgebra: dot, \

# The implementation here follows Cano et al., Phys. Rev. B 97, 035139 (2018)
# (https://doi.org/10.1103/PhysRevB.97.035139), specifically, Sections II.C-D

# ---------------------------------------------------------------------------------------- #

"""
    reduce_orbits_and_cosets(siteg::AbstractSiteGroup{D})

For an input site group, provided in conventional coordinates, reduce its cosets such
that the resulting orbit only contains Wyckoff positions that are not equivalent when
viewed in the primitive basis (as determined by the centering type `cntr`).
Additionally, the associated orbit will only contain positions whose _primitive_ 
coefficients lie in [0,1). I.e., the associated orbit lies in the canonical [0,1)ᴰ primitive
unit cell. The coset operations are adjusted accordingly.[^1]

[^1]: Note that this is only guaranteed for orbits without free parameters, since the
      presence of free parameters renders the choice ill-defined (we then only guarantee
      that the constant part of the orbit is in [0,1)ᴰ]).

The reduced site group is returned in a conventional basis, along with the reduced orbits,
also in a conventional basis.
"""
function reduce_orbits_and_cosets(
    siteg::AbstractSiteGroup{D}
) where D

    orbits = parent.(orbit(siteg))
    gαs    = copy(cosets(siteg))
    cntr   = centering(num(siteg), D)

    orbits′ = primitivize.(orbits, cntr) # primitive basis
    i = 1
    while i ≤ length(orbits)
        wp′ = parent(orbits′[i])
        wp′_r = RVec{D}(reduce_translation_to_unitrange(constant(wp′)), # coords. in [0,1)ᴰ
                        free(wp′))
        if isapproxin(wp′_r, (@view orbits′[1:i-1]), nothing, #=modw=# true)
            # `wp′` is equivalent to another position already in the orbit; delete it!
            deleteat!(orbits′, i)
            deleteat!(orbits, i)
            deleteat!(gαs, i)
            continue
        end
        # the position `wp′` is not equivalent to any position in `orbit[1:i-1]`: keep it!
        # we add `wp′_reduced` rather than `wp′`, because we want the coordinates to be
        # in [0,1)ᴰ for the primitive unit cell
        if !isapprox(wp′_r, wp′, nothing, #=modw=# false)
            orbits′[i] = wp′_r
            orbits[i] = conventionalize(wp′_r, cntr)

            # we also need to update the coset operation, cf. the additional translation
            # of `wp′_r` relative to `wp′`
            Δ′ = constant(wp′_r) - constant(wp′)
            g′ = primitivize(gαs[i], cntr, #=modw=# false)     # original coset operation (primitive basis)
            g′_r = compose(typeof(g′)(Δ′), g′, #=modτ=# false) # "reduced" coset operation (primitive basis)
            # NB: ↑ For `g′ = (W′|w′|U′)` (`U′` absent in spinless case), this simply gives
            #     `g′_r = (W′|w′+Δ′|U′)`; i.e. `g′` translated by `Δ′`; `compose` is used
            #     only because it lets us pass on `U′` in the spinful case (`DSymOperation`)
            #     without passing anything in the spinless case (`SymOperation`)
            g_r = conventionalize(g′_r, cntr, #=modw=# false)  # "reduced" coset operation (conventional basis)
            gαs[i] = g_r
        end
        i += 1 # process next element
    end

    wp = parent(position(siteg))
    wp_r = first(orbits)
    if wp != wp_r
        # we always assume that the first position in the orbits is a reference point for
        # the cosets, i.e., is the canonical Wyckoff position `wp`, such that we can
        # obtain the orbits by acting with the cosets on `wp` - but if the new orbits'
        # first element, i.e., `wp_r`, is a different Wyckoff position than `wp`, we can't
        # go ahead directly - the cosets are still relative to `wp`, but we'd like them to
        # be relative to `wp_r`; so, we need to adjust the cosets accordingly. To adjust,
        # we exploit that `wp` and `wp_r` differ by exactly `gαs[1]` (the first coset
        # operation) in the sense that `wp_r = compose(gαs[1], wp, false)`; to adjust, we
        # just apply the inverse of `gαs[1]` to all coset operations from the right:
        g = gαs[1]
        g⁻¹ = inv(g)
        for i in eachindex(gαs)
            gαs[i] = compose(gαs[i], g⁻¹, #=modw=# false) # adjust to new origin
        end
        # NB: this looks assymmetrical relative to the transformation below, but that's only
        #     because we've already applied the left-hand side transformation

        # the change of reference point will also affect the site group (i.e., the site 
        # group is now relative to `wp_r` rather than `wp`); we must adjust the operations
        # Deriviation of transformation below: 
        # (1) we have: wp_r = g*wp ⇔ wp = g⁻¹*wp_r
        # (2) existing site group: s*wp = wp        for s in sitegroup(wp)
        # (3) wanted site group:   s_r*wp_r = wp_r  for s_r in sitegroup(wp_r)
        # (4) combine: s*g⁻¹*wp_r = g⁻¹*wp_r ⇔ (g*s*g⁻¹)*wp_r = wp_r ⇒ s_r = g*s*g⁻¹
        ops_r = compose.(compose.(Ref(g), operations(siteg), false), Ref(g⁻¹), false)
        wp_r = WyckoffPosition{D}(siteg.wp.mult, siteg.wp.letter, wp_r)
        siteg_r = typeof(siteg)(num(siteg), wp_r, ops_r, gαs)
    else
        siteg_r = typeof(siteg)(num(siteg), position(siteg), operations(siteg), gαs)
    end

    return siteg_r, orbits
end

function reduce_orbits_and_cosets(siteir::AbstractSiteIrrep)
    siteg_r, _ = reduce_orbits_and_cosets(group(siteir))
    return _rebuild_irrep_with_modified_group(siteir, siteg_r)
end

# ---------------------------------------------------------------------------------------- #
# Bandrep related functions: induction/subduction

"""
    induce_bandrep(siteir::AbstractSiteIrrep, h::AbstractOperation, kv::KVec)

Return the band representation induced by the provided site symmetry irrep evaluated at `kv`
and for an operation `h` (a `DSymOperation` for a double-valued `DSiteIrrep`).

It is assumed that `group(siteir)` is not centering-reduced: i.e., a centering-reduction
attempt is always made; if the group cosets are already reduced in the sense of
[`reduce_orbits_and_cosets`](@ref), this makes no difference.
"""
function induce_bandrep(
    siteir::AbstractSiteIrrep{D},
    h::AbstractOperation{D},
    kv::KVec{D},
) where D
    
    siteg, orbits = reduce_orbits_and_cosets(group(siteir))
    return _induce_bandrep(characters(siteir), h, kv, siteg, orbits) # FIXME: `orbits` isa not right type
end

function _induce_bandrep(
        χs::Vector{ComplexF64},  # characters of site irrep
        h::AbstractOperation{D},
        kv::KVec{D},
        siteg::AbstractSiteGroup{D},               # (centering-reduced) site group
        orbits::Vector{RVec{D}}, # (centering-reduced) orbit of `siteg``
    ) where D
    kv′ = constant(h*kv) # <-- TODO: Why only constant part?
    gαs = cosets(siteg) # (centering-reduced) cosets of the site group
    # sum over all the (non-centering-equivalent) wyckoff positions/cosets in the orbit 
    χᴳₖ = zero(ComplexF64)
    for (wpα′, gα′) in zip(orbits, gαs)
        wpα′ = parent(wpα′)
        tα′α′ = constant(h*wpα′ - wpα′) # TODO: <-- explain why we only need constant part here?
        opᵗ   = typeof(h)(-tα′α′)

        gα′⁻¹     = inv(gα′)
        gα′⁻¹ggα′ = compose(gα′⁻¹, compose(opᵗ, compose(h, gα′, false), false ), false)

        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), siteg)
        if site_symmetry_index !== nothing
            χᴳₖ += cispi(2*dot(kv′, tα′α′)) * χs[site_symmetry_index]
            # NB: The sign in this `cis(...)` above is different from in Elcoro's. 
            #     I think this is consistent with our overall sign convention (see #12),
            #     however, and flipping the sign causes problems for the calculation of some
            #     subductions to `LGIrrep`s, which would be consistent with this. I.e.,
            #     I'm fairly certain this is consistent and correct given our phase
            #     conventions for `LGIrrep`s.
        end
    end
    return χᴳₖ
end

function subduce_onto_lgirreps(
        siteir_χs::AbstractVector{<:Number},
        siteg::AbstractSiteGroup{D},
        lgirs::AbstractVector{<:AbstractLGIrrep{D}}
    ) where D
    lg = group(first(lgirs))
    kv = position(lg)

    # characters of induced site representation and little group irreps (we use
    # `_induce_bandrep` below instead of `induce_bandrep` to avoid repeated calculation of
    # the centering-reduction of orbits and cosets of the site irrep/site group)
    orbits = parent.(orbit(siteg))
    site_χs  = _induce_bandrep.(Ref(siteir_χs), lg, Ref(kv), Ref(siteg), Ref(orbits))
    lgirs_χm = matrix(characters(lgirs))

    # little group irrep multiplicities, after subduction
    m  = lgirs_χm\site_χs
    m′ = round.(Int, real.(m)) # truncate to integers
    isapprox(m, m′, atol=DEFAULT_ATOL) || error(DomainError(m, "failed to convert to integers"))
    return m′
end
function subduce_onto_lgirreps(
    siteir::AbstractSiteIrrep{D}, lgirs::AbstractVector{<:AbstractLGIrrep{D}}
) where D
    return subduce_onto_lgirreps(characters(siteir), group(siteir), lgirs)
end

# ---------------------------------------------------------------------------------------- #

function calc_bandrep(
        siteir :: AbstractSiteIrrep{D},
        lgirsv :: AbstractVector{<:AbstractVector{<:AbstractLGIrrep{D}}},
        timereversal :: Bool
    ) where D

    # take the input site symmetry irrep, and reduce its little group cosets such that the
    # associated orbit lies in [0,1)ᴰ of the primitive unit cell
    siteir = reduce_orbits_and_cosets(siteir)
    siteg = group(siteir)
    siteir_χs = characters(siteir)
    multsv = [subduce_onto_lgirreps(siteir_χs, siteg, lgirs) for lgirs in lgirsv]
    
    occupation = sum(zip(first(multsv), first(lgirsv)); init=0) do (m, lgir)
        m * irdim(lgir)
    end
    n = SymmetryVector(lgirsv, multsv, occupation)

    return BandRep(siteir, n, timereversal)
end
function calc_bandrep(
        siteir :: AbstractSiteIrrep{D};
        timereversal :: Bool=true, 
        allpaths :: Bool=false
    ) where D
    lgirsd = lgirreps(num(siteir), Val(D); spinful=Val(isspinful(siteir)))
    allpaths || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && realify!(lgirsd)
    lgirsv = _collect_lgirsd_sorted(lgirsd)
    return calc_bandrep(siteir, lgirsv, timereversal)
end

# ---------------------------------------------------------------------------------------- #
"""
    bandreps(
        sgnum::Integer,
        ::Val{D}=Val(3);
        spinful::Union{Bool, Val{true}, Val{false}}=Val(false),
        timereversal::Bool=true,
        allpaths::Bool=false,
        explicitly_real::Bool=timereversal
    ) --> Collection{BandRep{D}}

    bandreps( # type-unstable convenience accessor
        sgnum::Integer,
        D::Integer;
        kws...
    ) --> Collection{BandRep{D}}

Compute the band representations of space group `sgnum` in dimension `D`.

If `spinful` is `Val(true)` (or `true`), the spinful band representations are computed
instead, induced from the double-valued site symmetry irreps (see [`siteirreps`](@ref)) and
subduced onto the double-valued little group irreps (currently available in 3D only). As
for `D`, the `Val` spelling keeps the return type inferrable and the `Bool` spelling does
not.

## Keyword arguments
- `timereversal` (default, `true`): whether the irreps used to induce the band
  representations are assumed to be time-reversal invariant (i.e., are coreps, see 
  [`realify`](@ref)).
- `allpaths` (default, `false`): whether the band representations are projected to all
  distinct **k**-points returned by `lgirreps` (`allpaths = false`), including high-symmetry
  **k**-lines and -plane, or only to the maximal **k**-points (`allpaths = true`), i.e.,
  just to high-symmetry points.
- `explicitly_real` (default, `timereversal`): whether, if `timereversal = true`, to
  ensure that the site symmetry irreps accompanying the band representations are chosen
  in the canonical form associated with time reversal (see [`physical_realify`](@ref)),
  i.e., explicitly real for spinless irreps and `J*conj(D)*J' = D` for spinful
  ones. This can be helpful for subsequent analysis involving the action of time-reversal
  symmetry.
- `include_nonmaximal` (default, `false`): whether to include band representations induced
  from site symmetry irreps of non-maximal Wyckoff positions. Passing as `true` will include
  band representations induced from all Wyckoff positions, regardless of maximality.

## Notes
All band representations associated with maximal Wyckoff positions are returned, 
irregardless of whether they are elementary (i.e., no regard is made to whether the band
representation is "composite"). As such, the returned band representations generally are
a superset of the set of elementary band representations (and so contain all elementary
band representations).

## Implementation
The implementation is based on Cano, Bradlyn, Wang, Elcoro, et al., [Phys. Rev. B **97**,
035139 (2018)](https://doi.org/10.1103/PhysRevB.97.035139), Sections II.C-D.
"""
function bandreps(
        sgnum::Integer,
        Dᵛ::Val{D} = Val(3);
        spinful = Val(false),
        timereversal::Bool = true,
        allpaths::Bool = false,
        explicitly_real::Bool = timereversal,
        include_nonmaximal::Bool = false,
    ) where D

    if explicitly_real && !timereversal
        error("`explicitly_real = true` is only meaningful for `timereversal = true`")
    end

    # get all the little group irreps that we want to subduce onto
    lgirsd = lgirreps(sgnum, Dᵛ; spinful)
    allpaths || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && realify!(lgirsd)
    lgirsv = _collect_lgirsd_sorted(lgirsd)

    # get the bandreps induced by every maximal site symmetry irrep
    sg = spacegroup(sgnum, Dᵛ; spinful)
    sitegs = sitegroups(sg)
    if !include_nonmaximal
        sitegs = findmaximal(sitegs)
    end
    brs = _bandrep_type(_spinfulval(spinful), Dᵛ)[]
    for siteg in sitegs
        siteirs = siteirreps(siteg; mulliken=true)
        if timereversal
            siteirs = realify(siteirs)
            explicitly_real && (siteirs = physical_realify(siteirs))
        end
        append!(brs, calc_bandrep.(siteirs, Ref(lgirsv), Ref(timereversal)))
    end

    return Collection(brs)
end
bandreps(sgnum::Integer, D::Integer; kws...) = bandreps(sgnum, Val(D); kws...)

# the band representation type induced by spinless or by spinful site symmetry irreps; keyed
# on `Val`s so that the type is fixed by dispatch, and so is propagated even if the
# dimension is not a compile-time constant
_bandrep_type(#=Val{S}=#::Val{false}, ::Val{D}) where D = BandRep{D, LGIrrep{D}, SiteIrrep{D}}
_bandrep_type(#=Val{S}=#::Val{true}, ::Val{D}) where D = BandRep{D, DLGIrrep{D}, DSiteIrrep{D}}

# ---------------------------------------------------------------------------------------- #

# performance optimization
function Base.stack(brs::Collection{<:BandRep})
    B = Matrix{Int}(undef, length(first(brs)), length(brs))
    @inbounds for (j, br) in enumerate(brs)
        i = 1
        for mults in multiplicities(br)
            n = length(mults)
            i′ = i + n - 1
            B[i:i′, j] = mults
            i = i′ + 1
        end
        B[i, j] = occupation(br)
    end
    return B
end
