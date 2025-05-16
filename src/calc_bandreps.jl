using LinearAlgebra: dot, \

# The implementation here follows Cano et al., Phys. Rev. B 97, 035139 (2018)
# (https://doi.org/10.1103/PhysRevB.97.035139), specifically, Sections II.C-D

# ---------------------------------------------------------------------------------------- #

"""
    reduce_dict_of_vectors([f=identity,] d::Dict{_, <:Vector})  -->  Vector{T}

Return the concatenated vector of all of vectors in `d` under the element-wise application
of `f`. Effectively flattens a `Dict` of `Vector`s to a single `Vector`.

Note that application of `f` to vector-elements of `d` must return a stable type `T`.
"""
function reduce_dict_of_vectors(f::F, d::Dict{<:Any, <:AbstractVector}) where F
    # get element type of application of `f` to elements of vectors in `d`; assumed fixed
    eltype = typeof(f(first(first(values(d)))))
    # get total number of elements across vectors in `d`
    N = sum(((_, v),) -> length(v), d)
    # preallocate output vector
    w     = Vector{eltype}(undef, N)
    start = 1
    for v in values(d)
        stop = start + length(v) - 1
        @inbounds for (i, j) in enumerate(start:stop)
            w[j] = f(v[i])
        end
        start = stop + 1
    end
    return w
end
reduce_dict_of_vectors(d::Dict{<:Any, <:Vector}) = reduce_dict_of_vectors(identity, d)

"""
    reduce_orbits_and_cosets(siteg::SiteGroup{D}))

For an input site group, provided in conventional coordinates, reduce its cosets such
that the resulting orbit only contains Wyckoff positions that are not equivalent when
viewed in the primitive basis (as determined by the centering type `cntr`). Additionally,
the associated orbit will only contain positions whose _primitive_ coefficients lie in
[0,1). I.e., the associated orbit lies in the canonical [0,1)ᴰ primitive unit cell. The
coset operations are adjusted accordingly.

The reduced site group is returned in a conventional basis, along with the reduced orbits,
also in a conventional basis.
"""
function reduce_orbits_and_cosets(
        siteg::SiteGroup{D}
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
            g′ = primitivize(gαs[i], cntr, #=modw=# false) # original coset operation (primitive basis)
            g′_r = SymOperation(g′.rotation, g′.translation + Δ′) # "reduced" coset operation (primitive basis)
            g_r = conventionalize(g′_r, cntr, #=modw=# false) # "reduced" coset operation (conventional basis)
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
        siteg_r = SiteGroup{D}(num(siteg), wp_r, ops_r, gαs)
    else
        siteg_r = SiteGroup{D}(num(siteg), position(siteg), operations(siteg), gαs)
    end

    return siteg_r, orbits
end

function reduce_orbits_and_cosets(siteir::SiteIrrep{D}) where D
    siteg_r, _ = reduce_orbits_and_cosets(group(siteir))
    siteir_r = SiteIrrep{D}(siteir.cdml, siteg_r, siteir.matrices, siteir.reality,
                           siteir.iscorep, siteir.pglabel)
    return siteir_r
end

# ---------------------------------------------------------------------------------------- #
# Bandrep related functions: induction/subduction

"""
    induce_bandrep(siteir::SiteIrrep, h::SymOperation, kv::KVec)

Return the band representation induced by the provided `SiteIrrep` evaluated at `kv` and for
a `SymOperation` `h`.

It is assumed that `group(siteir)` is not centering-reduced: i.e., a centering-reduction
attempt is always made; if the group cosets are already reduced in the sense of
[`reduce_orbits_and_cosets`](@ref), this makes no difference.
"""
function induce_bandrep(
    siteir::SiteIrrep{D},
    h::SymOperation{D},
    kv::KVec{D},
    ) where D   
    
    siteg, orbits = reduce_orbits_and_cosets(group(siteir))
    return _induce_bandrep(characters(siteir), h, kv, siteg, orbits) # FIXME: `orbits` isa not right type
end

function _induce_bandrep(
        χs::Vector{ComplexF64},  # characters of site irrep
        h::SymOperation{D},
        kv::KVec{D},
        siteg::SiteGroup{D},     # (centering-reduced) site group
        orbits::Vector{RVec{D}}, # (centering-reduced) orbit of `siteg``
    ) where D
    kv′ = constant(h*kv) # <-- TODO: Why only constant part?
    gαs = cosets(siteg) # (centering-reduced) cosets of the site group
    # sum over all the (non-centering-equivalent) wyckoff positions/cosets in the orbit 
    χᴳₖ = zero(ComplexF64)
    for (wpα′, gα′) in zip(orbits, gαs)
        wpα′ = parent(wpα′)
        tα′α′ = constant(h*wpα′ - wpα′) # TODO: <-- explain why we only need constant part here?
        opᵗ   = SymOperation(-tα′α′)

        gα′⁻¹     = inv(gα′)
        gα′⁻¹ggα′ = compose(gα′⁻¹, compose(opᵗ, compose(h, gα′, false), false ), false)

        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), siteg)
        if site_symmetry_index !== nothing
            χᴳₖ += cis(2π*dot(kv′, tα′α′)) * χs[site_symmetry_index]
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
        siteg::SiteGroup{D},
        lgirs::AbstractVector{LGIrrep{D}}
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
    siteir::SiteIrrep{D}, lgirs::AbstractVector{LGIrrep{D}}
) where D
    return subduce_onto_lgirreps(characters(siteir), group(siteir), lgirs)
end

# ---------------------------------------------------------------------------------------- #

function calc_bandrep(
        siteir :: SiteIrrep{D}, 
        lgirsv :: AbstractVector{<:AbstractVector{LGIrrep{D}}},
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
    
    spinful = false # NB: default; Crystalline currently doesn't have spinful irreps

    return NewBandRep(siteir, n, timereversal, spinful)
end
function calc_bandrep(
        siteir :: SiteIrrep{D}; 
        timereversal :: Bool=true, 
        allpaths :: Bool=false
    ) where D
    lgirsd = lgirreps(num(siteir), Val(D))
    allpaths || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && realify!(lgirsd)
    lgirsv = [lgirs for lgirs in values(lgirsd)]
    return calc_bandrep(siteir, lgirsv, timereversal)
end

# ---------------------------------------------------------------------------------------- #
"""
    calc_bandreps(sgnum::Integer, Dᵛ::Val{D}=Val(3);
                  timereversal::Bool=true,
                  allpaths::Bool=false,
                  explicitly_real::Bool=timereversal)

Compute the band representations of space group `sgnum` in dimension `D`, returning a
`BandRepSet`.

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
  to be explicitly real (or "physically" real; see [`physical_realify`](@ref)). This
  is helpful for subsequent analysis of the action of time-reversal symmetry.

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
function calc_bandreps(
        sgnum::Integer,
        Dᵛ::Val{D} = Val(3);
        timereversal::Bool = true,
        allpaths::Bool = false,
        explicitly_real::Bool = timereversal
    ) where D

    if explicitly_real && !timereversal
        error("`explicitly_real = true` is only meaningful for `timereversal = true`")
    end

    # get all the little group irreps that we want to subduce onto
    lgirsd = lgirreps(sgnum, Val(D))
    allpaths || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && realify!(lgirsd)
    lgirsv = [lgirs for lgirs in values(lgirsd)]

    # get the bandreps induced by every maximal site symmetry irrep
    sg     = spacegroup(sgnum, Dᵛ)
    sitegs = findmaximal(sitegroups(sg))
    brs    = NewBandRep{D}[]
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

# ---------------------------------------------------------------------------------------- #

# performance optimization
function Base.stack(brs::Collection{NewBandRep{D}}) where D
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
