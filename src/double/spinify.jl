"""
    spinify(ir::AbstractIrrep, irsˢ::Collection{<:AbstractIrrep}; kws...) --> Vector{<:AbstractIrrep}
    spinify(n::AbstractSymmetryVector, brsˢ::Collection{<:BandRep}) --> SymmetryVector
    spinify(n::AbstractSymmetryVector, lgirsvˢ::AbstractVector{<:Collection{<:DLGIrrep}}) --> SymmetryVector
    spinify(br::BandRep, brsˢ::Collection{<:BandRep}) --> CompositeBandRep
    spinify(cbr::CompositeBandRep, brsˢ::Collection{<:BandRep}) --> CompositeBandRep

Return the spinful counterpart of a spinless irrep `ir`, symmetry vector `n`, or (composite)
band representation `br` (`cbr`), in the absence of any coupling between spin and orbital
degrees of freedom: i.e., answer the question "if these spinless states were in fact doubly
degenerate but uncoupled spin-½ states, how would they transform under the double group?".

The answer is expressed in terms of the double-valued irreps `irsˢ` (or the double-valued
little group irreps of the spinful band representations `brsˢ`, or `lgirsvˢ`, ordered as in
`n`), which must be those of the double group of the same group as `ir`, `n`, or `br` (i.e.,
at the same **k**-point or Wyckoff position) and treat time reversal in the same way (i.e.,
either both or neither are [`realify`](@ref)'d).

## Background
Without coupling between spin and orbital degrees of freedom, spatial operations ``g ∈ G``
and spin rotations ``u ∈ \\mathrm{SU}(2)`` act _independently_: the symmetry group is
``G × \\mathrm{SU}(2)``, and a spinless irrep ``D`` of ``G``, once equipped with spin,
transforms as the irrep ``D ⊠ D^{(1/2)}`` of ``G × \\mathrm{SU}(2)``, i.e., as
``(g, u) ↦ D(g) ⊗ D^{(1/2)}(u)`` for any ``g ∈ G`` and any ``u ∈ \\mathrm{SU}(2)``, with
``D^{(1/2)}`` the spin-½ irrep of SU(2) (practically, the 2×2 matrix form of ``u``)[^1].
Spin–orbit coupling removes this independence, permitting only the operations that
transform orbital and spin degrees of freedom together: these are the operations
``\\bar{g}`` of the double group ``\\bar{G}``, each of which pairs a spatial operation
``g`` with one of the two spin rotations ``u`` associated with it (differing by sign;
unbarred vs. barred). Subduced onto ``\\bar{G}``, ``D ⊠ D^{(1/2)}`` becomes the generally
reducible representation ``\\bar{g} ↦ D(g) ⊗ U(\\bar{g})``, where
``U(\\bar{g}) = D^{(1/2)}(u)`` is the 2×2 matrix of that particular ``u`` (see
[`SU2`](@ref)). `spinify` decomposes it into the double-valued irreps ``\\bar{D}_j`` of
``\\bar{G}`` (provided by `irsˢ`, `brsˢ`, or `lgirsvˢ`), with multiplicities ``m_j``:

```math
(D ⊠ D^{(1/2)}) ↓ \\bar{G} ≅ \\bigoplus_j m_j \\bar{D}_j.
```

Since ``\\mathrm{tr}(D(g) ⊗ U(\\bar{g})) = \\mathrm{tr} D(g)\\,\\mathrm{tr} U(\\bar{g})``,
the multiplicities follow from the characters of ``D`` and ``\\bar{D}_j`` alone.

The decomposition generally involves more than a single double-valued irrep.
Accordingly, the irrep method returns the featured irreps ``\\bar{D}_j``, each repeated
``m_j`` times. Likewise for band representations: a spinified band representation
is a sum of spinful band representations induced from the _same_ Wyckoff position — but
not necessarily a single one; it is therefore returned as a [`CompositeBandRep`](@ref) over
`brsˢ`.
The occupation (or dimension) is always doubled.

[^1]: ``⊠`` denotes the _outer_ tensor product: for representations ``D`` of ``G``
      and ``E`` of ``H``, ``(D ⊠ E)(g, h) = D(g) ⊗ E(h)`` is a representation of
      ``G × H ∋ (g, h)``, whose two factors act on separate spaces.

## Keyword arguments
For the irrep method, keyword arguments `kws` are passed to [`find_representation`](@ref).

## Examples

The three-dimensional irrep T₁ᵤ of p-orbitals in point group m-3m decomposes into two
double-valued irreps, of dimension 2 and 4 (``j = 1/2`` and ``j = 3/2``):

```jldoctest
julia> pgirs  = pgirreps("m-3m", Val(3); mulliken=true);

julia> pgirsˢ = pgirreps("m-3m", Val(3); spinful=Val(true), mulliken=true);

julia> T₁ᵤ = pgirs[8]; label(T₁ᵤ)
"T₁ᵤ"

julia> spinify(T₁ᵤ, pgirsˢ)
2-element Vector{DPGIrrep{3}}:
 E₁ᵤˢ
 Fᵤˢ
```

For band representations, e.g. the p_z orbitals of graphene (plane group p6mm):
```jldoctest
julia> brs = bandreps(17, Val(2));

julia> brsˢ = bandreps(17, Val(2); spinful=Val(true));

julia> br = brs[5]
13-irrep BandRep{2} (spinless):
 (2b|A₁): [Γ₁+Γ₄, K₃, M₁+M₄] (2 bands)

julia> spinify(br, brsˢ)
7-irrep CompositeBandRep{2} (spinful):
 (2b|E₁ˢ) (4 bands)

julia> spinify(SymmetryVector(br), brsˢ)
7-irrep SymmetryVector{2} (spinful):
 [Γ₈ˢ+Γ₉ˢ, K₄ˢ+K₅ˢ+K₆ˢ, 2M₅ˢ] (4 bands)
```
"""
function spinify end

# ---------------------------------------------------------------------------------------- #
# internal implementation methods

# return the multiplicities of a single spinified irrep `ir` relative to a collection of
# spinful irreps `irsˢ`
function _spinify(
    ir :: AbstractIrrep{D},
    irsˢ :: Collection{<:AbstractIrrep{D}};
    representation_kws...
) where D
    # check that `ir` is spinless & `irsˢ` is spinful
    !isspinful(ir) || error("input `ir` must be a spinless irrep")
    isspinful(first(irsˢ)) || error("input `irsˢ` must be a collection of spinful irreps")
    
    # check that `ir` & `irsˢ` are of the same abstract kind (*LG*, *Site*, *PG*)
    if !((ir isa AbstractLGIrrep   && eltype(irsˢ) <: AbstractLGIrrep) ||
         (ir isa AbstractPGIrrep   && eltype(irsˢ) <: AbstractPGIrrep) ||
         (ir isa AbstractSiteIrrep && eltype(irsˢ) <: AbstractSiteIrrep))
        error("input `ir` and `irsˢ` must represent related kinds of irreps (e.g. point, site, or little group irreps)")
    end

    # check that `ir` and `irsˢ` refer to the same underlying group 
    g, gˢ = group(ir), group(irsˢ)
    num(g) == num(gˢ) || error("the input `ir` and `irsˢ` refer to different underlying groups")

    # check that `ir` & `irsˢ` are at the same k-point or Wyckoff position, if little group
    # or site irreps (comparing labels suffices, since the spinless & spinful irreps share
    # their k-point & Wyckoff position coordinates; point groups have empty labels)
    poslab, poslabˢ = positionlabel(g), positionlabel(gˢ)
    poslab ≠ poslabˢ && error(lazy"the input `ir` and `irsˢ` refer to different positions: $poslab vs. $poslabˢ")

    χs = characters(ir) # spinless characters of `ir` across operations of `g`
    χsˢ = _spinified_characters(χs, g, gˢ) # spinified characters across operations of `gˢ`

    ms = find_representation(χsˢ, irsˢ; representation_kws...)
    isnothing(ms) && error(lazy"failed to decompose the spinified irrep $(label(ir)) into the irreps `irsˢ`")
    
    return something(ms)
end

# characters of the (reducible) spinified representation D(g) ⊗ U(ḡ), over the operations
# ḡ of the double group `gˢ`, with g the spatial part of ḡ and U(ḡ) its SU(2) matrix: since
# tr(X⊗Y) = tr(X)tr(Y), these are simply products of the characters `χs` of D (over the
# operations of `g`) and tr U(ḡ)
function _spinified_characters(
    χs::AbstractVector{<:Number}, g::AbstractGroup{D}, gˢ::AbstractGroup{D}
) where D
    length(gˢ) == 2length(g) || error("`gˢ` must be the double group of `g`")
    χˢ = zeros(ComplexF64, length(gˢ))
    for (i, opˢ) in enumerate(gˢ)
        idx_in_g = findfirst(g) do op′
            isapprox(op′, SymOperation(opˢ), #=cntr=# nothing, #=modw=# false)
        end
        isnothing(idx_in_g) && error(lazy"failed to find matching spinless operation for $opˢ")
        χˢ[i] = χs[something(idx_in_g)] * tr(SU2(opˢ))
    end
    return χˢ
end

function _spinified_subduction_count_from_characters(
    g::Union{AbstractPointGroup{D}, AbstractSiteGroup{D}},
    χs::AbstractVector{<:Number},
    irˢ::Union{AbstractPGIrrep{D}, AbstractSiteIrrep{D}}
) where D
    gˢ = group(irˢ)
    χsˢ = _spinified_characters(χs, g, gˢ)
    s = dot(characters(irˢ), χsˢ) # ∑ χˢ^* χ: same "magic formula" as in `subduction_count`
    return _checked_subduction_sum_to_count(s, irˢ)
end

# ---------------------------------------------------------------------------------------- #
# :: AbstractIrrep{D} method

function spinify(
    ir :: AbstractIrrep{D},
    irsˢ :: Collection{<:AbstractIrrep{D}};
    representation_kws...
) where D
    ms = _spinify(ir, irsˢ; representation_kws...)
    return [irsˢ[i] for (i, m) in enumerate(ms) for _ in 1:m]
end

# ---------------------------------------------------------------------------------------- #
# :: AbstractSymmetryVector methods

function spinify(
    n::AbstractSymmetryVector{D, LGIrrep{D}},
    lgirsvˢ::AbstractVector{Collection{DLGIrrep{D}}},
) where D
    lgirsv = irreps(n)
    multsv = multiplicities(n)
    multsvˢ = map(zip(multsv, lgirsv, lgirsvˢ)) do (mults, lgirs, lgirsˢ)
        multsˢ = zeros(Int, length(lgirsˢ))
        for (m, lgir) in zip(mults, lgirs)
            iszero(m) && continue
            multsˢ += m * _spinify(lgir, lgirsˢ)
        end
        multsˢ
    end
    return SymmetryVector(lgirsvˢ, multsvˢ, 2*occupation(n))
end

function spinify(
    n::AbstractSymmetryVector{D, LGIrrep{D}},
    brsˢ::Collection{BandRep{D, DLGIrrep{D}, DSiteIrrep{D}}}
) where D
    lgirsvˢ = irreps(brsˢ)
    return spinify(n, lgirsvˢ)
end


# ---------------------------------------------------------------------------------------- #
# :: BandRep{D} & CompositeBandRep{D} methods

function spinify(
    br::BandRep{D, LGIrrep{D}, SiteIrrep{D}},
    brsˢ::Collection{BandRep{D, DLGIrrep{D}, DSiteIrrep{D}}}
) where D
    coefsˢ = zeros(Rational{Int}, length(brsˢ))
    _spinify!(coefsˢ, br, brsˢ)
    return CompositeBandRep(coefsˢ, brsˢ)
end

function spinify(
    cbr::CompositeBandRep{D, LGIrrep{D}, SiteIrrep{D}},
    brsˢ::Collection{BandRep{D, DLGIrrep{D}, DSiteIrrep{D}}}
) where D
    # spinifying is linear, so we can spinify each band representation separately
    coefsˢ = zeros(Rational{Int}, length(brsˢ))
    for (c, br) in zip(cbr.coefs, cbr.brs)
        iszero(c) && continue
        _spinify!(coefsˢ, br, brsˢ, c)
    end
    return CompositeBandRep(coefsˢ, brsˢ)
end

# add `c` times the spinification of `br`, in terms of `brsˢ`, to `coefsˢ`
function _spinify!(
    coefsˢ::AbstractVector{<:Real},
    br::BandRep{D, LGIrrep{D}, SiteIrrep{D}},
    brsˢ::Collection{BandRep{D, DLGIrrep{D}, DSiteIrrep{D}}},
    c :: Union{Nothing, <: Real} = nothing
) where D
    siteir = br.siteir
    siteg = group(siteir)
    χs = characters(siteir)

    μˢ = 0 # occupation of the spinification of `br`
    for (i, brˢ) in enumerate(brsˢ)
        position(br) == position(brˢ) || continue # must be at same Wyckoff position
        n = _spinified_subduction_count_from_characters(siteg, χs, brˢ.siteir)
        coefsˢ[i] += isnothing(c) ? n : c * n
        μˢ += n * occupation(brˢ)
    end
    if μˢ ≠ 2occupation(br)
        error(lazy"the spinified band representation has occupation $μˢ rather than the expected $(2occupation(br))")
    end
    return coefsˢ
end
