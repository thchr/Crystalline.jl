@deprecate get_littlegroups littlegroups
@deprecate get_lgirreps     lgirreps
@deprecate get_pgirreps     pgirreps
@deprecate get_wycks        wyckoffs
@deprecate WyckPos          WyckoffPosition

@deprecate(
    CharacterTable(irs::AbstractVector{<:AbstractIrrep},
                   αβγ::Union{AbstractVector{<:Real}, Nothing}),
    characters(irs, αβγ)
)
@deprecate(
    CharacterTable(irs::AbstractVector{<:AbstractIrrep}),
    characters(irs)
)

@deprecate kvec(lg::LittleGroup)     position(lg)
@deprecate kvec(lgir::LGIrrep)       position(lgir)
@deprecate wyck(g::SiteGroup)        position(g)
@deprecate wyck(wp::WyckoffPosition) parent(wp)
@deprecate wyck(BR::BandRep)         position(BR)

@deprecate(
    kstar(g::AbstractVector{SymOperation{D}}, kv::KVec{D}, cntr::Char) where D,
    orbit(g, kv, cntr)
    )
@deprecate kstar(lgir::LGIrrep) orbit(lgir)

@deprecate(
    SiteGroup(sg::SpaceGroup{D}, wp::WyckoffPosition{D}) where D,
    sitegroup(sg, wp)
)

@deprecate(
    SiteGroup(sgnum::Integer, wp::WyckoffPosition{D}) where D,
    sitegroup(sgnum, wp)
)

Base.@deprecate_binding( # cannot use plain `@deprecate` for types
    IrrepCollection,
    Collection{T} where T<:AbstractIrrep
)

@deprecate matrix(brs::BandRepSet; kws...) stack(brs::BandRepSet)

@deprecate classification indicator_group_as_string
@deprecate nontrivial_factors indicator_group