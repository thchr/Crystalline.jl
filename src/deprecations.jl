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

@deprecate kvec(lg::LittleGroup)     position
@deprecate kvec(lgir::LGIrrep)       position
@deprecate wyck(g::SiteGroup)        position
@deprecate wyck(wp::WyckoffPosition) position
@deprecate wyck(BR::BandRep)         position