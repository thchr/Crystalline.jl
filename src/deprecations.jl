@deprecate get_littlegroups littlegroups
@deprecate get_lgirreps     lgirreps
@deprecate get_pgirreps     pgirreps
@deprecate get_wycks        wyckoffs
@deprecate get_siteirreps   siteirreps
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