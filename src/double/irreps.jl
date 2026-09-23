# --- Double group irreps ---

"""
$(TYPEDEF)$(TYPEDFIELDS)

A double-valued (spinful) irrep of a little group, i.e., an irrep of a double little group
`g` in which ``\\bar{E}`` is represented by `-𝟙`.

The fields mirror those of [`LGIrrep`](@ref), with one matrix and one translation for each
of the `2|G|` operations of `g`.
"""
@struct_hash_equal struct DLGIrrep{D} <: AbstractLGIrrep{D}
    cdml::String
    g::DLittleGroup{D}
    matrices::Vector{Matrix{ComplexF64}}
    translations::Vector{SVector{D, Float64}}
    reality::Reality
    iscorep::Bool
    # constructor with automatic conversion & defaults
    function DLGIrrep{D}(
        cdml::AbstractString,
        g::DLittleGroup{D},
        matrices::AbstractVector{<:AbstractMatrix{<:Number}},
        translations::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing},
        reality::Reality,
        iscorep::Bool = false
    ) where {D}
        # `nothing` is a sentinel value for all-zero translations
        translations = if translations === nothing
            [zeros(SVector{D, Float64}) for _=OneTo(order(g))]
        else
            convert(Vector{SVector{D, Float64}}, translations)
        end
        if !(length(matrices) == length(translations) == order(g))
            error("length(matrices) (=$(length(matrices))), length(translations) \
                       (=$(length(translations))), & length(operations(g)) \
                       (=$(length(operations(g)))) must all be equal")
        end
        matrices = convert(Vector{Matrix{ComplexF64}}, matrices)
        return new{D}(String(cdml), g, matrices, translations, reality, iscorep)
    end
end
DLGIrrep(cdml::String, g::DLittleGroup{D}, args...) where D = DLGIrrep{D}(cdml, g, args...)

"""
$(TYPEDEF)$(TYPEDFIELDS)

A double-valued (spinful) irrep of a point group, i.e., an irrep of a double point group `g`
in which ``\\bar{E}`` is represented by `-𝟙`.

The fields mirror those of [`PGIrrep`](@ref), with one matrix for each of the `2|G|`
operations of `g`.
"""
@struct_hash_equal struct DPGIrrep{D} <: AbstractPGIrrep{D}
    cdml::String
    g::DPointGroup{D}
    matrices::Vector{Matrix{ComplexF64}}
    reality::Reality
    iscorep::Bool
end
function DPGIrrep{D}(
    cdml::String,
    pg::DPointGroup{D},
    matrices::Vector{Matrix{ComplexF64}},
    reality::Reality
) where D
    return DPGIrrep{D}(cdml, pg, matrices, reality, false)
end

"""
$(TYPEDEF)$(TYPEDFIELDS)

A double-valued (spinful) irrep of a site symmetry group, i.e., an irrep of a double site
symmetry group `g` in which ``\\bar{E}`` is represented by `-𝟙`.

The fields mirror those of [`SiteIrrep`](@ref), with one matrix for each of the `2|G|`
operations of `g`.
"""
@struct_hash_equal struct DSiteIrrep{D} <: AbstractSiteIrrep{D}
    cdml     :: String
    g        :: DSiteGroup{D}
    matrices :: Vector{Matrix{ComplexF64}}
    reality  :: Reality
    iscorep  :: Bool
    pglabel  :: String # label of point group that is isomorphic to the site group `g`
end

"""
    isspinful(x) -> Bool

Return whether `x` (an irrep, an irrep type, a character table, a symmetry vector, or a band
representation) is spinful, i.e., double-valued, as appropriate for half-integer angular
momentum. Otherwise, `x` is spinless, i.e., single-valued, as appropriate for integer
angular momentum.
"""
isspinful(ir::AbstractIrrep) = isspinful(typeof(ir))
isspinful(::Type{<:AbstractIrrep}) = false
isspinful(::Type{<:Union{DLGIrrep, DPGIrrep, DSiteIrrep}}) = true
# a character table is spinful if it is a table over a double group
isspinful(ct::AbstractCharacterTable) = eltype(operations(ct)) <: DSymOperation

# --- Loading (see `lgirreps` in /src/littlegroup_irreps.jl and `pgirreps` in
#     /src/pointgroup.jl) ---
# The data files store a double-valued irrep only on the operations of the ordinary little
# or point group; the barred operations follow them in the double group (see
# `doubled_operations`) and are represented by `D(Ēg) = -D(g)`
_doubled_matrices(Ps) = vcat(Ps, [-P for P in Ps])
_doubled_translations(::Nothing) = nothing
_doubled_translations(τs) = vcat(τs, τs)
