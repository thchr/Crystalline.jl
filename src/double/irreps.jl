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

# --- Time reversal ---
# For spinful systems, time reversal squares to -1, and `realify`'s pairing of irreps into
# co-representations does not apply
function realify(::AbstractVector{<:Union{DLGIrrep, DPGIrrep}}; kws...)
    error("co-representations of double-valued irreps are not yet implemented")
end

# --- Loading (see `lgirreps` in /src/littlegroup_irreps.jl and `pgirreps` in
#     /src/pointgroup.jl) ---
# The data file stores a double-valued irrep only on the operations of the ordinary little
# group; the barred operations follow them in the double little group (see
# `doubled_operations`) and are represented by `D(Ēg) = -D(g)`
function _lgirrep(cdml, lg::DLittleGroup{D}, P, τ, reality) where D
    return DLGIrrep{D}(cdml, lg, _doubled_matrices(P), _doubled_translations(τ), reality)
end
function _pgirrep(cdml, pg::DPointGroup{D}, P, reality) where D
    return DPGIrrep{D}(cdml, pg, _doubled_matrices(P), reality)
end
_doubled_matrices(Ps) = vcat(Ps, [-P for P in Ps])
_doubled_translations(::Nothing) = nothing
_doubled_translations(τs) = vcat(τs, τs)
