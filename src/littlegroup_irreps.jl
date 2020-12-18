const JldOrNothing = Union{Nothing,JLD2.JLDFile}

# Little group operations loading
"""
    get_littegroups(sgnum::Integer, D::Union{Val{Int}, Integer}=Val(3)) 
                                                        -> Dict{String, LittleGroup{D}}

For given space group number `sgnum` and dimension `D`, return the associated little groups
(`LittleGroups{D}`s) at high-symmetry k-points, lines, and planes (see also
[`get_lgirreps`](@ref)).

Returns a `Dict` with little group **k**-point labels as keys and vectors of
`LittleGroup{D}`s as values.

## Notes
A conventional crystallographic setting is assumed (as in [`spacegroup`](@ref)).

Unlike `spacegroup`, "centering"-copies of symmetry operations are not included in the
returned `LittleGroup`s; as an example, space group 110 (body-centered, with centering
symbol 'I') has a centering translation `[1/2,1/2,1/2]` in the conventional setting:
the symmetry operations returned by `spacegroup` thus includes e.g. both `{1|0}` and 
`{1|½,½,½}` while the symmetry operations returned by `get_littlegroups` only include
`{1|0}` (and so on).

Currently, only `D = 3` is supported.

## References
The underlying data is sourced from the ISOTROPY dataset: see also [`get_lgirreps`](@ref).
"""
function get_littlegroups(sgnum::Integer, ::Val{D}=Val(3),
                          jldfile::JldOrNothing=nothing) where D
     D ∉ (1,2,3) && _throw_invaliddim(D)
     (D==2 && sgnum ∈ (4,7,8,12)) && error("Nonsymmorphic groups not yet supported in 2D")

    sgops_str, klabs, kstrs, opsidxs = if isnothing(jldfile)
        JLD2.jldopen(pathof_littlegroups_data(D), "r") do jldfile
             _load_littlegroups_data(sgnum, jldfile)
        end
    else
        _load_littlegroups_data(sgnum, jldfile)
    end

    sgops = SymOperation{D}.(sgops_str)
    lgs = Dict{String, LittleGroup{D}}()
    @inbounds for (klab, kstr, opsidx) in zip(klabs, kstrs, opsidxs)
        lgs[klab] = LittleGroup{D}(sgnum, KVec{D}(kstr), klab, sgops[opsidx])
    end
    return lgs
end

function get_all_littlegroups(::Val{D}=Val(3)) where D
    JLD2.jldopen(pathof_littlegroups_data(D),"r") do lgfile
        return [get_littlegroups(sgnum, Val(D), lgfile) for sgnum in OneTo(MAX_SGNUM[D])]
    end
end
# convenience functions without Val(D) usage; avoid internally
get_littlegroups(sgnum::Integer, D::Integer, jldfile::JldOrNothing=nothing) = get_littlegroups(sgnum, Val(D), jldfile)
get_all_littlegroups(D::Integer) = get_all_littlegroups(Val(D))


#------------------------------------------------------------------------------------------
# LGIrrep loading
"""
    get_lgirreps(sgnum::Integer, D::Union{Val{Int}, Integer}=Val(3))
                                                    -> Dict{String, Vector{LGIrrep{D}}}

For given space group number `sgnum` and dimension `D`, return the associated little group
(or "small") irreps (`LGIrrep{D}`s) at high-symmetry k-points, lines, and planes. 

Returns a `Dict` with little group **k**-point labels as keys and vectors of `LGIrrep{D}`s
as values.

## Notes
- Currently, only `D = 3` is supported.
- The returned irreps are complex in general. Real irreps (as needed in time-reversal
  invariant settings) can subsequently be obtained with the [`realify`](@ref) method.
- Returned irreps are spinless.

## References
The underlying data is sourced from the ISOTROPY ISO-IR dataset. If used in research, please
cite the original reference material associated with ISO-IR:

- Stokes, Hatch, & Campbell, [ISO-IR, ISOTROPY Software Suite](https://stokes.byu.edu/iso/irtables.php)
- Stokes, Campbell, & Cordes, [Acta Cryst. A. **69**, 388-395 (2013)](https://doi.org/10.1107/S0108767313007538).

The ISO-IR dataset is occasionally missing some **k**-points that lie outside the basic
domain but still resides in the representation domain (i.e. **k**-points with postscripted
'A', 'B', etc. labels, such as 'ZA'). In such cases, the missing irreps may instead have
been manually sourced from the Bilbao Crystallographic Database.
"""
function get_lgirreps(sgnum::Integer, Dᵛ::Val{D}=Val(3), lgs_jldfile::JldOrNothing=nothing,
                      irs_jldfile::JldOrNothing=nothing) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)
    (D==2 && sgnum ∈ (4,7,8,12)) && error("Nonsymmorphic groups not yet supported in 2D")
  
    lgs = get_littlegroups(sgnum, Dᵛ, lgs_jldfile)

    Ps_list, τs_list, types_list, cdmls_list = if isnothing(irs_jldfile)
        JLD2.jldopen(pathof_lgirreps_data(D), "r") do irs_jldfile
            _load_lgirreps_data(sgnum, irs_jldfile)
        end
    else
        _load_lgirreps_data(sgnum, irs_jldfile)
    end

    lgirsd = Dict{String, Vector{LGIrrep{D}}}()
    for (Ps, τs, types, cdmls) in zip(Ps_list, τs_list, types_list, cdmls_list)
        klab = klabel(first(cdmls))
        lg   = lgs[klab]
        lgirsd[klab] = [LGIrrep{D}(cdml, lg, P, τ, type) for (P, τ, type, cdml) in zip(Ps, τs, types, cdmls)]
    end
    
    return lgirsd
end
function get_lgirreps(sgnum::Integer, D::Integer, lgs_jldfile::JldOrNothing=nothing, 
                      irs_jldfile::JldOrNothing=nothing)
    get_lgirreps(sgnum, Val(D), lgs_jldfile, irs_jldfile)
end

function get_all_lgirreps(Dᵛ::Val{D}=Val(3)) where D
    JLD2.jldopen(pathof_littlegroups_data(D),"r") do lgfile
        JLD2.jldopen(pathof_lgirreps_data(D),"r") do irfile
            return [get_lgirreps(sgnum, Dᵛ, lgfile, irfile) for sgnum in OneTo(MAX_SGNUM[D])]
        end
    end
end
get_all_lgirreps(D::Integer) = get_all_lgirreps(Val(D))

# ===== utility functions (loads raw data from the harddisk) =====
const DATA_LITTLEGROUPS_FILENAME = "littlegroups_data.jld2"
const DATA_LGIRREPS_FILENAME     = "irreps_data.jld2"
@inline pathof_littlegroups_data(D::Integer) =
    joinpath(dirname(@__DIR__), "data/lgirreps", string(D)*"d", DATA_LITTLEGROUPS_FILENAME)
@inline pathof_lgirreps_data(D::Integer) =
    joinpath(dirname(@__DIR__), "data/lgirreps", string(D)*"d", DATA_LGIRREPS_FILENAME)

function _load_littlegroups_data(sgnum::Integer, jldfile::JLD2.JLDFile)   
    jldgroup = jldfile[string(sgnum)]
    sgops_str::Vector{String}      = jldgroup["sgops"]
    klabs::Vector{String}          = jldgroup["klab_list"]
    kstrs::Vector{String}          = jldgroup["kstr_list"]
    opsidxs::Vector{Vector{Int16}} = jldgroup["opsidx_list"]

    return sgops_str, klabs, kstrs, opsidxs
end

function _load_lgirreps_data(sgnum::Integer, jldfile::JLD2.JLDFile)
    jldgroup = jldfile[string(sgnum)] 
    # ≈ 70% of the time in loading all irreps is spent in getting Ps_list and τs_list
    Ps_list::Vector{Vector{Vector{Matrix{ComplexF64}}}}             = jldgroup["matrices_list"]
    τs_list::Vector{Vector{Union{Nothing,Vector{Vector{Float64}}}}} = jldgroup["translations_list"]
    types_list::Vector{Vector{Int64}}                               = jldgroup["type_list"]
    cdmls_list::Vector{Vector{String}}                              = jldgroup["cdml_list"]

    return Ps_list, τs_list, types_list, cdmls_list
end



# unexported character table convenience constructors (see also CharacterTable(::AbstractVector{<:AbstractIrrep})))
# TODO: Move these to types.jl and fix inconsistent method naming?
function chartable(klab::String, sgnum::Integer, Dᵛ::Val, αβγ=nothing)
    lgirsd = get_lgirreps(sgnum, Dᵛ)
    CharacterTable(lgirsd[klab], αβγ)
end
chartable(klab::String, sgnum::Integer, D::Integer=3, αβγ=nothing) = 
    chartable(klab, sgnum, Val(D), αβγ)


function chartable(kv::KVec, sgnum::Integer, Dᵛ::Val, αβγ=nothing)
    lgirsd = get_lgirreps(sgnum, Dᵛ)
    kidx = findfirst(x->kvec(first(x))==kv, lgirsd)
    if kidx === nothing
        throw(DomainError(kv, "Could not find input `kv` in the requested space group"))
    else
        return CharacterTable(lgirsd[kidx], αβγ)
    end
end
chartable(kv::KVec, sgnum::Integer, D::Integer=3, αβγ=nothing) = 
    chartable(kv, sgnum, Val(D), αβγ)


# old attempt at trying to have the data files open all the time, whenever the 
# module is called, and then closed afterwards. Ultimately, this only worked 
# rather sporadically and seemed quite buggy (though it was faster, since
# we didn't have to open the file every time we wanted to read from it (which,
# apparently, is quite expensive in JLD and JLD2))
# we want to keep the irrep files open whenever Crystalline is brought into play
# otherwise, we have to pay a large price to locate it etc.
#= 
   const IRREPS_DATA_FILE_3D = JLD2.jldopen((@__DIR__)*"/../data/lgirreps/3d/irreps_data.jld2", "r")
   atexit(()->close(IRREPS_DATA_FILE_3D))
   const LGS_DATA_FILE_3D = JLD2.jldopen((@__DIR__)* "/../data/lgirreps/3d/littlegroups_data.jld2", "r")
   atexit(()->close(LGS_DATA_FILE_3D))
=#