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
                          jldfile::JLD2.JLDFile=LGS_JLDFILES[D]) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)

    sgops_str, klabs, kstrs, opsidxs = _load_littlegroups_data(sgnum, jldfile)

    sgops = SymOperation{D}.(sgops_str)
    lgs = Dict{String, LittleGroup{D}}()
    @inbounds for (klab, kstr, opsidx) in zip(klabs, kstrs, opsidxs)
        lgs[klab] = LittleGroup{D}(sgnum, KVec{D}(kstr), klab, sgops[opsidx])
    end
    return lgs
end
# convenience functions without Val(D) usage; avoid internally
get_littlegroups(sgnum::Integer, D::Integer) = get_littlegroups(sgnum, Val(D))

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
- The returned irreps are complex in general. Real irreps (as needed in time-reversal
  invariant settings) can subsequently be obtained with the [`realify`](@ref) method.
- Returned irreps are spinless.
- The irrep labelling follows CDML conventions.
- Irreps along lines or planes may depend on free parameters `αβγ` that parametrize the
  **k** point. To evaluate the irreps at a particular value of `αβγ` and return the
  associated matrices, use `(lgir::LGIrrep)(αβγ)`. If `αβγ` is an empty tuple in this call,
  the matrices associated with `lgir` will be evaluated assuming `αβγ = [0,0,...]`.

## References
The underlying data is sourced from the ISOTROPY ISO-IR dataset. Please cite original
reference material associated with ISO-IR:

1. Stokes, Hatch, & Campbell, 
   [ISO-IR, ISOTROPY Software Suite](https://stokes.byu.edu/iso/irtables.php).
2. Stokes, Campbell, & Cordes,
   [Acta Cryst. A. **69**, 388-395 (2013)](https://doi.org/10.1107/S0108767313007538).

The ISO-IR dataset is occasionally missing some **k**-points that lie outside the basic
domain but still resides in the representation domain (i.e. **k**-points with postscripted
'A', 'B', etc. labels, such as 'ZA'). In such cases, the missing irreps may instead have
been manually sourced from the Bilbao Crystallographic Database.
"""
function get_lgirreps(sgnum::Integer, Dᵛ::Val{D}=Val(3),
                      lgs_jldfile::JLD2.JLDFile=LGS_JLDFILES[D],
                      irs_jldfile::JLD2.JLDFile=LGIRREPS_JLDFILES[D]) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)
  
    lgs = get_littlegroups(sgnum, Dᵛ, lgs_jldfile)

    Ps_list, τs_list, realities_list, cdmls_list = _load_lgirreps_data(sgnum, irs_jldfile)

    lgirsd = Dict{String, Vector{LGIrrep{D}}}()
    for (Ps, τs, realities, cdmls) in zip(Ps_list, τs_list, realities_list, cdmls_list)
        klab = klabel(first(cdmls))
        lg   = lgs[klab]
        lgirsd[klab] = [LGIrrep{D}(cdml, lg, P, τ, Reality(reality)) for (P, τ, reality, cdml) in zip(Ps, τs, realities, cdmls)]
    end
    
    return lgirsd
end
get_lgirreps(sgnum::Integer, D::Integer) = get_lgirreps(sgnum, Val(D))


# ===== utility functions (loads raw data from the harddisk) =====
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
    realities_list::Vector{Vector{Int8}}                            = jldgroup["realities_list"]
    cdmls_list::Vector{Vector{String}}                              = jldgroup["cdml_list"]

    return Ps_list, τs_list, realities_list, cdmls_list
end


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