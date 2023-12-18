# Little group operations loading
"""
    littlegroups(sgnum::Integer, D::Union{Val{Int}, Integer}=Val(3)) 
                                                        -> Dict{String, LittleGroup{D}}

For given space group number `sgnum` and dimension `D`, return the associated little groups
(`LittleGroups{D}`s) at high-symmetry k-points, lines, and planes (see also
[`lgirreps`](@ref)).

Returns a `Dict` with little group **k**-point labels as keys and vectors of
`LittleGroup{D}`s as values.

## Notes
A conventional crystallographic setting is assumed (as in [`spacegroup`](@ref)).

Unlike `spacegroup`, "centering"-copies of symmetry operations are not included in the
returned `LittleGroup`s; as an example, space group 110 (body-centered, with centering
symbol 'I') has a centering translation `[1/2,1/2,1/2]` in the conventional setting:
the symmetry operations returned by `spacegroup` thus includes e.g. both `{1|0}` and 
`{1|½,½,½}` while the symmetry operations returned by `littlegroups` only include
`{1|0}` (and so on).

Currently, only `D = 3` is supported.

## References
The underlying data is sourced from the ISOTROPY dataset: see also [`lgirreps`](@ref).
"""
function littlegroups(sgnum::Integer, ::Val{D}=Val(3),
                          jldfile::JLD2.JLDFile=LGS_JLDFILES[D][]) where D
    D ∉ (1,2,3) && _throw_invalid_dim(D)

    sgops_str, klabs, kstrs, opsidxs = _load_littlegroups_data(sgnum, jldfile)

    sgops = SymOperation{D}.(sgops_str)
    lgs = Dict{String, LittleGroup{D}}()
    @inbounds for (klab, kstr, opsidx) in zip(klabs, kstrs, opsidxs)
        lgs[klab] = LittleGroup{D}(sgnum, KVec{D}(kstr), klab, sgops[opsidx])
    end
    return lgs
end
# convenience functions without Val(D) usage; avoid internally
littlegroups(sgnum::Integer, D::Integer) = littlegroups(sgnum, Val(D))

#------------------------------------------------------------------------------------------
# LGIrrep loading
"""
    lgirreps(sgnum::Integer, D::Union{Val{Int}, Integer}=Val(3))
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
function lgirreps(sgnum::Integer, Dᵛ::Val{D}=Val(3),
                      lgs_jldfile::JLD2.JLDFile=LGS_JLDFILES[D][],
                      irs_jldfile::JLD2.JLDFile=LGIRREPS_JLDFILES[D][]) where D
    D ∉ (1,2,3) && _throw_invalid_dim(D)
  
    lgs = littlegroups(sgnum, Dᵛ, lgs_jldfile)

    Ps_list, τs_list, realities_list, cdmls_list = _load_lgirreps_data(sgnum, irs_jldfile)

    lgirsd = Dict{String, Vector{LGIrrep{D}}}()
    for (Ps, τs, realities, cdmls) in zip(Ps_list, τs_list, realities_list, cdmls_list)
        klab = klabel(first(cdmls))
        lg   = lgs[klab]
        lgirsd[klab] = [LGIrrep{D}(cdml, lg, P, τ, Reality(reality)) for (P, τ, reality, cdml) in zip(Ps, τs, realities, cdmls)]
    end
    
    return lgirsd
end
lgirreps(sgnum::Integer, D::Integer) = lgirreps(sgnum, Val(D))


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

# ---------------------------------------------------------------------------------------- #
# Misc functions with `LGIrrep`



function ⊕(lgir1::LGIrrep{D}, lgir2::LGIrrep{D}) where D
    if position(lgir1) ≠ position(lgir2) || num(lgir1) ≠ num(lgir2) || order(lgir1) ≠ order(lgir2)
        error("The direct sum of two LGIrreps requires identical little groups")
    end
    if lgir1.translations ≠ lgir2.translations
        error("The provided LGIrreps have different translation-factors and cannot be \
               combined within a single translation factor system")
    end
    
    cdml = label(lgir1)*"⊕"*label(lgir2)
    g   = group(lgir1)
    T   = eltype(eltype(lgir1.matrices))
    z12 = zeros(T, irdim(lgir1), irdim(lgir2))
    z21 = zeros(T, irdim(lgir2), irdim(lgir1))
    matrices = [[m1 z12; z21 m2] for (m1, m2) in zip(lgir1.matrices, lgir2.matrices)]
    translations = lgir1.translations
    reality = UNDEF
    iscorep = lgir1.iscorep || lgir2.iscorep

    return LGIrrep{D}(cdml, g, matrices, translations, reality, iscorep)
end