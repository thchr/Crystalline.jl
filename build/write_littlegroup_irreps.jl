using Crystalline, JLD2

"""
    __write_littlegroupirreps(LGIRS::Vector{Dict{Vector{LGIrrep}}})
                                                    --> ::String, ::String

Write all little group/small irreps, i.e. `LGIrrep`s, in the input to disk, as JLD2 files,
in order to ease subsequent loading of `LGIrrep`s. Input is of the type

    `LGIRS::Vector{Dict{Vector{LGIrrep}}}`

intended as vector-indexed across space group number, then dict-indexed across distinct
**k**-point labels, and finally vector-indexed across distinct irreps; in practice, calling 

    `__write_littlegroupirreps()`

will load `LGIRS` via `parselittlegroupirreps` and write **all** the `LGIrrep`s in ISOTROPY.
There is generally no reason for a user to **ever** do this.

Returns the filepath of the saved .jld2 files.
"""
function __write_littlegroupirreps(
            LGIRS::Vector{Dict{String, <:AbstractVector{LGIrrep{D}}}}) where D

    savepath = (@__DIR__)*"/../data/irreps/lgs/"*string(D)*"d"
    filename_lgs    = joinpath(savepath, "littlegroups_data.jld2")
    filename_irreps = joinpath(savepath, "irreps_data.jld2")
    
    JLD2.jldopen(filename_lgs, "w") do littlegroups_file
    JLD2.jldopen(filename_irreps, "w") do irreps_file

    for lgirsd in LGIRS # fixed sgnum: lgirsd has structure lgirsd[klab][iridx]
        sgnum = num(first(lgirsd["Γ"]))
        Nk = length(lgirsd) 

        # ==== little groups data ====
        sgops = operations(first(lgirsd["Γ"]))
        klab_list = Vector{String}(undef, Nk)
        kstr_list = Vector{String}(undef, Nk)
        opsidx_list  = Vector{Vector{Int16}}(undef, Nk) # Int16 because number of ops ≤ 192; ideally, would use UInt8
        for (kidx, (klab, lgirs)) in enumerate(lgirsd) # lgirs is a vector of LGIrreps, all at the same 𝐤-point
            lgir = first(lgirs) # 𝐤-info is the same for each LGIrrep in vector lgirs
            klab_list[kidx] = klab
            kstr_list[kidx] = filter(!isspace, chop(string(position(lgir)); head=1, tail=1))
            opsidx_list[kidx] = map(y->findfirst(==(y), sgops), operations(lgir))
        end
    
        # ==== irreps data ====
        matrices_list = [[lgir.matrices for lgir in lgirs] for lgirs in values(lgirsd)]
        #translations_list = [[lgir.translations for lgir in lgirs] for lgirs in values(lgirsd)]
        # don't want to save a bunch of zeros if all translations are zero: 
        # instead, save `nothing` as a sentinel value
        translations_list = [Union{Nothing, Vector{Vector{Float64}}}[
                                    all(iszero, Crystalline.translations(lgir)) ? nothing : Crystalline.translations(lgir) 
                                    for lgir in lgirs] for lgirs in values(lgirsd)] # dreadful generator, but OK...
        realities_list = [[Integer(reality(lgir)) for lgir in lgirs] for lgirs in values(lgirsd)]
        cdml_list      = [[label(lgir) for lgir in lgirs] for lgirs in values(lgirsd)]
        
        # ==== save data ====
        # little groups
        littlegroups_file[string(sgnum)*"/sgops"]       = xyzt.(sgops)
        littlegroups_file[string(sgnum)*"/klab_list"]   = klab_list
        littlegroups_file[string(sgnum)*"/kstr_list"]   = kstr_list
        littlegroups_file[string(sgnum)*"/opsidx_list"] = opsidx_list

        # irreps
        irreps_file[string(sgnum)*"/matrices_list"]     = matrices_list
        irreps_file[string(sgnum)*"/translations_list"] = translations_list
        irreps_file[string(sgnum)*"/realities_list"]    = realities_list # ::Vector{Int8}
        irreps_file[string(sgnum)*"/cdml_list"]         = cdml_list
    end # end of loop

    end # close irreps_file
    end # close littlegroups_file

    return filename_lgs, filename_irreps
end


# ---------------------------------------------------------------------------------------- #
# 1D line groups: little groups irreps written down manually

LGIRS_1D = [Dict{String, Vector{LGIrrep{1}}}() for _ in 1:2]
function make_1d_lgirrep(sgnum::Integer, klab::String, cdml_suffix::String,
                         kx, ops::Vector{SymOperation{1}}, scalars::Vector{<:Number}=[1.0,],
                         reality_type::Reality=REAL)
    cdml = klab*cdml_suffix
    lg   = LittleGroup{1}(sgnum, KVec{1}(string(kx)), klab, ops)
    matrices = [fill(ComplexF64(v), 1,1) for v in scalars]
    translations = [zeros(1) for _ in scalars]
    return LGIrrep{1}(cdml, lg, matrices, translations, reality_type, false)
end

# Line group 1
LGIRS_1D[1]["Γ"] = [make_1d_lgirrep(1, "Γ", "₁", 0,    [S"x"],        [1.0])]
LGIRS_1D[1]["X"] = [make_1d_lgirrep(1, "X", "₁", 0.5,  [S"x"],        [1.0])]
LGIRS_1D[1]["Ω"] = [make_1d_lgirrep(1, "Ω", "₁", "u",  [S"x"],        [1.0], COMPLEX)]

# Line group 2
LGIRS_1D[2]["Γ"] = [make_1d_lgirrep(2, "Γ", "₁⁺", 0,   [S"x", S"-x"], [1.0, 1.0]),   # even
                    make_1d_lgirrep(2, "Γ", "₁⁻", 0,   [S"x", S"-x"], [1.0, -1.0])]  # odd
LGIRS_1D[2]["X"] = [make_1d_lgirrep(2, "X", "₁⁺", 0.5, [S"x", S"-x"], [1.0, 1.0]),   # even
                    make_1d_lgirrep(2, "X", "₁⁻", 0.5, [S"x", S"-x"], [1.0, -1.0])]  # odd
LGIRS_1D[2]["Ω"] = [make_1d_lgirrep(2, "Ω", "₁", "u",  [S"x"],        [1.0])]

# ---------------------------------------------------------------------------------------- #
# 2D plane groups: little groups irreps extracted for symmorphic groups via point groups

include("setup_2d_littlegroup_irreps.jl") # defines the variable LGIRS_2D′

# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# actually write .jld files for 1D and 3D
# (to do this, we must close `LGIRREPS_JLDFILES` and `LGS_JLDFILES` - which are opened upon
# initialization of Crystalline - before writing new content to them; to that end, we simply
# close the files below (make sure you don't have Crystalline loaded in another session
# though..!))
foreach(jldfile -> close(jldfile[]), Crystalline.LGIRREPS_JLDFILES)
foreach(jldfile -> close(jldfile[]), Crystalline.LGS_JLDFILES)

# 3D (from ISOTROPY)
include(joinpath((@__DIR__), "ParseIsotropy.jl")) # load the ParseIsotropy module
using Main.ParseIsotropy                                # (exports `parselittlegroupirreps`)
LGIRS_3D = parselittlegroupirreps()
# ... ISOTROPY is missing several irreps; we bring those in below, obtained from manual
# transcription of irreps from Bilbao; script below defines `LGIRS_add` which stores these
# manual additions (a Dict with `sgnum` keys)
include(joinpath((@__DIR__), "..", "data/irreps/lgs/manual_lgirrep_additions.jl"))
for (sgnum, lgirsd_add) in LGIRS_add # merge Bilbao additions with ISOTROPY irreps
    merge!(LGIRS_3D[sgnum], lgirsd_add)
end
__write_littlegroupirreps(LGIRS_3D)

# 2D (from point group matching)
__write_littlegroupirreps(LGIRS_2D′)

# 1D (manual)
__write_littlegroupirreps(LGIRS_1D)