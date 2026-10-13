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
# actually write .jld files for 3D, then 2D and 1D
# (to do this, we must close `LGIRREPS_JLDFILES` and `LGS_JLDFILES` - which are opened upon
# initialization of Crystalline - before writing new content to them; to that end, we simply
# close the files below (make sure you don't have Crystalline loaded in another session
# though..!))
close(Crystalline.LGIRREPS_JLDFILES[3][])
close(Crystalline.LGS_JLDFILES[3][])

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
filename_lgs, filename_irreps = __write_littlegroupirreps(LGIRS_3D)
# reopen the 3D files, from which the 2D and 1D irreps are derived below
Crystalline.LGS_JLDFILES[3][] = JLD2.jldopen(filename_lgs, "r")
Crystalline.LGIRREPS_JLDFILES[3][] = JLD2.jldopen(filename_irreps, "r")

# 2D and 1D (little groups from tabulated k-points; irreps derived from parent space groups,
# single- and double-valued)
include(joinpath(@__DIR__, "line_plane_littlegroups.jl"))
include(joinpath(@__DIR__, "line_plane_irreps.jl"))
using .LinePlaneLittleGroups, .LinePlaneIrreps
for D in (2, 1)
    write_littlegroups(Val(D))
    write_lgirreps(Val(D); spinful=Val(false))
    write_lgirreps(Val(D); spinful=Val(true))
end