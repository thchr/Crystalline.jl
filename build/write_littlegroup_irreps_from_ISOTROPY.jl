using Crystalline, JLD2

"""
    __write_littlegroupirreps(LGIRS::Vector{Dict{Vector{LGIrrep}}})
                                                    --> ::String, ::String

Write all little group/small irreps, i.e. `LGIrrep`s, in the input to disk, as JSON files,
in order to ease subsequent loading of `LGIrrep`s. Input is of the type

    `LGIRS::Vector{Dict{Vector{LGIrrep}}}`

intended as vector-indexed across space group number, then dict-indexed across distinct
k-point labels, and finally vector-indexed across distinct irreps; in practice, calling 

    `__write_littlegroupirreps()`

will load `LGIRS` via `parselittlegroupirreps` and write **all** the `LGIrrep`s in ISOTROPY.
There is generally no reason for a user to **ever** do this.

Returns the filepath of the saved .jld2 files.
"""
function __write_littlegroupirreps(LGIRS)
    savepath = (@__DIR__)*"/../data/lgirreps/3d/"
    filename_lgs    = savepath*"/littlegroups_data"
    filename_irreps = savepath*"/irreps_data"
    
    JLD2.jldopen(filename_lgs*".jld2", "w") do littlegroups_file
    JLD2.jldopen(filename_irreps*".jld2", "w") do irreps_file

    for lgirsd in LGIRS # fixed sgnum: lgirsd has structure lgirsd[klab][iridx]
        sgnum = num(first(first(lgirsd)))
        Nk = length(lgirsd) 

        # ==== little groups data ====
        sgops = operations(first(lgirsd["Γ"]))
        klab_list = Vector{String}(undef, Nk)
        kstr_list = Vector{String}(undef, Nk)
        opsidx_list  = Vector{Vector{Int16}}(undef, Nk) # Int16 because number of ops ≤ 192; ideally, would use UInt8
        for (kidx, (klab, lgirs)) in enumerate(lgirsd) # lgirs is a vector of LGIrreps, all at the same 𝐤-point
            lgir = first(lgirs) # 𝐤-info is the same for each LGIrrep in vector lgirs
            klab_list[kidx] = klab
            kstr_list[kidx] = filter(!isspace, chop(string(kvec(lgir)); head=1, tail=1))
            opsidx_list[kidx] = map(y->findfirst(==(y), sgops), operations(lgir))
        end
    
        # ==== irreps data ====
        matrices_list = [[lgir.matrices for lgir in lgirs] for lgirs in values(lgirsd)]
        #translations_list = [[lgir.translations for lgir in lgirs] for lgirs in values(lgirsd)]
        # don't want to save a bunch of zeros if all translations are zero: 
        # instead, save `nothing` as a sentinel value
        translations_list = [Union{Nothing, Vector{Vector{Float64}}}[
                                    all(iszero, translations(lgir)) ? nothing : translations(lgir) 
                                    for lgir in lgirs] for lgirs in values(lgirsd)] # dreadful generator, but OK...
        type_list = [[type(lgir) for lgir in lgirs] for lgirs in values(lgirsd)]
        cdml_list = [[label(lgir) for lgir in lgirs] for lgirs in values(lgirsd)]
        
        # ==== save data ====
        # little groups
        littlegroups_file[string(sgnum)*"/sgops"]       = xyzt.(sgops)
        littlegroups_file[string(sgnum)*"/klab_list"]   = klab_list
        littlegroups_file[string(sgnum)*"/kstr_list"]   = kstr_list
        littlegroups_file[string(sgnum)*"/opsidx_list"] = opsidx_list

        # irreps
        irreps_file[string(sgnum)*"/matrices_list"] = matrices_list
        irreps_file[string(sgnum)*"/translations_list"] = translations_list
        irreps_file[string(sgnum)*"/type_list"] = type_list
        irreps_file[string(sgnum)*"/cdml_list"] = cdml_list
    end # end of loop

    end # close irreps_file
    end # close irreps_file

    return filename_lgs, filename_irreps
end
__write_littlegroupirreps() = __write_littlegroupirreps(parselittlegroupirreps())

__write_littlegroupirreps()