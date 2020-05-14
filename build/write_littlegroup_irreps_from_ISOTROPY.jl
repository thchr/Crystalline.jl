using Crystalline, JLD2

"""
    __write_littlegroupirreps(LGIRS::Vector{Vector{NTuple{N,LGIrrep}}})
                                                    --> Nothing

Write all little group small irreps associated with all specific space 
group to disk, as JSON files, to ease subsequent loading of little group 
small irreps. Takes a vector of a vector of little group small irreps of 
the sort
    `LGIRS::Vector{Vector{NTuple{N,LGIrrep}}}`
i.e., vector-indexed across space group number, then vector-indexed across 
distinct k-points, and finally tuple-indexed across distinct irreps; in 
practice, calling 
    `__write_littlegroupirreps()`
will load LGIRS via `parselittlegroupirreps` write **all** the little group
irreps to disk. There is generally no reason for a user to **ever** do this.

Returns the filepath of the saved .jld2 files
"""
function __write_littlegroupirreps(LGIRS)
    savepath = (@__DIR__)*"/../data/lgirreps/3d/"
    filename_lgs    = savepath*"/littlegroups_data"
    filename_irreps = savepath*"/irreps_data"
    
    JLD2.jldopen(filename_lgs*".jld2", "w") do littlegroups_file
    JLD2.jldopen(filename_irreps*".jld2", "w") do irreps_file

    for lgirsvec in LGIRS # fixed sgnum: lgirsvec has structure lgirsvec[kidx][iridx]
        sgnum = num(first(first(lgirsvec)))
        Nk = length(lgirsvec) 

        # ==== little groups data ====
        sgops = operations(first(first(lgirsvec)))
        klab_list = Vector{String}(undef, Nk)
        kstr_list = Vector{String}(undef, Nk)
        opsidx_list  = Vector{Vector{Int16}}(undef, Nk) # Int16 because number of ops ≤ 192; ideally, would use UInt8
        for (kidx, lgirs) in enumerate(lgirsvec) # lgirs is a tuple of LGIrreps, all at the same 𝐤-point
            lgir = first(lgirs) # 𝐤-info is the same for each LGIrrep in tuple lgirs
            klab_list[kidx] = klabel(lgir)
            kstr_list[kidx] = filter(!isspace, chop(string(kvec(lgir)); head=1, tail=1))
            opsidx_list[kidx] = map(y->findfirst(==(y), sgops), operations(lgir))
        end
    
        # ==== irreps data ====
        matrices_list = [[lgir.matrices for lgir in lgirs] for lgirs in lgirsvec]
        #translations_list = [[lgir.translations for lgir in lgirs] for lgirs in lgirsvec]
        # don't want to save a bunch of zeros if all translations are zero: 
        # instead, save `nothing` as a sentinel value
        translations_list = [Union{Nothing, Vector{Vector{Float64}}}[
                                    all(iszero, translations(lgir)) ? nothing : translations(lgir) 
                                    for lgir in lgirs] for lgirs in lgirsvec] # dreadful generator, but OK...
        type_list = [[type(lgir) for lgir in lgirs] for lgirs in lgirsvec]
        cdml_list = [[label(lgir) for lgir in lgirs] for lgirs in lgirsvec]
        
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