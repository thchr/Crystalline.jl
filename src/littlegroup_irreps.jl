function get_littlegroups(sgnum::Integer, dim::Integer=3,
                          jldfile::Union{Nothing,JLD2.JLDFile}=nothing)
    sgops, klabs, kstrs, opsidxs = _load_littlegroups_data(sgnum, dim, jldfile)
    Nk = length(klabs)
    lgs = Vector{LittleGroup{dim}}(undef, Nk)
    @inbounds for kidx in Base.OneTo(Nk)
        lgs[kidx] = LittleGroup{dim}(sgnum, KVec(kstrs[kidx]), klabs[kidx], sgops[opsidxs[kidx]])
    end
    return lgs
end

function get_all_littlegroups(dim::Integer=3)
    if dim == 3
        sgnums = 1:230
    else 
        throw_2d_not_yet_implemented()
    end
    JLD2.jldopen(SGOps.DATA_PATH_LITTLEGROUPS_3D,"r") do lgfile;
        return [get_littlegroups(sgnum, dim, lgfile) for sgnum in sgnums]; 
    end
end

function get_lgirreps(sgnum::Integer, dim::Integer=3,
                      lgs_jldfile::Union{Nothing,JLD2.JLDFile}=nothing,
                      irs_jldfile::Union{Nothing,JLD2.JLDFile}=nothing)
    lgs = get_littlegroups(sgnum, dim, lgs_jldfile)
    Ps_list, τs_list, type_list, cdml_list = _load_irreps_data(sgnum, dim, irs_jldfile)

    lgirsvec = Vector{Vector{LGIrrep{dim}}}(undef, length(lgs))
    @inbounds for (kidx, lg) in enumerate(lgs)
        Nirr = length(type_list[kidx])
        lgirsvec[kidx] = Vector{LGIrrep{dim}}(undef, Nirr)
        @inbounds for iridx in Base.OneTo(Nirr)
            lgirsvec[kidx][iridx] = LGIrrep{dim}(cdml_list[kidx][iridx],
                                                 lg, 
                                                 Ps_list[kidx][iridx], 
                                                 τs_list[kidx][iridx], 
                                                 type_list[kidx][iridx])
        end
    end
    
    return lgirsvec
end

function get_all_lgirreps(dim::Integer=3)
    if dim == 3
        sgnums = 1:230
    else 
        throw_2d_not_yet_implemented()
    end

    JLD2.jldopen(SGOps.DATA_PATH_LITTLEGROUPS_3D,"r") do lgfile;
        JLD2.jldopen(SGOps.DATA_PATH_LGIRREPS_3D,"r") do irfile;
            return [get_lgirreps(sgnum, dim, lgfile, irfile) for sgnum in sgnums]; 
        end
    end
end

# ===== utility functions (does the actual loading from the harddisk) =====
const DATA_PATH_LITTLEGROUPS_3D = (@__DIR__)*"/../data/lgirreps/3d/littlegroups_data.jld2"
const DATA_PATH_LGIRREPS_3D = (@__DIR__)*"/../data/lgirreps/3d/irreps_data.jld2"
function _load_littlegroups_data(sgnum::Integer, dim::Integer, ::Nothing)    
    JLD2.jldopen(DATA_PATH_LITTLEGROUPS_3D, "r") do jldfile
        return _load_littlegroups_data(sgnum, dim, jldfile)
    end   
end
function _load_littlegroups_data(sgnum::Integer, dim::Integer, jldfile::JLD2.JLDFile)
    dim ≠ 3 && throw_2d_not_yet_implemented()
    
    jldgroup = jldfile[string(sgnum)]
    sgops_str::Vector{String}      = jldgroup["sgops"]
    klabs::Vector{String}          = jldgroup["klab_list"]
    kstrs::Vector{String}          = jldgroup["kstr_list"]
    opsidxs::Vector{Vector{Int16}} = jldgroup["opsidx_list"]

    return SymOperation.(sgops_str), klabs, kstrs, opsidxs
end

function _load_irreps_data(sgnum::Integer, dim::Integer, ::Nothing)
    JLD2.jldopen(DATA_PATH_LGIRREPS_3D, "r") do jldfile
        return _load_irreps_data(sgnum, dim, jldfile)
    end
end
function _load_irreps_data(sgnum::Integer, dim::Integer, jldfile::JLD2.JLDFile)
    dim ≠ 3 && throw_2d_not_yet_implemented()

    jldgroup = jldfile[string(sgnum)] 
    # ≈ 70% of the time in loading all irreps is spent in getting Ps_list and τs_list
    Ps_list::Vector{Vector{Vector{Matrix{ComplexF64}}}}             = jldgroup["matrices_list"]
    τs_list::Vector{Vector{Union{Nothing,Vector{Vector{Float64}}}}} = jldgroup["translations_list"]
    type_list::Vector{Vector{Int64}}                                = jldgroup["type_list"]
    cdml_list::Vector{Vector{String}}                               = jldgroup["cdml_list"]

    return Ps_list, τs_list, type_list, cdml_list
end

throw_2d_not_yet_implemented() = throw(DomainError("Only 3D little groups are supported at this time"))

# character table construction
function chartable(lgirs::AbstractVector{LGIrrep{D}}) where D
    table = Array{ComplexF64}(undef, length(lgirs), order(first(lgirs)))
    for (i,row) in enumerate(eachrow(table))
        row .= characters(lgirs[i])
    end
    tag = join(["#", string(num(first(lgirs)))])
    return CharacterTable{D}(operations(first(lgirs)), label.(lgirs), table, tag)
end

function chartable(klab::String, sgnum::Integer, dim::Integer)
    lgirsvec = get_lgirreps(sgnum, dim)
    kidx = findfirst(x->klabel(first(x))==klab, lgirsvec)
    if kidx === nothing
        throw(DomainError(klab, "Could not find the input klabel `klab` in the requested space group"))
    else
        return chartable(lgirsvec[kidx])
    end
end

function chartable(kv::KVec, sgnum::Integer, dim::Integer)
    lgirsvec = get_lgirreps(sgnum, dim)
    # TODO: Implement matching to generic KVec format, so that 
    #       we can specify `kv` at concrete non-special momenta 
    #       and still match
    kidx = findfirst(x->kvec(first(x))==kv, lgirsvec)
    if kidx === nothing
        throw(DomainError(kv, "Could not find the input k-vector `kv` in the requested space group"))
    else
        return chartable(lgirsvec[kidx])
    end
end


# old attempt at trying to have the data files open all the time, whenever the 
# module is called, and then closed afterwards. Ultimately, this only worked 
# rather sporadically and seemed quite buggy (though it was faster, since
# we didn't have to open the file every time we wanted to read from it (which,
# apparently, is quite expensive in JLD and JLD2))
# we want to keep the irrep files open whenever SGOps is brought into play
# otherwise, we have to pay a large price to locate it etc.
#= 
   const IRREPS_DATA_FILE_3D = JLD2.jldopen((@__DIR__)*"/../data/lgirreps/3d/irreps_data.jld2", "r")
   atexit(()->close(IRREPS_DATA_FILE_3D))
   const LGS_DATA_FILE_3D = JLD2.jldopen((@__DIR__)* "/../data/lgirreps/3d/littlegroups_data.jld2", "r")
   atexit(()->close(LGS_DATA_FILE_3D))
=#