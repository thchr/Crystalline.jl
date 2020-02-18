const JldOrNothing = Union{Nothing,JLD2.JLDFile}

# Little group operations loading
function get_littlegroups(sgnum::Integer, ::Val{D},
                          jldfile::JldOrNothing=nothing) where D
    D ≠ 3 && throw_2d_not_yet_implemented(D)

    sgops_str, klabs, kstrs, opsidxs = if isnothing(jldfile)
        JLD2.jldopen(DATA_PATH_LITTLEGROUPS_3D, "r") do jldfile
             _load_littlegroups_data(sgnum, jldfile)
        end
    else
        _load_littlegroups_data(sgnum, jldfile)
    end

    sgops = SymOperation{D}.(sgops_str)
    Nk = length(klabs)
    lgs = Vector{LittleGroup{D}}(undef, Nk)
    @inbounds for kidx in Base.OneTo(Nk)
        lgs[kidx] = LittleGroup{D}(sgnum, KVec(kstrs[kidx]), klabs[kidx], 
                                   sgops[opsidxs[kidx]])
    end
    return lgs
end

function get_all_littlegroups(::Val{D}) where D
    JLD2.jldopen(SGOps.DATA_PATH_LITTLEGROUPS_3D,"r") do lgfile
        return [get_littlegroups(sgnum, Val(D), lgfile) for sgnum in Base.OneTo(MAX_SGNUM[D])]
    end
end
# convenience functions without Val(D) usage; avoid internally
get_littlegroups(sgnum::Integer, D::Integer=3, jldfile::JldOrNothing=nothing) = get_littlegroups(sgnum, Val(D), jldfile)
get_all_littlegroups(D::Integer=3) = get_all_littlegroups(Val(D))

# Little group irrep loading
function get_lgirreps(sgnum::Integer, Dᵛ::Val{D}, lgs_jldfile::JldOrNothing=nothing,
                      irs_jldfile::JldOrNothing=nothing) where D
    D ≠ 3 && throw_2d_not_yet_implemented(D)
  
    lgs = get_littlegroups(sgnum, Dᵛ, lgs_jldfile)

    Ps_list, τs_list, type_list, cdml_list = if isnothing(irs_jldfile)
        JLD2.jldopen(DATA_PATH_LGIRREPS_3D, "r") do irs_jldfile
            _load_irreps_data(sgnum, irs_jldfile)
        end
    else
        _load_irreps_data(sgnum, irs_jldfile)
    end

    lgirsvec = Vector{Vector{LGIrrep{D}}}(undef, length(lgs))
    @inbounds for (kidx, lg) in enumerate(lgs)
        Nirr = length(type_list[kidx])
        lgirsvec[kidx] = Vector{LGIrrep{D}}(undef, Nirr)
        @inbounds for iridx in Base.OneTo(Nirr)
            lgirsvec[kidx][iridx] = LGIrrep{D}(cdml_list[kidx][iridx],
                                               lg, 
                                               Ps_list[kidx][iridx], 
                                               τs_list[kidx][iridx], 
                                               type_list[kidx][iridx])
        end
    end
    
    return lgirsvec
end
function get_lgirreps(sgnum::Integer, D::Integer=3, lgs_jldfile::JldOrNothing=nothing, 
                      irs_jldfile::JldOrNothing=nothing)
    get_lgirreps(sgnum, Val(D), lgs_jldfile, irs_jldfile)
end

function get_all_lgirreps(Dᵛ::Val{D}) where D
    JLD2.jldopen(SGOps.DATA_PATH_LITTLEGROUPS_3D,"r") do lgfile;
        JLD2.jldopen(SGOps.DATA_PATH_LGIRREPS_3D,"r") do irfile;
            return [get_lgirreps(sgnum, Dᵛ, lgfile, irfile) for sgnum in Base.OneTo(MAX_SGNUM[D])]; 
        end
    end
end
get_all_lgirreps(D::Integer=3) = get_all_lgirreps(Val(D))

# ===== utility functions (loads raw data from the harddisk) =====
const DATA_PATH_LITTLEGROUPS_3D = (@__DIR__)*"/../data/lgirreps/3d/littlegroups_data.jld2"
const DATA_PATH_LGIRREPS_3D = (@__DIR__)*"/../data/lgirreps/3d/irreps_data.jld2"
function _load_littlegroups_data(sgnum::Integer, jldfile::JLD2.JLDFile)   
    jldgroup = jldfile[string(sgnum)]
    sgops_str::Vector{String}      = jldgroup["sgops"]
    klabs::Vector{String}          = jldgroup["klab_list"]
    kstrs::Vector{String}          = jldgroup["kstr_list"]
    opsidxs::Vector{Vector{Int16}} = jldgroup["opsidx_list"]

    return sgops_str, klabs, kstrs, opsidxs
end

function _load_irreps_data(sgnum::Integer, jldfile::JLD2.JLDFile)
    jldgroup = jldfile[string(sgnum)] 
    # ≈ 70% of the time in loading all irreps is spent in getting Ps_list and τs_list
    Ps_list::Vector{Vector{Vector{Matrix{ComplexF64}}}}             = jldgroup["matrices_list"]
    τs_list::Vector{Vector{Union{Nothing,Vector{Vector{Float64}}}}} = jldgroup["translations_list"]
    type_list::Vector{Vector{Int64}}                                = jldgroup["type_list"]
    cdml_list::Vector{Vector{String}}                               = jldgroup["cdml_list"]

    return Ps_list, τs_list, type_list, cdml_list
end

throw_2d_not_yet_implemented(D) = throw(DomainError(D, "Only D=3 little groups and irreps are supported at this time"))




# character table construction
function chartable(lgirs::AbstractVector{LGIrrep{D}}) where D
    table = Array{ComplexF64}(undef, length(lgirs), order(first(lgirs)))
    for (i,row) in enumerate(eachrow(table))
        row .= characters(lgirs[i])
    end
    tag = join(["#", string(num(first(lgirs)))])
    return CharacterTable{D}(operations(first(lgirs)), label.(lgirs), table, tag)
end

function chartable(klab::String, sgnum::Integer, Dᵛ::Val)
    lgirsvec = get_lgirreps(sgnum, Dᵛ)
    kidx = findfirst(x->klabel(first(x))==klab, lgirsvec)
    if kidx === nothing
        throw(DomainError(klab, "Could not find the input klabel `klab` in the requested space group"))
    else
        return chartable(lgirsvec[kidx])
    end
end
chartable(klab::String, sgnum::Integer, D::Integer=3) = chartable(klab, sgnum, Val(D))


function chartable(kv::KVec, sgnum::Integer, Dᵛ::Val)
    lgirsvec = get_lgirreps(sgnum, Dᵛ)
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
chartable(kv::KVec, sgnum::Integer, D::Integer=3) = chartable(kv, sgnum, Val(D))


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