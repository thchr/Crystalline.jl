using LinearAlgebra
using Crystalline
using StaticArrays
import Crystalline:(get_lgirreps)

function get_lgirreps(sgnum::Int64, D::Int64, timereversal::Bool)
    if timereversal==false
        get_lgirreps(sgnum, D) 
    else 
        return Dict(klab => realify(lgirs) for (klab, lgirs) in get_lgirreps(sgnum, D))
    end
end

function group_representative(symop::SymOperation, sg::SpaceGroup)
    findfirst(≈(symop), sg)
end

function _find_isomorphic_parent_pointgroup(G)
    D = dim(G)
    ctᴳ = MultTable(G).table
    @inbounds for iuclab in PGS_IUCs[D]
        P = pointgroup(iuclab, D)
        ctᴾ = MultTable(P).table
        if ctᴳ == ctᴾ # bit sloppy; would miss ismorphisms that are "concealed" by row/column swaps
            return P
        end
    end
    return nothing # in case we didn't find any isomorphic parent
end

"This function gives the site_symmetry irreps for the Site Symmetry Group"
function site_symmetry_irreps(sitegroup::SiteGroup{D}, timereversal::Bool) where D
    site_parent_pointgroup=_find_isomorphic_parent_pointgroup(sitegroup)

    pgirs = get_pgirreps(site_parent_pointgroup.label, Val(D))
    timereversal==false ? pgirs : realify(pgirs)
end

"The function below gives the character table for the site symmetry group provided"
function site_symmetry_irreps_characters(sitegroup::SiteGroup{D}, timereversal::Bool) where D
    irreps=site_symmetry_irreps(sitegroup, timereversal)
    character_table=Matrix{Complex{Float64}}(undef, size(irreps[1].matrices)[1], size(irreps)[1] )
    size(character_table)
    for j in 1:size(irreps)[1]
        for i in 1:size(irreps[1].matrices)[1]
            character_table[i, j]=tr(irreps[j].matrices[i]) 
        end
    end
character_table
end


"This function provides the induced band representation"
function induced_band_representation(sitegroup::SiteGroup{D}, h::SymOperation{D}, ρ::Int64, k::KVec, timereversal::Bool) where D
    orbit_of_sitegroup=orbit(sitegroup, sitegroup.wp);
    length_of_orbit=length(orbit_of_sitegroup)
    #Now we sum over all the wyckoff positions in the orbit  
    gα=cosets(sitegroup)
    χ_G_k_G=0;
    wp=sitegroup.wp;
    wp_rvec=wp.qv.cnst
    gα_inverse=gα[1]

    transformed_k=transform(k, rotation(h)).cnst
    irrep_matrices=site_symmetry_irreps(sitegroup, timereversal)[ρ].matrices
    for α′ in 1:length_of_orbit 

        gα′=gα[α′]
        tα′α′=rotation(h)*orbit_of_sitegroup[α′].qv.cnst+translation(h)-orbit_of_sitegroup[α′].qv.cnst
        symop_t = SymOperation(-tα′α′)
        gα_inverse=inv(gα′)    

        gα′⁻¹ggα′ = compose(gα_inverse, compose(symop_t, compose(h,  gα[α′], false), false ), false)
        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), sitegroup)

        if site_symmetry_index !== nothing
            χ_G_k_G += cis(2π*dot(transformed_k, tα′α′))*tr(irrep_matrices[site_symmetry_index])
        end
    end
    return χ_G_k_G
end

function subduce_onto_littlegroups(k::String, kvec::KVec, ρ::Int64, sitegroup::SiteGroup{D}, timereversal::Bool) where D
    sgnum=num(sitegroup)
    sg=spacegroup(sgnum, D)
    lgirs=get_lgirreps(sgnum, D, timereversal)
    character_table=CharacterTable(lgirs[k]).chartable

    XVector=induced_band_representation.(Ref(sitegroup), group(first(lgirs[k])) , ρ, Ref(kvec), timereversal)

    return round.(character_table\XVector)

end

function return_irvec(sgnum::Int64, ρ::Int64, sitegroup::SiteGroup{D}, maximal::Bool, timereversal::Bool) where D
    maximal ? klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_maximal_lgs(sgnum, D)) : klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_littlegroups(sgnum, D))
    irvec = Int[]
    for (klab, kv) in klab2kv
        kv_array   = subduce_onto_littlegroups(klab, kv, ρ, sitegroup, timereversal)
        append!(irvec, Int.(kv_array))
    end
    return irvec
end

function return_irvec2d(sgnum::Integer, D::Integer, maximal::Bool, timereversal::Bool)
    sg=spacegroup(sgnum, D)
    wps=get_wycks(sgnum, D)
    A=Vector{Vector{Int}}(undef, 0)

    sitegroups=SiteGroup.(Ref(sg), wps)

    maximal_wyckoffs=Array{WyckPos{2}, 1}(undef, 0)

    for maximal_sitegroup in findmaximal(sitegroups)
        push!(maximal_wyckoffs, maximal_sitegroup.wp)
    end

    for wp in maximal_wyckoffs

        sitegroup=SiteGroup(sg, wp)
        ρ_max=size(site_symmetry_irreps(sitegroup, timereversal))[1]
        for ρ in 1:ρ_max
            push!(A, return_irvec(sgnum, ρ, sitegroup, maximal, timereversal))
        end
    end
    return hcat(A...)
end


function get_maximal_lgirreps(sgnum::Integer, D::Integer, timereversal::Bool)
    maximal_lgirreps=Dict{String, Array{LGIrrep{D}, 1}}()
    for (klab, lgirrep) in get_lgirreps(sgnum, D, timereversal)
        isspecial(kvec(lgirrep[1].lg)) || continue
            push!(maximal_lgirreps, klab => lgirrep)
    end
    return maximal_lgirreps
end

function get_maximal_lgs(sgnum::Integer, D::Integer)
    maximal_lgs=Dict{String, LittleGroup{D}}()
    for (klab, lg) in get_littlegroups(sgnum, D)
        isspecial(kvec(lg)) || continue
            push!(maximal_lgs, klab => lg)
    end
    return maximal_lgs
end

function make_bandrep_set(sgnum::Integer; D::Integer=2, maximal::Bool, timereversal::Bool)

    maximal ? klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_maximal_lgs(sgnum, D)) : klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_littlegroups(sgnum, D))

    irvec2d=return_irvec2d(sgnum, D, maximal, timereversal)
    klabs=Array{String, 1}(undef, 0)
    kvs=Array{KVec{D}, 1}(undef, 0)
    for (klab, kv) in klab2kv
        push!(klabs, klab)
        push!(kvs, kv)
    end

    spinful, allpaths, timeinvar=false, false, true
    
    lgirreps= maximal ? get_maximal_lgirreps(sgnum, D, timereversal) : get_lgirreps(sgnum, D, timereversal)

    irlabs=[([[lgg.cdml for lgg in lg] for (klab, lg) in lgirreps]...)...]
    irdims=[([[size(lgg.matrices[1])[1] for lgg in lg] for (klab, lg) in lgirreps]...)...]

    wps=get_wycks(sgnum, D)
    sg=spacegroup(sgnum, D)

    wps_special=String[]
    sitesyms=String[]
    bandreps_array=Array{BandRep, 1}(undef, 0)

    vec_index=1

    first_k=irlabs[1][1]

    dim_array=[i[1]==first_k for i in irlabs]

    dim_array=dim_array.*irdims
    
    sitegroups=SiteGroup.(Ref(sg), wps)

    maximal_wyckoffs=Array{WyckPos{2}, 1}(undef, 0)

    for maximal_sitegroup in findmaximal(sitegroups)
        push!(maximal_wyckoffs, maximal_sitegroup.wp)
    end

    for wp in maximal_wyckoffs
        sitegroup=SiteGroup(spacegroup(sgnum, D), wp)

        sitesym=_find_isomorphic_parent_pointgroup(sitegroup).label

        sitegroup_symmetry_irreps=site_symmetry_irreps(sitegroup, timereversal)
        ρ_max=size(sitegroup_symmetry_irreps)[1]
        for ρ in 1:ρ_max
            mulliken_label=mulliken(sitegroup_symmetry_irreps[ρ])
            dim=dot(irvec2d[:, vec_index], dim_array)
            band_rep=BandRep(string( wp.mult, wp.letter), sitesym, string( wp.mult, wp.letter) , dim, false, false, irvec2d[:, vec_index], irlabs)
            push!(wps_special, string( wp.mult))
            push!(sitesyms, sitesym)
            push!(bandreps_array, band_rep )
            vec_index+=1
        end
    end
    return BandRepSet(sgnum, bandreps_array, kvs, klabs, irlabs, allpaths, spinful, timeinvar)
end

#Below, we make an array of bandrepsets in 2 dimensions for maximal k points
BRS_VEC_ELEMENTARY_MAXIMAL = make_bandrep_set.(1:17, D=2, maximal=true, timereversal=false);

#Below, we make an array of bandrepsets in 2 dimensions for all k points
BRS_VEC_ELEMENTARY_ALL = make_bandrep_set.(1:17, D=2, maximal=false, timereversal=false);

#Below, we make an array of bandrepsets in 2 dimensions for maximal k points
BRS_VEC_ELEMENTARYTR_MAXIMAL= make_bandrep_set.(1:17, D=2, maximal=true, timereversal=true);

#Below, we make an array of bandrepsets in 2 dimensions for all k points
BRS_VEC_ELEMENTARYTR_ALL=make_bandrep_set.(1:17, D=2, maximal=false, timereversal=true);



