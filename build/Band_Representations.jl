using LinearAlgebra
using Crystalline
using StaticArrays
import Crystalline:(get_lgirreps)

function get_lgirreps(sgnum::Int64, D::Int64, timereversal::Bool)
    if timereversal==false
        get_lgirreps(sgnum, D) 
    else 
        tr_lgirreps=Dict{String, Array{LGIrrep{D}, 1}}()
        for (klab, lg) in get_lgirreps(sgnum, D)
            push!(tr_lgirreps, klab => realify(lg))
        end
        return tr_lgirreps
    end
end

function group_representative(symop::SymOperation, sg::SpaceGroup)
    symop_xyzt=xyzt(symop)
    for i in 1:size(sg)[1]
        xyztiter=xyzt(sg[i])
        if xyztiter==symop_xyzt
            return i
        end
    end
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

function in_site_symmetry(sitegroup::SiteGroup, symop::SymOperation)
    #index= findfirst(x -> x==rotation(symop), rotation.(sitegroup)) 
    index= findfirst(x -> x==xyzt(symop), xyzt.(sitegroup)) 
    if index == nothing
        return 0 
    else 
        return index
    end
end

"This function gives the label for spacegroup(i, j)"
function get_label(i::Int64, j::Int64)
    return (label∘spacegroup)(i, j) 
end

"This function gives the site_symmetry irreps for the Site Symmetry Group"
function site_symmetry_irreps(sitegroup::SiteGroup{D}, timereversal::Bool) where D
    site_parent_pointgroup=_find_isomorphic_parent_pointgroup(sitegroup)

    timereversal==false ? get_pgirreps(site_parent_pointgroup.label, Val(D)) : realify(get_pgirreps(site_parent_pointgroup.label, Val(D)))
end

"The function below gives the character table for the site symmetry group provided"
function site_symmetry_irreps_characters(sitegroup::SiteGroup{D}, timereversal::Bool) where D
    irreps=site_symmetry_irreps(sitegroup, timereversal)
    character_table=Matrix{Complex{Float64}}(undef, size(irreps[1].matrices)[1], size(irreps)[1] )
    size(character_table)
    for j in 1:size(irreps)[1]
        for i in 1:size(irreps[1].matrices)[1]
            character_table[i, j]=tr(irreps[j].matrices[i]) +0.00im +0.00
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
    for α′ in 1:length_of_orbit 

        gα′=gα[α′]
        tα′α′=rotation(h)*orbit_of_sitegroup[α′].qv.cnst+translation(h)-orbit_of_sitegroup[α′].qv.cnst
        symop_t=SymOperation(one(Crystalline.SquareStaticMatrices.SqSMatrix{D,Float64}), -tα′α′ )

        gα_inverse=inv(gα′)    

        site_symmetry_index=in_site_symmetry(sitegroup, compose(gα_inverse, compose(symop_t, compose(h,  gα[α′], false), false ), false)   )
      
        if site_symmetry_index >0
            χ_G_k_G=χ_G_k_G+exp(2im*π*dot(transformed_k, tα′α′))*tr(site_symmetry_irreps(sitegroup, timereversal)[ρ].matrices[site_symmetry_index])
        end
    end
    return χ_G_k_G
end

function subduce_onto_littlegroups(k::String, kvec::KVec, ρ::Int64, sitegroup::SiteGroup{D}, timereversal::Bool) where D
    sgnum=sitegroup.num
    num_lgirreps= size(get_lgirreps(sgnum, D, timereversal)[k])[1]
    
    size_lgirreps=size(get_lgirreps(sgnum, D, timereversal)[k][1].matrices)[1]

    character_table=Matrix{Complex{Float64}}(undef, size_lgirreps, num_lgirreps )

    for j in 1:num_lgirreps
        for i in 1:size_lgirreps
            character_table[i, j]=tr(get_lgirreps(sgnum, D, timereversal)[k][j].matrices[i])
        end
    end

    XVector=Matrix{Complex{Float64}}(undef, (size_lgirreps, 1))

    for i in 1:size_lgirreps
        rep_index=group_representative(group(get_lgirreps(sgnum, D, timereversal)[k][1])[i], spacegroup(sgnum, D)) 
        iter=induced_band_representation(sitegroup, spacegroup(sgnum, D)[rep_index], ρ, kvec, timereversal)
        XVector[i, 1]= iter
    end

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

    irvec2d=return_irvec2d(sgnum, 2, maximal, timereversal)
    klabs=Array{String, 1}(undef, 0)
    kvs=Array{KVec{D}, 1}(undef, 0)
    for (klab, kv) in klab2kv
        push!(klabs, klab)
        push!(kvs, kv)
    end

    spinful, allpaths, timeinvar=false, false, true
    
    maximal ? lgirreps=get_maximal_lgirreps(sgnum, D, timereversal) : lgirreps=get_lgirreps(sgnum, D, timereversal)

    irlabs=[([[lgg.cdml for lgg in lg] for (klab, lg) in lgirreps]...)...]
    irdims=[([[size(lgg.matrices[1])[1] for lgg in lg] for (klab, lg) in lgirreps]...)...]

    wps=get_wycks(sgnum, D)
    sg=spacegroup(sgnum, D)

    wps_special=[]
    sitesyms=[]
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

BRS_VEC_ELEMENTARY_MAXIMAL=BandRepSet[];
BRS_VEC_ELEMENTARY_ALL=BandRepSet[];

BRS_VEC_ELEMENTARYTR_MAXIMAL=BandRepSet[];
BRS_VEC_ELEMENTARYTR_ALL=BandRepSet[];

#Below, we make an array of bandrepsets in 2 dimensions for maximal k points
for br_num in 1:17
    push!(BRS_VEC_ELEMENTARY_MAXIMAL, make_bandrep_set(br_num, D=2, maximal=true, timereversal=false))
end

#Below, we make an array of bandrepsets in 2 dimensions for all k points
for br_num in 1:17
    push!(BRS_VEC_ELEMENTARY_ALL, make_bandrep_set(br_num, D=2, maximal=false, timereversal=false))
end

#Below, we make an array of bandrepsets in 2 dimensions for maximal k points
for br_num in 1:17
    push!(BRS_VEC_ELEMENTARYTR_MAXIMAL, make_bandrep_set(br_num, D=2, maximal=true, timereversal=true))
end

#Below, we make an array of bandrepsets in 2 dimensions for all k points
for br_num in 1:17
    push!(BRS_VEC_ELEMENTARYTR_ALL, make_bandrep_set(br_num, D=2, maximal=false, timereversal=true))
end


