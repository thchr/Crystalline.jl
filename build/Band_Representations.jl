using LinearAlgebra
using Crystalline
using Markdown
using StaticArrays
#inv for symmetry operations
#xyzt(group(get_lgirreps(17, 2)["K"][3])[3])
#seitz(group(site_symmetry_irreps(siteg2b)[3])[1])

"Return the representative of the little group irrep representative"
function group_representative(lgirrep::SymOperation, sg::SpaceGroup)
    lgxyzt=xyzt(lgirrep)
    for i in 1:size(sg)[1]
        xyztiter=xyzt(sg[i])
        if xyztiter==lgxyzt
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


#The function below 
function in_site_symmetry(A::SiteGroup, b::SymOperation)
    wp=A.wp;
    wp_rvec=wp.qv.cnst
    in_symmetry=true;
    try
        Int.(rotation(b)*wp_rvec-wp_rvec);
    catch
        in_symmetry=false;
    end
    #print(in_symmetry)

    if in_symmetry==false
        return 0
    else 
        return findfirst(x -> x==rotation(b), rotation.(A))
    end
end

"This function gives the label for spacegroup(i, j)"
function get_label(i::Int64, j::Int64)
    return (label∘spacegroup)(i, j) 
end

"This function gives the site_symmetry irreps for the Site Symmetry Group"
function site_symmetry_irreps(sitegroup::SiteGroup{D}) where D
   # site_parent_pointgroup=Crystalline.find_parent_pointgroup(sitegroup)
    site_parent_pointgroup=_find_isomorphic_parent_pointgroup(sitegroup)

    get_pgirreps(site_parent_pointgroup.label, Val(D))
end

"The function below gives the character table for the site symmetry group provided"
function site_symmetry_irreps_characters(sitegroup::SiteGroup{D}) where D
    irreps=site_symmetry_irreps(sitegroup)
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
function induced_band_representation(sitegroup::SiteGroup{D}, h::SymOperation{D}, ρ::Int64, k::KVec) where D
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
        #print("here1\n\n\n")    

        gα′=gα[α′]
        tα′α′=rotation(h)*orbit_of_sitegroup[α′].qv.cnst+translation(h)-orbit_of_sitegroup[α′].qv.cnst
        symop_t=SymOperation(one(Crystalline.SquareStaticMatrices.SqSMatrix{D,Float64}), -tα′α′ )
        
      #  print("here2\n\n\n")    

        #println(tα′α′)
        #println(transformed_k)
        #println(exp(1im*dot(transformed_k,tα′α′)))
        #find g_alpha_inverse

      # for possible_inverse in 1:length_of_orbit 
       #     if  gα[α′] ∘ gα[possible_inverse]== spacegroup(17, D)[1] 
       #        gα_inverse=gα[possible_inverse];
       #    end
       #end

        gα_inverse=inv(gα′)    

        #site_symmetry_index=in_site_symmetry(sitegroup, gα_inverse ∘ symop_t ∘ h ∘ gα[α′])
       site_symmetry_index=in_site_symmetry(sitegroup, compose(gα_inverse, compose(symop_t, compose(h,  gα[α′], false), false ), false)   )
      
     #  print("here3\n\n\n")    

       #print(site_symmetry_index, "\n\n\n")

        if site_symmetry_index >0
            #println(tα′α′)
            #print(exp(2*pi*1im*dot(transformed_k,tα′α′)))
            χ_G_k_G=χ_G_k_G+exp(2im*π*dot(transformed_k, tα′α′))*tr(site_symmetry_irreps(sitegroup)[ρ].matrices[site_symmetry_index])
            #print(transformed_k)
            #print(dot(transformed_k, tα′α′))
        end
    end
    return χ_G_k_G
end


function subduce_onto_littlegroups(k::String, kvec::KVec, ρ::Int64, sitegroup::SiteGroup{D}) where D
    sgnum=sitegroup.num
    num_lgirreps= size(get_lgirreps(sgnum, D)[k])[1]
    
    size_lgirreps=size(get_lgirreps(sgnum, D)[k][1].matrices)[1]

    character_table=Matrix{Complex{Float64}}(undef, size_lgirreps, num_lgirreps )

    for j in 1:num_lgirreps
        for i in 1:size_lgirreps
            character_table[i, j]=tr(get_lgirreps(sgnum, D)[k][j].matrices[i])
        end
    end
    #return character_table
    
    #Now we obtain the vector that we require: 
    XVector=Matrix{Complex{Float64}}(undef, (size_lgirreps, 1))

    for i in 1:size_lgirreps
        rep_index=group_representative(group(get_lgirreps(sgnum, D)[k][1])[i], spacegroup(sgnum, D)) 
        iter=induced_band_representation(sitegroup, spacegroup(sgnum, D)[rep_index], ρ, kvec)
       XVector[i, 1]= iter
    end

    #return round.(character_table\XVector)
    #return character_table, XVector
    return round.(character_table\XVector)

end

#StringtoKVec=Dict("Γ"=>KVec("0, 0"), "K"=>KVec("1/3, 1/3") , "M"=>KVec("1/2, 0"))

function return_irvec(sgnum::Int64, ρ::Int64, sitegroup::SiteGroup{D}) where D
    klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_littlegroups(sgnum, D))
    irvec = Int[]
    for (klab, kv) in klab2kv
        isspecial(kv) || continue
        kv_array   = subduce_onto_littlegroups(klab, kv, ρ, sitegroup)
        append!(irvec, Int.(kv_array))
    end
    return irvec
end

function return_irvec2d(sgnum::Integer, D::Integer=2)
    sg=spacegroup(sgnum, D);
    wps=get_wycks(sgnum, D)
    A=Vector{Vector{Int}}(undef, 0)

    sitegroups=SiteGroup.(Ref(sg), wps)

    maximal_wyckoffs=Array{WyckPos{2}, 1}(undef, 0)

    for maximal_sitegroup in findmaximal(sitegroups)
        push!(maximal_wyckoffs, maximal_sitegroup.wp)
    end

    #for wp in wps
    for wp in maximal_wyckoffs
       # print(wp)
        #isspecial(qvec(wp)) || continue 
        #isspecial(wp) || continue 

        sitegroup=SiteGroup(sg, wp)
        ρ_max=size(site_symmetry_irreps(sitegroup))[1]
        #print(ρ_max)
        for ρ in 1:ρ_max
            push!(A, return_irvec(sgnum, ρ, sitegroup))
        end
    end
    return hcat(A...)
end


function get_maximal_lgirreps(sgnum::Integer, D::Integer)
    maximal_lgirreps=Dict{String, Array{LGIrrep{D}, 1}}()
    for (klab, lgirrep) in get_lgirreps(sgnum, D)
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


function make_bandrep_set(sgnum::Integer; D::Integer=2)
    #klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_littlegroups(sgnum, D))
    klab2kv = Dict(klab => kvec(lg) for (klab, lg) in get_maximal_lgs(sgnum, D))

    irvec2d=return_irvec2d(sgnum, 2)
    klabs=Array{String, 1}(undef, 0)
    kvs=Array{KVec{D}, 1}(undef, 0)
    for (klab, kv) in klab2kv
        push!(klabs, klab)
        push!(kvs, kv)
       # if isspecial(kv)==true 
      #  end
    end

    spinful, allpaths, timeinvar=false, false, true
    
    #irlabs=[([[lgg.cdml for lgg in lg] for (klab, lg) in get_lgirreps(sgnum, D)]...)...]
    #irdims=[([[size(lgg.matrices[1])[1] for lgg in lg] for (klab, lg) in get_lgirreps(sgnum, D)]...)...]

    irlabs=[([[lgg.cdml for lgg in lg] for (klab, lg) in get_maximal_lgirreps(sgnum, D)]...)...]
    irdims=[([[size(lgg.matrices[1])[1] for lgg in lg] for (klab, lg) in get_maximal_lgirreps(sgnum, D)]...)...]


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

            sitegroup_symmetry_irreps=site_symmetry_irreps(sitegroup)
            ρ_max=size(sitegroup_symmetry_irreps)[1]
            for ρ in 1:ρ_max
                mulliken_label=mulliken(sitegroup_symmetry_irreps[ρ])

               # print(mulliken_label)
                dim=dot(irvec2d[:, vec_index], dim_array)
                 band_rep=BandRep(string( wp.mult, wp.letter), sitesym, string( wp.mult, wp.letter) , dim, false, false, irvec2d[:, vec_index], irlabs)
                 push!(wps_special, string( wp.mult))
                 push!(sitesyms, sitesym)
                 push!(bandreps_array, band_rep )
                 #print(band_rep)
                 vec_index+=1
        end
    end
    #bandreps=Array{BandRep, 1}(undef, 0)
   # print(bandreps_array)

    return BandRepSet(sgnum, bandreps_array, kvs, klabs, irlabs, allpaths, spinful, timeinvar)

end

BRS_VEC=[];

BRS_NUMS=[];

for br_num in 1:17
    if br_num ∉[4, 7, 8, 12]
        push!(BRS_VEC, make_bandrep_set(br_num, D=2))
        push!(BRS_NUMS, br_num )
    end
end


