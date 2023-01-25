#using AbstractTrees
using Crystalline
using Crystalline: PG_ORDERs, _is_setting_matched_subgroup, _subgroup_fastpath_checks

# TODO: maybe worth following approach in https://doi.org/10.1103/PhysRevB.31.2908

function find_minimal_supergroups(H::PointGroup{D}) where D
    Gs = PointGroup{D}[] # supergroups
    for iuclabᴳ in PG_IUCs[D]
        # check if H < G
        bool, G = issubgroup_thorough(iuclabᴳ, H; check_minimal=true)
        bool || continue
        G::PointGroup{D}
        
        # check if G < G′ for any G′ in Gs; otherwise, could be a minimal supergroup
        if !any(G′ -> issubgroup_thorough(G′, G), Gs)   
            push!(Gs, G)
        end
    end

    # check if we accidentally added a group to Gs that is a subgroup of another group in Gs
    kill_idxs = Int[]
    for (i, G) in enumerate(Gs)
        if any(G′ -> issubgroup_thorough(G′, G), @view Gs[1:i-1])  ||
           any(G′ -> issubgroup_thorough(G′, G), @view Gs[i+1:end])
           push!(kill_idxs, i)
        end
    end
    deleteat!(Gs, kill_idxs)

    return Gs
end

function issubgroup_thorough(G::PointGroup{D}, H::PointGroup{D}) where D
    Nᴳ, Nᴴ = order(G), order(H)
    _subgroup_fastpath_checks(Nᴳ, Nᴴ) || return false

    # must now check more carefully if `iuclab` is a minimal supergroup
    return issubgroup_thorough_slow(G, H)
end

function issubgroup_thorough(iuclabᴳ::AbstractString, H::PointGroup{D};
                             check_minimal::Bool=false) where D
    Nᴳ, Nᴴ = PG_ORDERs[iuclabᴳ], order(H)
    _subgroup_fastpath_checks(Nᴳ, Nᴴ) || return false, nothing

    if check_minimal && div(Nᴳ, Nᴴ) > D+1
        # minimal supergroup index of a _PointGroup_ is bounded to D+1 (cf. p. 795 in ITA5)
        return false, nothing
    end

    # must now check more carefully if `iuclab` is a minimal supergroup
    G = pointgroup(iuclabᴳ, Val(D)) # candidate supergroup
    return issubgroup_thorough_slow(G, H), G
end

function issubgroup_thorough_slow(G::PointGroup{D}, H::PointGroup{D}) where D
    # check if H < G, assuming that G and H are in the same setting
    _is_setting_matched_subgroup(G, H) && return true

    # check if H < G, allowing for setting differences

    # definitely proved that H is not a subgroup of G
    return false
end

# find the subgroups of `G` of index `i`
function subgroups(G::PointGroup, i::Int)
    
end

label.(find_minimal_supergroups(pointgroup("1")))