using LinearAlgebra
using Crystalline
using Crystalline: irdim, constant
import Crystalline: get_lgirreps
using StaticArrays

# The implementation here follows Elcoro et al., Phys. Rev. B 97, 035139 (2018)
# (https://doi.org/10.1103/PhysRevB.97.035139), specifically, Sections II.C-D

# FIXME: The code is unnecessarily slow because it calls `getters` (e.g. `get_lgirreps`,
#        `site_symmetry_irreps`, etc.) redundantly & repetitvely. Ought to restructure.

function get_lgirreps(sgnum::Int64, Dᵛ::Val, timereversal::Bool)
    if !timereversal
        return get_lgirreps(sgnum, Dᵛ)
    else 
        return Dict(klab => realify(lgirs) for (klab, lgirs) in get_lgirreps(sgnum, Dᵛ))
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

"""
    site_symmetry_irreps(sitegroup::SiteGroup, timereversal) --> Vector{PGIrrep}

Return the site symmetry irreps of the provided `SiteGroup` as a `Vector{PGIrrep}`.
"""
function site_symmetry_irreps(sitegroup::SiteGroup{D}, timereversal::Bool) where D
    site_parent_pointgroup = _find_isomorphic_parent_pointgroup(sitegroup)
    pgirs = get_pgirreps(site_parent_pointgroup.label, Val(D))

    return timereversal ? realify(pgirs) : pgirs
end

"""
    induced_band_representation(sitegroup, h, ρ, kv, timereversal)

Return the band representation induced by the provided `SiteGroup`'s `ρ`th irrep evaluated
at `kv` and for a `SymOperation` `h`.
"""
function induced_band_representation(sitegroup::SiteGroup{D}, h::SymOperation{D}, ρ::Int64, 
            kv::KVec, timereversal::Bool) where D

    wp     = wyck(sitegroup)   
    orbits = orbit(sitegroup, wp)
    gαs    = cosets(sitegroup)
    kv′    = constant(h∘kv) # <-- TODO: Why only constant part?
    χs     = characters(site_symmetry_irreps(sitegroup, timereversal)[ρ])
    
    # sum over all the wyckoff positions/cosets in the orbit 
    χ_G_k_G = zero(ComplexF64)
    for (wpα′, gα′) in zip(orbits, gαs)
        tα′α′ = constant(vec(h∘wpα′) - vec(wpα′)) # TODO: <-- explain why we only need constant part here?
        opᵗ   = SymOperation(-tα′α′)

        gα′⁻¹     = inv(gα′)
        gα′⁻¹ggα′ = compose(gα′⁻¹, compose(opᵗ, compose(h, gα′, false), false ), false)

        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), sitegroup)
        if site_symmetry_index !== nothing
            χ_G_k_G += cis(2π*dot(kv′, tα′α′)) * χs[site_symmetry_index]
            # TODO: The sign in this `cis(...)` above is different from in Elcoro's. Why?
        end
    end
    return χ_G_k_G
end

function subduce_onto_littlegroups(lg::LittleGroup{D}, ρ::Int64, sitegroup::SiteGroup{D},
            timereversal::Bool) where D

    sgnum = num(sitegroup)
    lgirs = get_lgirreps(sgnum, Val(D), timereversal)[klabel(lg)]
    χs    = CharacterTable(lgirs).chartable

    # TODO: Meaningful variable name!?
    XVector = induced_band_representation.(Ref(sitegroup), lg, ρ, Ref(kvec(lg)), timereversal)

    # TODO: Check correctness of this Float->Int conversion..!
    return round.(χs\XVector)
end

function return_irvec(sgnum::Int64, ρ::Int64, sitegroup::SiteGroup{D}, allpaths::Bool,
            timereversal::Bool) where D

    lgs = get_littlegroups(sgnum, Val(D))
    allpaths || filter!(((_, lg),) -> isspecial(kvec(lg)), lgs)
    
    # TODO: This is not safe: iteration-order of dicts is not a stable thing. Fix.
    irvec = Int[]
    for (klab, lg) in lgs
        kv_array = subduce_onto_littlegroups(lg, ρ, sitegroup, timereversal)
        append!(irvec, Int.(kv_array))
    end
    return irvec
end

function return_irvec2d(sgnum::Integer, Dᵛ::Val{D}, allpaths::Bool, timereversal::Bool) where D
    sg  = spacegroup(sgnum, Dᵛ)
    wps = get_wycks(sgnum, Dᵛ)

    sitegroups = findmaximal(SiteGroup.(Ref(sg), wps))
    A = Vector{Int}[]
    for sitegroup in sitegroups
        Nⁱʳˢ = length(site_symmetry_irreps(sitegroup, timereversal))
        for ρ in 1:Nⁱʳˢ
            push!(A, return_irvec(sgnum, ρ, sitegroup, allpaths, timereversal))
        end
    end
    return hcat(A...)
end

# ---------------------------------------------------------------------------------------- #
# TODO: Doc-string
function make_bandrep_set(sgnum::Integer, Dᵛ::Val{D}=Val(2); allpaths::Bool=true, 
            timereversal::Bool=true) where D

    lgirsd = get_lgirreps(sgnum, Dᵛ, timereversal)
    allpaths || filter!(((_, lgirs),) -> isspecial(kvec(first(lgirs))), lgirsd)
    lgs = Dict(klab => group(first(lgirs)) for (klab, lgirs) in lgirsd)

    irlabs = vcat((label.(lgirs) for lgirs in values(lgirsd))...)
    irdims = vcat((irdim.(lgirs) for lgirs in values(lgirsd))...)

    wps = get_wycks(sgnum, Dᵛ)
    sg  = spacegroup(sgnum, Dᵛ)

    irdims_Γmask = [irlab[1]=='Γ' for irlab in irlabs] .* irdims # masked w/ Γ-irreps only

    sitegroups   = findmaximal(SiteGroup.(Ref(sg), wps))
    spinful      = false
    decomposable = false # NB: placeholder, because we don't know the true answer presently
    vec_idx      = 1

    irvec2d = return_irvec2d(sgnum, Dᵛ, allpaths, timereversal)
    brs = BandRep[]
    for sitegroup in sitegroups
        sitesym = label(_find_isomorphic_parent_pointgroup(sitegroup))
        wplab   = label(wyck(sitegroup))

        sitegroup_irs = site_symmetry_irreps(sitegroup, timereversal)
        for sitegroup_ir in sitegroup_irs
            irvec = irvec2d[:, vec_idx]
            mulliken_label = mulliken(sitegroup_ir)
            brdim = dot(irvec, irdims_Γmask)
            br    = BandRep(wplab, sitesym, mulliken_label*"↑G", brdim, decomposable, 
                            spinful, irvec, irlabs)
            push!(brs, br)
            vec_idx += 1
        end
    end

    # TODO: There's potentially some kind of bad implicit dependence on iteration order of
    #       `Dict`s here
    #       --> Fix.
    klabs = collect(keys(lgs))
    kvs   = kvec.(values(lgs))

    return BandRepSet(sgnum, brs, kvs, klabs, irlabs, allpaths, spinful, timereversal)
end