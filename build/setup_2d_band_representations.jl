using LinearAlgebra
using Crystalline
using Crystalline: irdim, constant, AbstractIrrep, AbstractGroup, iscorep, _mulliken
import Crystalline: label, mulliken, realify, group
using StaticArrays

# The implementation here follows Elcoro et al., Phys. Rev. B 97, 035139 (2018)
# (https://doi.org/10.1103/PhysRevB.97.035139), specifically, Sections II.C-D

# ---------------------------------------------------------------------------------------- #
# Site symmetry irreps

struct SiteIrrep{D} <: AbstractIrrep{D}
    cdml     :: String
    g        :: SiteGroup{D}
    matrices :: Vector{Matrix{ComplexF64}}
    reality  :: Reality
    iscorep  :: Bool
    pglab    :: String
end


function _find_isomorphic_parent_pointgroup(g::AbstractGroup{D}) where D
    ctᴳ = MultTable(g).table
    @inbounds for iuclab in PGS_IUCs[D]
        pg  = pointgroup(iuclab, Val(D))
        ctᴾ = MultTable(pg).table
        if ctᴳ == ctᴾ # bit sloppy; would miss ismorphisms that are "concealed" by row/column swaps
            return pg
        end
    end
    return nothing # in case we didn't find any isomorphic parent
end

"""
    get_siteirreps(sitegroup::SiteGroup) --> Vector{PGIrrep}

Return the site symmetry irreps associated with the provided `SiteGroup`, derived from a
search over isomorphic point groups.
"""
function get_siteirreps(siteg::SiteGroup{D}) where D
    parent_pg = _find_isomorphic_parent_pointgroup(siteg)
    pglab     = label(parent_pg)
    pgirs     = get_pgirreps(pglab, Val(D))
    
    siteirs = map(pgirs) do pgir
        SiteIrrep{D}(label(pgir), siteg, pgir.matrices, reality(pgir), pgir.iscorep, pglab)
    end
    return siteirs
end
pglabel(siteir::SiteIrrep)  = siteir.pglab # associated point group label
mulliken(siteir::SiteIrrep) = _mulliken(pglabel(siteir), label(siteir), iscorep(siteir))

# ---------------------------------------------------------------------------------------- #
# Misc utility functions 

function realify(lgirsd::Dict{<:Any, <:Vector{<:LGIrrep}})
    return Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd)
end

"""
    reduce_dict_of_vectors([f=identity,] d::Dict{_, <:Vector})  -->  Vector{T}

Return the concatenated vector of all of vectors in `d` under the element-wise application
of `f`. Effectively flattens a `Dict` of `Vector`s to a single `Vector`.

Note that application of `f` to vector-elements of `d` must return a stable type `T`.
"""
function reduce_dict_of_vectors(f::F, d::Dict{<:Any, <:Vector}) where F
    # get element type of application of `f` to elements of vectors in `d`; assumed fixed
    eltype = typeof(f(first(first(values(d)))))
    # get total number of elements across vectors in `d`
    N = sum(((_, v),) -> length(v), d)
    # preallocate output vector
    w      = Vector{eltype}(undef, N)
    start  = 1
    for v in values(d)
        stop = start + length(v) - 1
        @inbounds for (i, j) in enumerate(start:stop)
            w[j] = f(v[i])
        end
        start = stop + 1
    end
    return w
end
reduce_dict_of_vectors(d::Dict{<:Any, <:Vector}) = reduce_dict_of_vectors(identity, d)

# ---------------------------------------------------------------------------------------- #
# Bandrep related functions: induction/subduction

"""
    induced_band_representation(siteir, h, kv)

Return the band representation induced by the provided `SiteIrrep` evaluated at `kv` and for
a `SymOperation` `h`.
"""
function induce_bandrep(siteir::SiteIrrep{D}, h::SymOperation{D}, kv::KVec{D}) where D

    siteg  = group(siteir)
    wp     = wyck(siteg)   
    orbits = orbit(siteg, wp)
    gαs    = cosets(siteg)
    kv′    = constant(h∘kv) # <-- TODO: Why only constant part?
    χs     = characters(siteir)
    
    # sum over all the wyckoff positions/cosets in the orbit 
    χᴳₖ = zero(ComplexF64)
    for (wpα′, gα′) in zip(orbits, gαs)
        tα′α′ = constant(vec(h∘wpα′) - vec(wpα′)) # TODO: <-- explain why we only need constant part here?
        opᵗ   = SymOperation(-tα′α′)

        gα′⁻¹     = inv(gα′)
        gα′⁻¹ggα′ = compose(gα′⁻¹, compose(opᵗ, compose(h, gα′, false), false ), false)

        site_symmetry_index = findfirst(≈(gα′⁻¹ggα′), siteg)
        if site_symmetry_index !== nothing
            χᴳₖ += cis(2π*dot(kv′, tα′α′)) * χs[site_symmetry_index]
            # TODO: The sign in this `cis(...)` above is different from in Elcoro's. Why?
        end
    end
    return χᴳₖ
end

function subduce_onto_lgirreps(siteir::SiteIrrep{D}, lgirs::Vector{LGIrrep{D}}) where D
    lg = group(first(lgirs))
    kv = kvec(lg)

    # characters of induced site representation and little group irreps
    site_χs  = induce_bandrep.(Ref(siteir), lg, Ref(kv))
    lgirs_χm = CharacterTable(lgirs).chartable

    # little group irrep multiplicities, after subduction
    m  = lgirs_χm\site_χs
    m′ = round.(Int, real.(m)) # truncate to integers
    isapprox(m, m′, atol=1e-12) || error("failure to convert to integers")
    return m′
end

function calc_bandrep(siteir::SiteIrrep{D}, lgirsd::Dict{String, Vector{LGIrrep{D}}};
            irlabs::Vector{String} = reduce_dict_of_vectors(label, lgirsd),
            irdims::Vector{Int}    = reduce_dict_of_vectors(irdim, lgirsd)) where D

    irvec = Int[]
    for (_, lgirs) in lgirsd
        kv_array = subduce_onto_lgirreps(siteir, lgirs)
        append!(irvec, Int.(kv_array))
    end    
    
    irdims_Γmask = [irlab[1]=='Γ' for irlab in irlabs] .* irdims # masked w/ Γ-irreps only
    brdim   = dot(irvec, irdims_Γmask)
    wycklab = label(wyck(group(siteir)))
    spinful      = false # NB: default; Crystalline currently doesn't have spinful irreps
    decomposable = false # NB: placeholder, because we don't know the true answer presently

    return BandRep(wycklab, pglabel(siteir), mulliken(siteir)*"↑G", brdim, decomposable,
                   spinful, irvec, irlabs)
end

# ---------------------------------------------------------------------------------------- #
# TODO: Doc-string
function calc_bandreps(sgnum::Integer, Dᵛ::Val{D}=Val(2);
                       allpaths::Bool=false, timereversal::Bool=true) where D

    # get all the little group irreps that we want to subduce onto
    lgirsd = get_lgirreps(sgnum, Val(D))
    allpaths     || filter!(((_, lgirs),) -> isspecial(first(lgirs)), lgirsd)
    timereversal && (lgirsd = realify(lgirsd))

    irlabs = reduce_dict_of_vectors(label, lgirsd)
    irdims = reduce_dict_of_vectors(irdim, lgirsd)

    # get the bandreps induced by every maximal site symmetry irrep
    wps    = get_wycks(sgnum, Dᵛ)
    sg     = spacegroup(sgnum, Dᵛ)
    sitegs = findmaximal(SiteGroup.(Ref(sg), wps))
    brs    = BandRep[]
    for siteg in sitegs
        siteirs = get_siteirreps(siteg)
        timereversal && (siteirs = realify(siteirs))
        append!(brs, calc_bandrep.(siteirs, Ref(lgirsd); irlabs=irlabs, irdims=irdims))
    end

    # TODO: There's potentially some kind of bad implicit dependence/overreliance on fixed
    #       on iteration order of `Dict`s here and related methods; probably OK, but unwise
    klabs  = collect(keys(lgirsd))
    kvs    = kvec.(first.(values(lgirsd)))
    spinful      = false # NB: default; Crystalline currently doesn't have spinful irreps
    decomposable = false # NB: placeholder, because we don't know the true answer presently

    return BandRepSet(sgnum, brs, kvs, klabs, irlabs, allpaths, spinful, timereversal)
end