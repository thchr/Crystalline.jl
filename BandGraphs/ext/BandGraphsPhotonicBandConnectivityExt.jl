module BandGraphsPhotonicBandConnectivityExt

# ---------------------------------------------------------------------------------------- #

if isdefined(Base, :get_extension)
    using PhotonicBandConnectivity
else
    using ..PhotonicBandConnectivity
end
const PBC = PhotonicBandConnectivity
using PhotonicBandConnectivity: is_vrep

using BandGraphs: SubductionTable
using Crystalline: lgirreps, label, klabel, realify!, 
                   LGIrrep, REAL, find_representation, group,
                   SymmetryVector, irreps, multiplicities, occupation,
                   Collection
using LinearAlgebra: I

import BandGraphs: add_transverse_vrep!, exfiltrate_transverse_vrep

# ---------------------------------------------------------------------------------------- #

function add_transverse_vrep!(
    subts::AbstractVector{SubductionTable{D}},
    vrep::LGIrrep{D} # virtual, transverse irrep
    ) where D
    add_transverse_vrep!(subts, label(vrep))
end
function add_transverse_vrep!(
    subts::AbstractVector{SubductionTable{D}},
    vrep_lab::AbstractString
    ) where D
    # include the transverse & singular gamma irrep in the subduction tables at Γ
    irlabs, coefs = _parse_composite_rep(vrep_lab)
    for (i, subt) in enumerate(subts)
        kᴳlab = string(subt.c.kᴳ.label)
        if kᴳlab == "Γ" || kᴳlab == "Γ′"
            row = zeros(Int, size(subt.table, 2))
            for (c, irlab) in zip(coefs, irlabs)
                irlab_monodromy = replace(irlab, "Γ"=>"Γ′")
                idx = findfirst(irlabᴳ->irlabᴳ==irlab || irlabᴳ==irlab_monodromy, subt.irlabsᴳ)
                if isnothing(idx)
                    error(lazy"could not find $irlab in subduction table; contains $(subt.irlabsᴳ)")
                end
                row .+= c * subt.table[idx, :]
            end
            @assert all(≥(0), row) # we cannot operate with negative subduction counts
            table′ = vcat(subt.table, permutedims(row))
            vrep_lab′ = kᴳlab == "Γ′" ? replace(vrep_lab, "Γ"=>"Γ′") : vrep_lab
            irlabsᴳ′ = vcat(subt.irlabsᴳ, vrep_lab′)
            subt′ = SubductionTable{D}(subt.num, subt.c, irlabsᴳ′, subt.irlabsᴴ, table′, 
                                       subt.monodromy)
            subts[i] = subt′
        end
    end
    return subts
end

function _parse_composite_rep(Γlab)
    i = 1
    coefs = Int[]
    irlabs = String[]
    Γlab = strip(Γlab, ('(', ')'))
    while true
        i′ = findnext('Γ', Γlab, i)
        isnothing(i′) && error("could not find next irrep entry; ill-formed composite label $Γlab")
        coef_str = Γlab[i:prevind(Γlab, i′)]
        j = findnext(r"\+|\-", Γlab, i′)
        if isnothing(j)
            j′ = lastindex(Γlab)
        else
            j′ = prevind(Γlab, only(j))
        end
        coef = if isempty(coef_str) || coef_str == "+"
            1
        elseif coef_str == "-"
            -1
        else
            parse(Int, coef_str)
        end
        push!(coefs, coef)
        push!(irlabs, Γlab[i′:j′])
        isnothing(j) && break
        i = only(j)
    end
    return irlabs, coefs
end

# ---------------------------------------------------------------------------------------- #

# split the composite vrep up, removing it from the symmetry vector `n`, and adding the
# multiplicities of its constituent irreps to `n`
function exfiltrate_transverse_vrep(n::SymmetryVector)
    Γidx = something(findfirst(lgirs -> klabel(first(lgirs)) == "Γ", irreps(n)))
    iridx = findfirst(lgir -> is_vrep(lgir), irreps(n)[Γidx])
    isnothing(iridx) && error(lazy"provided symmetry vector $n does not contain a vrep")
    vrep = irreps(n)[Γidx][iridx]
    exfiltrate_transverse_vrep(n, label(vrep), Γidx, iridx)
end
function exfiltrate_transverse_vrep(
    n::SymmetryVector,
    vrep_lab::AbstractString,
    Γidx::Int,
    iridx::Int
    )
    irlabs, coefs = _parse_composite_rep(vrep_lab)

    lgirsΓ = irreps(n)[Γidx]
    multsΓ = multiplicities(n)[Γidx]
    multsΓ′ = copy(multsΓ)
    for (i, lgir) in enumerate(lgirsΓ)
        i′ = findfirst(==(label(lgir)), irlabs)
        isnothing(i′) && continue
        multsΓ′[i] += coefs[i′]
    end
    
    lgirsΓ′ = Collection(copy(parent(lgirsΓ)))
    deleteat!(parent(lgirsΓ′), iridx)
    deleteat!(multsΓ′, iridx)
    mults′ = [i == Γidx ? multsΓ′ : mults for (i, mults) in enumerate(multiplicities(n))]
    lgirsv′ = [i == Γidx ? lgirsΓ′ : lgirs for (i, lgirs) in enumerate(irreps(n))]

    return SymmetryVector(lgirsv′, mults′, occupation(n))
end

end # module 