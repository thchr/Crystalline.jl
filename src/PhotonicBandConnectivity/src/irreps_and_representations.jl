function compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    # Find the Hilbert basis that respects the compatibility relations
    sb, _, BRS = compatibility_bases(sgnum, spinful=false, timereversal=timereversal)
    nsᴴ = matrix(sb)
    # Find the indices of the Γ irreps in `BRS::BandRepSet` (and hence in `nsᴴ`), and how  
    # they map to the corresponding irrep indices in `lgirs`
    irlabs_brs = irreplabels(BRS)
    irlabs_lgirs = Crystalline.formatirreplabel.(label.(lgirs))
    Γidxs = map(irlab->findfirst(==(irlab), irlabs_brs), irlabs_lgirs)

    return sb, Γidxs
end

function get_lgirreps_at_Γ(sgnum::Integer, Dᵛ::Val=Val(3)) # small irreps at Γ
   lgirs = first(get_lgirreps(sgnum,  Dᵛ))
   kv = kvec(first(lgirs))
   @assert all(iszero, kv.k₀) && isspecial(kv) # Make sure that lgirs indeed is sampled at Γ

   return lgirs
end

# irrep-expansions/representation at Γ for the transverse (2T), longitudinal (1L), and triad
# (2T+1L) plane wave branches that touch ω=0 at Γ
"""
    find_representation²ᵀ⁺¹ᴸ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation²ᵀ⁺¹ᴸ(sgnum::Integer; timereversal::Bool=false)
"""
function find_representation²ᵀ⁺¹ᴸ end
"""
    find_representation¹ᴸ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation¹ᴸ(sgnum::Integer; timereversal::Bool=false)
"""
function find_representation¹ᴸ    end
"""
    find_representation²ᵀ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation²ᵀ(sgnum::Integer; timereversal::Bool=false)
"""
function find_representation²ᵀ    end

for postfix in ("²ᵀ⁺¹ᴸ", "¹ᴸ", "²ᵀ")
    f = Symbol("find_representation"*postfix) # method to be defined
    symvals_fun = Symbol("get_symvals"*postfix)

    # "root" accessors via lgirs
    @eval function $f(lgirs::AbstractVector{LGIrrep{3}})
        lg = group(first(lgirs))
        symvals = $symvals_fun(lg)

        return find_representation(symvals, lgirs)
    end

    # convenience accessors via 
    @eval function $f(sgnum::Integer; timereversal::Bool=false)
        lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
        timereversal && (lgirs = realify(lgirs))

        return $f(lgirs)
    end
end