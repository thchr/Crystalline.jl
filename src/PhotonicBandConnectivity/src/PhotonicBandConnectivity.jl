module PhotonicBandConnectivity

using Crystalline, LinearAlgebra 
import Crystalline: rotation

include(pkgdir(Crystalline)*"/src/SymmetryBases/SymmetryBases.jl")
using .SymmetryBases

export minimal_expansion_of_zero_freq_bands, sum_symbases, sum_symbases!

# ------------------------------------------------------------------------------------------

include("planewave_symvals.jl")
include("irreps_and_representations.jl")
include("symbasis_utils.jl")
include("constrained_expansions.jl")
#include("src/legacy_constrained_expansions.jl")

# ------------------------------------------------------------------------------------------

function minimal_expansion_of_zero_freq_bands(sgnum::Integer; 
                                              timereversal::Bool=true, verbose::Bool=true,
                                              shuffle_1Lpick::Bool=false)

    # Irreps at Γ, irrep-multiplicities of ω=0 2T bands, and symmetry operations
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    timereversal && (lgirs = realify(lgirs))
    lg = group(first(lgirs))
    rotvals = map(op->(W=rotation(op); Crystalline.rotation_order_3d(det(W), tr(W))), lg)

    # 2T irreps; check if "simple treatment"/fast-path is applicable
    ms²ᵀ = find_representation²ᵀ(lgirs)
    has_nonmirror_improper = any(∈((-1, -3, -4, -6)), rotvals)
    is_regular²ᵀ = all(≥(0), ms²ᵀ)
    if !has_nonmirror_improper && is_regular²ᵀ 
        # Scenario (easy case): All symvals known & regular 2T irrep

        # If there are no non-mirror improper rotations, we can directly infer the irrep of
        # the 2T branches. If that irrep is regular (i.e. has no negative coefficients), we
        # don't need to invoke 1L at all, and can solve for just 2T alone.
        return find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                               safetychecks=true, verbose=verbose)

    else 
        # Two possible scenarios (treat with same approach):
        #   - All symvals known & irregular 2T and regular 1L irreps
        #   - Not all symvals known; multiple irrep options
        ms¹ᴸ = find_representation¹ᴸ(lgirs)
        ms   = ms²ᵀ .+ ms¹ᴸ
        @assert all(ms .== ms²ᵀ .+ ms¹ᴸ)
        @assert all(≥(0), ms¹ᴸ)                        # check: 1L irrep regular (Γ₁)
        @assert ms == find_representation²ᵀ⁺¹ᴸ(lgirs)  # →: [2T+1L] = 2T+1L

        return find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                                               verbose=verbose, 
                                               shuffle_1Lpick=shuffle_1Lpick)
    end
    # TODO: The returned coefficients cⁱs do not necessarily each describe different 
    #       symmetry vectors n, since the Hilbert basis is not linearly independent.
    #       We should consider whether it would be better to only return expansions that 
    #       result in unique symmetry vectors, i.e. in general a subset of cⁱs
end

function find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                         verbose::Bool=true, safetychecks::Bool=false)
    verbose && println("SG ", sgnum)

    sb, Γidxs = compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    νsᴴ = fillings(sb)
    νᴴₘₐₓ = maximum(νsᴴ)

    # We seek an expansion with coefficients cᵢ≥0 such that
    #   P(Γ) ∑ᵢ cᵢ 𝐧ᴴᵢ ≥ 𝐦(Γ)
    # where P(Γ) projects out the Γ-irreps from the Hilbert bases 𝐧ᴴᵢ. In code, letting
    # `nsᴴ = matrix(sb)`, this means we seek a solution with `nsᴴ[Γidxs,:]*c ≥ ms`. 
    # Finally, we impose a filling
    # constraint, such that the overall number of bands is at most ν. In code, this requires
    # that `nsᴴ[end,:]*c == ν`. Moreover, all 𝐧ᴴᵢ that does not have at least one nonzero
    # entry matching `ms` will not help us in fulfilling these constraints in a nontrivial
    # way, so we can ignore those (would correspond to just stacking on some bands).
    # Finally, we can restrict the sum to contain at most two 𝐧ᴴᵢ (same or different): if we
    # have more, then at least one of them isn't actually needed to fulfil ``𝐧(Γ) ≥ 𝐦(Γ)``,
    # and can then be considered a trivial stacking.

    # the "nontrivial" parts of `nᴴ` must have at least one positive element for the same 
    # irrep as a nonzero index of `ms`; we can ignore all the others
    ntidxs²ᵀ = find_symmetry_constrained_bases(sb, ms²ᵀ, Γidxs)

    cⁱs = Vector{Int}[]
    maxterms = 2
    for ν²ᵀᵗ in 2:2νᴴₘₐₓ # target filling for 2T branches (≥2)
        cⁱs = filling_symmetry_constrained_expansions(ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs, 
                                                ntidxs²ᵀ, # include only "nontrivial" bases
                                                maxterms) # limit to two basis terms

        if !isempty(cⁱs)
            verbose      && println("   ⇒ νᵀ = ", ν²ᵀᵗ, ": ", length(cⁱs), " solutions")            
            safetychecks && safetycheck²ᵀ(cⁱs, ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
            
            return cⁱs, ν²ᵀᵗ, sb, nothing
        end
    end
    throw("Found no valid expansions consistent with constraints")
end

function find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                                         verbose::Bool=false, shuffle_1Lpick::Bool=false)
    verbose && print("SG ", sgnum)

    sb, Γidxs = compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    Nⁱʳʳ = length(first(sb))
    notΓidxs = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]
    νsᴴ = fillings(sb)
    νᴴₘᵢₙ, νᴴₘₐₓ = extrema(νsᴴ)
    
    # Here, the irrep of 1L is regular (Γ₁) and the irrep of 2T is irregular (i.e. has 
    # negative coefficients). As a result, it is impossible to expand 2T's irrep in the
    # Hilbert basis since it has strictly positive elements and coefficients. We can still
    # can try to find an expansion for 2T+1L simultaneously.
    ntidxs¹ᴸ  = find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
    _, pick¹ᴸ = findmin(νsᴴ[ntidxs¹ᴸ])
    if shuffle_1Lpick 
        if length(ntidxs¹ᴸ) > pick¹ᴸ
            pick¹ᴸ += 1
        else 
            @warn "Unable to shuffle 1L pick"
        end
    end
    idx¹ᴸ     = ntidxs¹ᴸ[pick¹ᴸ] # TODO: Test that resulting expansions for 2T are invariant wrt. to this choice
    nᴸ = sb[idx¹ᴸ]
    νᴸ = νsᴴ[idx¹ᴸ]
    
    verbose && println(" (νᴴₘᵢₙ = ", νᴴₘᵢₙ, ", νᴸ = ", νᴸ, ")")
    # find _all_ feasible solutions to ms constraints for fixed and minimal νᵗ; we use a 
    # recursive looping construct to find candidate expansions
    max_patience_νᵗ = max(4*νᴴₘₐₓ, 8)
    for νᵗ in 3:max_patience_νᵗ # target filling (≥3) (function returns from loop)
        verbose && print("   … νᵗ = ", νᵗ, ": ")

        # Find the solutions to c₁νᴴ₁ + c₂νᴴ₂ + ... = νᵗ subject to the 2T+1L ms constraint
        # The below basically uses recursion to do a nested set `max_terms` loops which
        # solves linear Diophantine equation and checks symmetry constraints as well; the
        # maximum number of included bases in a valid expansion is div(νᵗ, νᴴₘᵢₙ, RoundDown)
        cⁱs = filling_symmetry_constrained_expansions(νᵗ, ms, νsᴴ, sb, Γidxs)
        verbose && println(length(cⁱs), " candidates")
        
        # Proceed to check combinations of nᴸ and n=sum(sb[cⁱ])
        cⁱs_valid = Vector{Int}[]
        n = similar(first(sb)) # candidate solution buffer     
        for cⁱ in cⁱs # 2T+1L constraints
            sum_symbases!(n, sb, cⁱ) # compute new candidate vector from cⁱ indices
            # test 1: n(∉Γ)-nᴸ(∉Γ) ≥ 0
            if all(≥(0), @views n[notΓidxs] .- nᴸ[notΓidxs])
                # test 2: [n(Γ)-n¹ᴸ⁺²ᵀₚᵢₙ(Γ)] - [nᴸ(Γ)-n¹ᴸₚᵢₙ(Γ)] ≥ 0
                if all(≥(0), (n[Γidxs] .- ms) .- (nᴸ[Γidxs] .- ms¹ᴸ)) 
                    push!(cⁱs_valid, cⁱ) # found a valid solution; push to storage
                end
            end
        end

        if !isempty(cⁱs_valid)
            νᵀ = νᵗ - νᴸ
            verbose && println("   ⇒ νᵀ = ", νᵀ, ": ", length(cⁱs_valid), " solutions")
            return cⁱs_valid, νᵀ, sb, idx¹ᴸ
        end
    end
    
    throw("Found no valid expansions consistent with constraints")
end

# -----------------------------------------------------------------------------------------

end # module