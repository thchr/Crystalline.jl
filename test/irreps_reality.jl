using Crystalline, Test, LinearAlgebra
using Crystalline: iscorep

if !isdefined(Main, :LGIRS)
    LGIRS = get_all_lgirreps(Val(3))
end

@testset "Coreps/physically real irreps" begin

# Compare evaluation of the Herring criterion with the tabulated realities in ISOTROPY
@testset "Herring criterion" begin
    #= error_count = 0 =#       # for debugging
    for (sgnum, lgirsd) in enumerate(LGIRS)
        sgops = operations(first(lgirsd["Γ"])) # this is important: we may _not_ use trivial repeated sets, i.e. spacegroup(..) would not work generally
        for lgirs in values(lgirsd)
            for lgir in lgirs
                iso_reality = reality(lgir)
                #= try =#       # for debugging
                @test iso_reality == calc_reality(lgir, sgops)
                #= catch err    # for debugging
                    if true
                        println(sgnum, ", ", 
                                Crystalline.subscriptify.(label(lgir)), ", ", 
                                centering(sgnum,3), ", ",
                                issymmorph(sgnum))
                    end
                    error_count += 1
                    println(err)
                end =#
            end
        end
    end

    #=                          # for debugging
    if error_count>0
        println(Crayon(foreground=:red, bold=true), "Outright errors: ", error_count)
    else
        println(Crayon(foreground=:green, bold=true), "No errors")
    end
    =#
end # @testset "Herring criterion" 

@testset "Corep orthogonality" begin
αβγ = Crystalline.TEST_αβγ # for evaluating characters/irreps at non-special points
# compute all coreps (aka "physically real" irreps, i.e. incorporate time-reversal sym) via
# the ordinary irreps and `realify`
LGIRS′ = [Dict(klab=>realify(lgirs) for (klab,lgirs) in lgirsd) for lgirsd in LGIRS]

function corep_orthogonality_factor(lgir::AbstractIrrep)
    if iscorep(lgir)
        r = reality(lgir)
        r == PSEUDOREAL && return 4
        r == COMPLEX    && return 2
        error(DomainError(r, "unreachable; invalid combination of iscorep=true and reality type"))
    else
        return 1
    end
end

## 2ⁿᵈ orthogonality theorem (characters) [automatically includes 1ˢᵗ orthog. theorem also]: 
#       ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦfNₒₚ⁽ᵃ⁾  
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same little group (with i running over the 
# Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements).
# Because we are dealing with coreps here, the orthogonality relations are modified
# somewhat. Assuming that we still only sum over the unitary operations (i.e. no gray 
# operations, i.e. products with time-inversion), we should include a multiplicative factor
# `f` in the orthogonality relations, as in Bradley & Cracknell Eq. 7.4.10. `f` equals:
#       real => 1; pseudoreal => 4; complex => 2
# see `corep_orthogonality_factor(lgir)` above.
@testset "1ˢᵗ & 2ⁿᵈ orthogonality theorem" begin
    for lgirsd in LGIRS             # lgirsd: dict of little group irrep collections
        for lgirs in values(lgirsd) # lgirs:  vector of distinct little group irreps
            Nₒₚ = order(first(lgirs))
            # compute coreps (aka "physically real" irreps, i.e. incorporate time-reversal)
            lgcoreps = realify(lgirs)
            for (a, lgir⁽ᵃ⁾) in enumerate(lgcoreps)
                χ⁽ᵃ⁾ = characters(lgir⁽ᵃ⁾, αβγ)
                f = corep_orthogonality_factor(lgir⁽ᵃ⁾)
                for (β, lgir⁽ᵝ⁾) in enumerate(lgcoreps)
                    χ⁽ᵝ⁾ = characters(lgir⁽ᵝ⁾, αβγ)
                    orthog2nd = dot(χ⁽ᵃ⁾, χ⁽ᵝ⁾) # ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾
                    @test (orthog2nd ≈ (a==β)*f*Nₒₚ)  atol=1e-12
                end
            end
        end
    end
end # @testset "2ⁿᵈ orthogonality theorem"

# TODO: Great orthogonality theorem of coreps:
#       It's not immediately obvious how to do this properly. E.g., if we just try:
#
#               ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ = δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
#
#       with the sum only running over the unitary elements, then we run into trouble when
#       (nm) and (jk) are referencing different blocks in the block diagonal form of the 
#       "realified" matrices. Of course, when they reference the same blocks, everything is
#       fine, since it's then just as if we had never "realified" them in the first place.
#       Introducing a multiplicative factor `f` as in the character orthogonality theorem
#       also isn't helpful here.
#       The issue probably is that we really do need to use the entire gray group here, i.e.
#       include the antiunitary elements in the i-sum. That's more bothersome, and we don't
#       have a use for the irreps of the antiunitary elements otherwise, so it doesn't seem
#       worthwhile...

end # @testset "Corep orthogonality"

end # @testset "Coreps/physically real irreps"
