# Legacy implementations of `filling_symmetry_constrained_expansions`, now superceded by a 
# recursive implementation 

const PyNormaliz = SymmetryBases.PyNormaliz # TODO: Drop when/if Normaliz.jl matures

"""
    filling_symmetry_constrained_expansions_loop(νᵗ, ms, νsᴴ, sb, Γidxs)

Legacy method: see `filling_symmetry_constrained_expansions`. 

Limited to expansions of at most 4 Hilbert bases elements.
"""
function filling_symmetry_constrained_expansions_loop(νᵗ, ms, νsᴴ, sb, Γidxs)
    args = (νᵗ, ms, νsᴴ, sb, Γidxs)
    cⁱs = Vector{Int}[] # solution vector storage
    ms′ = similar(ms)   # buffer
    for i in eachindex(sb)
        ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,), args...)
        ν ≥ νᵗ && continue # if ν already meets or overflows νᵗ we cannot add more
        for j in i:length(sb)
            ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,j), args...)
            ν ≥ νᵗ && continue
            for k in j:length(sb)
                ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,j,k), args...)
                ν ≥ νᵗ && continue
                for l in k:length(sb)
                    test_expansion_add_if_valid!(cⁱs, ms′, (i,j,k,l), args...)
                end
            end
        end
    end
    return cⁱs
end

"""
    filling_symmetry_constrained_expansions_normaliz(νᵗ, ms, νsᴴ, sb, Γidxs)

Legacy method: see `filling_symmetry_constrained_expansions`. 
Uses Normaliz to solve the linear Diophantine equation that defines the filling constraint.

Is at least about 20-250× slower than the recursive and looping implementations; sometimes
much slower. The surface implementation, however, is considerably simpler.
"""
function filling_symmetry_constrained_expansions_normaliz(νᵗ, ms, νsᴴ, sb, Γidxs)
    cⁱs = Vector{Int}[]
    # all possible solutions to the filling constraint
    cⁱs_νᵗ = filling_constrained_expansions(νsᴴ, νᵗ, verbose=false)
    n = similar(first(sb))
    for cⁱ in cⁱs_νᵗ
        sum_symbases!(n, sb, cⁱ) # set n = nsᴴ_Γ*c
        if all((@view n[Γidxs]) .≥ ms)
            push!(cⁱs, cⁱ)
        end
    end

    return cⁱs
end