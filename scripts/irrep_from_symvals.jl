using Crystalline
import Crystalline: DEFAULT_ATOL, rotation, irdim
using LinearAlgebra #JuMP, GLPK , Nemo

if !isdefined(Main, :SymmetryBases)
    include((@__DIR__)*"/../src/SymmetryBases/SymmetryBases.jl")
end
using Main.SymmetryBases
const PyNormaliz = Main.SymmetryBases.PyNormaliz

# ------------------------------------------------------------------------------------------

function minimal_expansion_of_zero_freq_bands(sgnum::Integer; 
                                              timereversal::Bool=true, verbose::Bool=true)

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
                                               verbose=verbose)
    end
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
                    verbose=verbose)
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
    idx¹ᴸ   = ntidxs¹ᴸ[pick¹ᴸ] # TODO: Test that resulting expansions for 2T are invariant wrt. to this choice
    nᴸ = sb[idx¹ᴸ]
    νᴸ = νsᴴ[idx¹ᴸ]
    
    verbose && println(" (νᴴₘᵢₙ = ", νᴴₘᵢₙ, ", νᴸ = ", νᴸ, ")")
    # find _all_ feasible solutions to ms constraints for fixed and minimal νᵗ; can at 
    # include any number of Hilbert bases in general - we have a fast path for less than 
    # 4 bases, and otherwise fall back to PyNormaliz
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
            return cⁱs_valid, νᵀ, idx¹ᴸ, sb
        end
    end
    
    throw("Found no valid expansions consistent with constraints")
end

# -----------------------------------------------------------------------------------------
"""
    find_symmetry_constrained_bases(sb::SymBasis, ms::AbstractVector{Int},
                                    Γidxs::AbstractVector{Int})

From a Hilbert basis, specified via `sb::SymBasis` whose elements are the Hilbert basis
vectors, find those that has at least one positive element in overlap with a set irrep 
multiplicities `ms`, whose indices in the rows of the Hilbert basis vectors are specified by
`Γidxs`.
Returns an array of indices into the the vectors of `sb`.
"""
function find_symmetry_constrained_bases(sb::SymBasis, ms::AbstractVector{Int},
                                         Γidxs::AbstractVector{Int})
    ntidxsᴴ = Int[]
    for (idx, nᴴ) in enumerate(sb)
        if has_mutual_positive_elements((@view nᴴ[Γidxs]), ms)
            push!(ntidxsᴴ, idx)
        end
    end
    return ntidxsᴴ
end

# ≡ any(x>0 & y>0) w/o allocations
has_mutual_positive_elements(x, y) = any(xy -> (xy[1] > 0) & (xy[2] > 0), zip(x,y))

function add_solution!(cⁱs::Vector{Vector{Int}}, ijks::NTuple{N, Int}) where N
    # push `ijks` to solution storage `cⁱs` as a vector of indices
    push!(cⁱs, [idx for idx in ijks])
end

function coef2idxs(c::AbstractVector{Int})
    N = sum(c)
    cⁱ = Vector{Int}(undef, N)
    pos₁, pos₂, idx = 0, 0, 0
    while true
        idx  = findnext(≠(0), c, idx+1)
        pos₁ = pos₂+1
        pos₂ = pos₂+c[idx]
        cⁱ[pos₁:pos₂] .= idx
        pos₂ == N && break
    end
    return cⁱ
end

function idxs2coef(cⁱ, Nᴴ) # Nᴴ = length(sb)
    c = zeros(Int, Nᴴ)
    for i in cⁱ
        c[i] += 1
    end
    return c
end

"""
    filling_symmetry_constrained_expansions(νᵗ::Integer, ms::AbstractVector{<:Integer}, νsᴴ,
                                            sb::SymBasis, Γidxs, 
                                            ntidxs=eachindex(sb),
                                            maxdepth=div(νᵗ, minimum(νsᴴ), RoundDown))

Given a compatibility basis `sb` with Hilbert bases ``[𝐧₁ᴴ, 𝐧₂ᴴ, ...]`` with associated
fillings `vsᴴ` ``= [ν₁ᴴ, ν₂ᴴ, ...]``, find all expansions `cⁱs` that (a) satisfy the filling
constraint (a linear Diophantine equation)

``c₁ν₁ᴴ + c₂v₂ᴴ + ... =`` `νᵗ`

with non-negative, integer coefficients {cᵢ∈ℕ} and (b) satisfy the symmetry constraint

``(𝐧 = c₁𝐧₁ᴴ + c₂𝐧₂ᴴ + ...)(Γ) ≥`` `ms`

evaluated only at the Γ-point, whose indices into the ``𝐧ᵢᴴ`` vector are specified by
`Γidxs`.

Optionally, if the caller wants to restrict the expansion to a subset of the bases in `sb`,
the argument `ntidxs` can provide an indexing into allowable bases of `sb`.

# Implementation
Recursion is used to build a nested set of for loops, of depth `maxdepth`, corresponding 
to the inclusion of at most `maxdepth` Hilbert bases (this limits the maximum meaningful 
value of `maxdepth` to `div(νᵗ, minimum(νsᴴ), RoundDown)`; its default value). 

# Note
See also the `*_loop` and `*_normaliz` legacy methds that achieve the same goal by different
means. They are retained in the codebase, despite being less capable or much slower, 
respectively, in the belief that they might more clearly illustrate the approach.
"""
function filling_symmetry_constrained_expansions(νᵗ::Integer, ms::AbstractVector{<:Integer},
                                        νsᴴ, sb::SymBasis, Γidxs,
                                        ntidxs=eachindex(sb),
                                        maxdepth::Integer=div(νᵗ, minimum(νsᴴ), RoundDown))

    νᵗ > 0 || throw(DomainError(νᵗ, "must be positive"))

    cⁱs = Vector{Int}[] # solution vector storage
    ms′ = similar(ms)   # buffer
    _filling_symmetry_constrained_expansions!(cⁱs, ms′, (), νᵗ, ms, νsᴴ, sb, Γidxs, 
                                              1, length(ntidxs), 1, maxdepth, ntidxs)
end
function _filling_symmetry_constrained_expansions!(cⁱs, ms′, ijks, νᵗ, ms, νsᴴ, 
                sb::SymBasis, Γidxs, startidx, stopidx, depth, maxdepth, ntidxs)
    depth > maxdepth && return cⁱs
    for idxᵢ in startidx:stopidx
        i = ntidxs[idxᵢ]
        ν = test_expansion_add_if_valid!(cⁱs, ms′, (ijks...,i), νᵗ, ms, νsᴴ, sb, Γidxs)
        ν ≥ νᵗ && continue # matched/overflowed νᵗ constraint; nothing more to add

        # did not yet match/overflow filling constraint: add more Hilbert basis vectors
        _filling_symmetry_constrained_expansions!(cⁱs, ms′, (ijks...,i), νᵗ, ms, 
                νsᴴ, sb, Γidxs, idxᵢ, stopidx, depth+1, maxdepth, ntidxs)
    end
    return cⁱs
end

function test_expansion_add_if_valid!(cⁱs, ms′, # push to cⁱs; use ms′ as an updating buffer
                                      ijks::NTuple{N,Int}, νᵗ, ms, νsᴴ, sb, Γidxs) where N

    ν = _sum_fillings(ijks, νsᴴ)                   # accumulate band fillings
    ν ≠ νᵗ && return ν                             # return early if ν overflows νᵗ
    _update_symmetry_constraints!(ms′, ijks, ms, sb, Γidxs) # update Γ-constraints in ms′

    # check if nᴴᵢ+nᴴⱼ+nᴴₖ+... fulfil symmetry constraints from `ms`
    if all(≤(0), ms′) # check if nᴴᵢ+nᴴⱼ+nᴴₖ+... fulfill `ms` constraints
        add_solution!(cⁱs, ijks) # push a solution "i+j+k+..." to storage `cⁱs`
    end

    return ν # return filling associated with `ijks` expansion
end

# equivalent of ν = νsᴴ[i] + νsᴴ[j] + νsᴴ[k] + ... for i,j,k, in ijks, recursively
_sum_fillings(ijks::NTuple{1,Int}, νsᴴ) = νsᴴ[first(ijks)]
function _sum_fillings(ijks::NTuple{N,Int}, νsᴴ) where N
    νsᴴ[first(ijks)] + _sum_fillings(Base.tail(ijks), νsᴴ)
end

# update Γ-constraints, assigning to ms′
@inline function _update_symmetry_constraints!(ms′, ijks::NTuple{N,Int}, ms, sb::SymBasis, Γidxs) where N
    if N == 1
        i, = ijks
        @views ms′ .= ms .- sb[i][Γidxs]
    elseif N == 2
        i,j = ijks
        @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs]
    elseif N == 3
        i,j,k = ijks
        @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs] .- sb[k][Γidxs]
    elseif N == 4
        i,j,k,l = ijks
        @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs] .- sb[k][Γidxs] .- sb[l][Γidxs]
    else # fall back to looping
        ms′ .= ms
        for ijk in ijks 
            @views ms′ .-= sb[ijk][Γidxs]
        end
    end
    return ms′
end

# we bother to optimize this, as it can be a bottleneck; much faster than a naive 
# implementation like `sum(sb[idxs])`
function sum_symbases!(n, sb::SymBasis, idxs)
    Nⁱʳʳ = length(n)
    n .= sb[first(idxs)] # peel off 1st iter & ensure invariance to n's inititialization
    @inbounds for idx in @view idxs[2:end]
        nᴴ = sb[idx]
        for i in 1:Nⁱʳʳ
            n[i] += nᴴ[i]
        end
    end
    return n
end
sum_symbases(sb::SymBasis, idxs) = sum_symbases!(similar(first(sb)), sb, idxs)

# -----------------------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------------------

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

function safetycheck²ᵀ(cⁱs, ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
    # check that all solutions are valid and unique
    all(cⁱ->isvalid_solution(cⁱ, ν²ᵀᵗ, ms²ᵀ, sb, Γidxs), cⁱs) || throw("Found invalid solutions")
    allunique(cⁱs) || throw("Found repeated solutions, unexpectedly")

    # Check that it didn't matter whether we excluded "trivial" basis elements or not
    cⁱs′ = filling_symmetry_constrained_expansions(ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
    Set(cⁱs) ≠ Set(cⁱs′) && throw("Did not obtain equivalent solution sets")
end             

function isvalid_solution(cⁱ::Vector{Int}, νᵗ::Int, ms::Vector{Int}, sb::SymBasis, Γidxs)
    n = sum_symbases(sb, cⁱ)
    return all(n[Γidxs] .≥ ms) && n[end] == νᵗ
end

"""
    filling_constrained_expansions(νsᴴ::AbstractVector{<:Int}, νᵗ::Int)

Find all non-negative integer solutions ``{cᵢ}`` to the linear Diophantine equation

``c₁ν₁ᴴ + c₂v₂ᴴ + ... =`` `νᵗ`

with `νsᴴ` ``= [ν₁ᴴ, ν₂ᴴ, ...]`` denoting the fillings associated with a Hilbert basis.

Solutions are returned as a `::Vector{Vector{Int}}`. Uses PyNormaliz to solve the integral
polytope defined by the above inhomogeneous equation.

Optionally prints number of solutions, if the kwarg `verbose::Bool=false` is set to `true`.
"""
function filling_constrained_expansions(νsᴴ::AbstractVector{Int}, νᵗ::Int; 
                                        verbose::Bool=false)

    νᵗ > 0 || throw(DomainError(νᵗ, "must be positive"))
    
    # We want to avoid including terms where νᵢᴴ > vᵗ since they cannot feature in a valid
    # solution anyway and actually end up slowing down the calculation significantly
    nt_idxs = findall(≤(νᵗ), νsᴴ)
    # Specify linear Diophantine equation via PyNormaliz's Cone constructor
    inhom_eqs = reshape([@view νsᴴ[nt_idxs]; -νᵗ], 1, length(nt_idxs)+1)
    #inhom_eqs = reshape([νsᴴ; -νᵗ], 1, length(νsᴴ)+1)
    P = PyNormaliz.Cone(inhom_equations = inhom_eqs)
    # Find non-negative integer solutions to the above integral polytope
    normaliz_sols = P.LatticePoints() # distinct solutions across rows

    # last column of `normaliz_sols` is a multiplier on ``-νᵗ``: should be 1, otherwise it 
    # corresponds to finding a solution that has a filling equal to a multiple of νᵗ. We 
    # filter out these solutions below.
    #cⁱs = [coef2idxs(c′[1:end-1]) for c′ in eachrow(normaliz_sols) if isone(c′[end])]
    cⁱs = [nt_idxs[coef2idxs(c′[1:end-1])] for c′ in eachrow(normaliz_sols) if isone(c′[end])]
    
    if verbose 
        println("   νᵗ = ", νᵗ, ": ", length(cⁱs), " νᵗ-constrained candidate solutions = ")
        if length(cⁱs) ≠ size(normaliz_sols, 1) 
            println("      DISCARDED \"MULTIPLES\"-SOLUTIONS W/ MULTIPLICITY = ",
                    filter(≠(1), unique(normaliz_sols[:,end])))
        end
    end

    return cⁱs
end

# ------------------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------------------
# (sum of) symmetry eigenvalues of ω=0 branches

# two transverse and two longitudinal plane waves (2T+1L)
function get_symvals²ᵀ⁺¹ᴸ(ops::AbstractVector{SymOperation{3}})
    symvals = Vector{ComplexF64}(undef, length(ops))
    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        n = abs(rotval)
        # This covers every case, including rotations, mirrors, rotoinversions, & inversion
        θ = 2π/n
        symvals[i] = sign(rotval)* (cis(θ) + cis(-θ) + one(ComplexF64))
    end
    
    return symvals
end

# single longitudinal plane wave (1L)
function get_symvals¹ᴸ(ops::AbstractVector{SymOperation{3}})
    symvals = ones(ComplexF64, length(ops)) 
end

# two transverse plane waves (2T)
function get_symvals²ᵀ(ops::AbstractVector{SymOperation{3}})
    symvals = Vector{ComplexF64}(undef, length(ops))

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        
        n = abs(rotval) # rotation order 

        if !signbit(rotval)                                     # ← Proper rotation
            # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
            # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
            θ = 2π/n
            symvals[i] = 2cos(θ) # eⁱᶿ + e⁻ⁱᶿ = 2cos(θ)

        else                                                    # ← Improper rotation
            # It is not generally possible to infer the all the symmetry eigenvalues of 
            # roto-inversions with rotval = (-1, -3, -4, -6) for the two transverse 
            # plane waves (2T) in isolation. This is because there are there no lines of
            # symmetry from Γ along which 2T could be symmetry-allowed eigenfunctions 
            # for the rotoinversions.
            # Instead, we pick a _possible_ choice for 1L and infer _possible_ symmetry
            # values from [2T+1L] - 1L

            # It _is_ possible for a simple mirror though (i.e., rotation followed by
            # inversion, i.e. -2 === m): the right choice is to pick the symmetry
            # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
            # each high-symmetry k-vector)
            if rotval == -2                     # ← Mirror
                symvals[i] = zero(ComplexF64) # (+1) + (-1)

            elseif rotval ∈ (-1, -3, -4, -6)    # ← Roto-inversions & inversion
                θ = 2π/n
                # In general, we have: 
                #   [2T+1L] - 1L = [(-eⁱᶿ) + (-e⁻ⁱᶿ) + (-1)] - (+1) = -2cos(θ) - 2.0
                # For inversion specifically, this is:
                #   [2T+1L] - 1L = [(-1) + (-1) + (-1)] - (+1) = -4
                symvals[i] = -2cos(θ) - 2.0 
                # SGs w/ non-mirror rotoinversions are 81:88, 111:142, 147:148, 162:167, 
                # 174:176, 187:194, 200:206, and 215:230
            end
        end
    end

    return symvals
end

# convenience accessors via space/little groups, ensuring primitive basis
for f in (:get_symvals²ᵀ⁺¹ᴸ, :get_symvals¹ᴸ, :get_symvals²ᵀ)
    @eval $f(sg::Union{LittleGroup{3}, SpaceGroup{3}}) = $f(operations(primitivize(sg)))
end

# ------------------------------------------------------------------------------------------

if false
    νs = getindex.(minimal_expansion_of_zero_freq_bands.(1:230, timereversal=true, verbose=false),
                   2)

    # Compare with Watanabe & Lu
    Base.ndigits(::Nothing) = 1 # haaaack
    include("scripts/watanabelu_results.jl") # loads Watanabe & Lu data (in `Msᵂᴸ`)
    Q = [[sg, M, Mbound] for (sg, M, Mbound) ∈ zip(1:230, νs, getindex.(Msᵂᴸ, 2))]
    Q′ = filter(x-> x[2]!==nothing, Q) # filter out those sgs that are not currently implemented (i.e. allow only regular 2T)
    issues = map(x->x[2]===nothing ? "─" : (x[2]≥(x[3]) ? " " : "!"), Q)
    differences = map(x->x[2]===nothing ? "─" : (x[2]==(x[3]) ? " " : "≠"), Q)

    foreach(vcat.(Q, issues, differences)) do x
        println("|", " "^(4-ndigits(x[1])), x[1], " |", " "^(3-ndigits(x[2])),  # SG no.
                x[2] === nothing ? "─" : x[2], " | ",  # our M predictions
                x[3] == 2 ? "=" : "≥", x[3], " | ",    # M-bounds from Watanabe & Lu
                x[4], " | ",                           # bound violations
                x[5], " |"                             # differences from bound?
        )
    end
end