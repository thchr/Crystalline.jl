using Crystalline
import Crystalline: DEFAULT_ATOL, rotation, irdim
using LinearAlgebra, Test, JuMP, GLPK #, Nemo

if !isdefined(Main, :SymmetryBases)
    include((@__DIR__)*"/../src/SymmetryBases/SymmetryBases.jl")
end
using Main.SymmetryBases
const PyNormaliz = Main.SymmetryBases.PyNormaliz

# ------------------------------------------------------------------------------------------

function find_bandrepresentation_lowest_bands(
            sgnum::Integer; 
            timereversal::Bool=true, 
            verbose::Bool=true
            )
    # Irreps at Γ, irrep-multiplicities of ω=0 2T bands, and symmetry operations
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    timereversal && (lgirs = realify(lgirs))
    lg = group(first(lgirs))
    rotvals = map(op->(W=rotation(op); Crystalline.rotation_order_3d(det(W), tr(W))), lg)

    # 2T irreps; check if "simple treatment"/fast-path is applicable
    ms²ᵀ = find_zero_freq_gamma_transverse_representation(lgirs)
    has_nonmirror_improper = any(∈((-1, -3, -4, -6)), rotvals)
    is_regular²ᵀ = all(≥(0), ms²ᵀ)
    if !has_nonmirror_improper && is_regular²ᵀ # All symvals known & regular 2T irrep

        # Simple case: if there are no non-mirror improper rotations, we can directly infer 
        # the irrep of the 2T branches. If that irrep is regular (i.e. has no negative
        # coefficients), we don't need to invoke 1L at all, and can solve for just 2T alone.
        msˣ = ms²ᵀ
        return find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                              safetychecks=true, verbose=verbose)

    else#if !has_nonmirror_improper # All symvals known & irregular 2T and regular 1L irreps
        ms¹ᴸ = find_zero_freq_gamma_longitudinal_representation(lgirs)
        ms = find_zero_freq_gamma_representation(lgirs)
        @assert all(ms .== ms²ᵀ .+ ms¹ᴸ)
        @assert all(≥(0), ms¹ᴸ)                              # check: 1L irrep regular (Γ₁)
        @assert ms == find_zero_freq_gamma_representation(lgirs) # →: [2T+1L] = 2T+1L

        return find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                                               verbose=verbose)

    #else # Not all symvals known; multiple irrep options
        # TODO: implement non-fast path treatment
        # Main idea is to explore every choice for 1L and 2T that is _feasible_ given the 
        # existing symmetry constraints.
    #    return nothing, nothing, nothing
    end
    #=
    sb, Γidxs = compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    nsᴴ = matrix(sb)
    maycontain = allowable_nsᴴ_idxs(msˣ, nsᴴ, Γidxs)
    Nᴴ = size(nsᴴ, 2)

    # Try to find a nonnegative expansion in minimal "band-numbers" of basis of `nsᴴ` 
    # subject to constraints from `ms`; define as a constrained feasibility problem
    m = Model(GLPK.Optimizer)
    @variable(m, c[1:Nᴴ] >= 0, Int)
    @constraint(m, nsᴴ[Γidxs,:]*c .>= ms)

    # Define `filling_constraint` as a variable, so its value can be changed on the fly
    @variable(m, filling_constraint)
    @constraint(m, nsᴴ[end,:]'*c == filling_constraint)
    for filling = 2:MAX_FILL
        # Impose constraints on the number of included bands (band filling)
        fix(filling_constraint, filling) # Change/set `filling_constraint` value
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found
            verbose && println("   Found feasible solution with ν = ", filling, " bands")
            return m, filling
        end

        # I think this might be doable instead (or optimizable) by just looking at the 
        # individual bases, and excluding those that do not have any elements at _some_ of 
        # the required Γ irreps. Might still need an optimization step though. Wonder if we
        # can try to find all the solutions by excluding a solution once we've found it?
    end
    return m, -1
    =#
end

function find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                        verbose::Bool=true, safetychecks::Bool=false)
    verbose && println("SG ", sgnum)

    sb, Γidxs = compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    nsᴴ = matrix(sb)
    νsᴴ = fillings(sb)
    νᴴₘₐₓ = maximum(νsᴴ)
    Nᴴ = length(sb)

    # We seek an expansion with coefficients cᵢ≥0 such that
    #   P(Γ) ∑ᵢ cᵢ 𝐧ᴴᵢ ≥ 𝐦(Γ)
    # where P(Γ) projects out the Γ-irreps from the Hilbert bases 𝐧ᴴᵢ. In code, this 
    # means we seek a solution with `nsᴴ[Γidxs,:]*c ≥ ms`. Finally, we impose a filling
    # constraint, such that the overall number of bands is at most ν. In code, this requires
    # that `nsᴴ[end,:]*c == ν`. Moreover, all 𝐧ᴴᵢ that does not have at least one nonzero
    # entry matching `ms` will not help us in fulfilling these constraints in a nontrivial
    # way, so we can ignore those (would correspond to just stacking on some bands).
    # Finally, we can restrict the sum to contain at most two 𝐧ᴴᵢ (same or different): if we
    # have more, then at least one of them isn't actually needed to fulfil ``𝐧(Γ) ≥ 𝐦(Γ)``,
    # and can then be considered a trivial stacking.

    # the "nontrivial" parts of `nᴴ` must have at least one positive element for the same 
    # irrep as a nonzero index of `ms`; we can ignore all the others
    ntidxs²ᵀ_nsᴴ = find_symmetry_constrained_bases(sb, ms²ᵀ, Γidxs)

    cⁱs = Vector{Int}[]
    ms′ = similar(ms²ᵀ) # buffer
    νᵗ = 2 # target filling (≥2)
    while isempty(cⁱs) && νᵗ ≤ 2νᴴₘₐₓ
        for idxᵢ in eachindex(ntidxs²ᵀ_nsᴴ)
            i = ntidxs²ᵀ_nsᴴ[idxᵢ]
            # add nᴴᵢ+nᴴⱼ to `cⁱs` if consistent w/ filling and symmetry constraints
            ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,), νᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
            ν ≥ νᵗ && continue # if ν already meets or overflows νᵗ we cannot add more

            # try to add one more Hilbert basis vector, i.e. test nᴴᵢ+nᴴⱼ
            for idxⱼ in idxᵢ:length(ntidxs²ᵀ_nsᴴ)
                j = ntidxs²ᵀ_nsᴴ[idxⱼ]
                # add nᴴᵢ+nᴴⱼ to `cⁱs` if consistent w/ filling and symmetry constraints
                test_expansion_add_if_valid!(cⁱs, ms′, (i,j), νᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
            end
        end
        νᵗ += 1 # increment target filling
    end
    isempty(cⁱs) && throw("Found no valid expansions consistent with constraints")
    νₘᵢₙ = νᵗ - 1

    if safetychecks
        # check that all solutions are valid and unique
        for cⁱ in cⁱs
            isvalid_solution(cⁱ, sb, ms²ᵀ, νₘᵢₙ, Γidxs) || throw("Found invalid solution")
        end
        allunique(cⁱs) || throw("Found repeated solutions, unexpectedly")
    
        # We can also compare our solutions with a complementary approach:
        cⁱs′ = solve_from_linear_diophantine_eqs(νsᴴ, νₘᵢₙ, sb, ms²ᵀ, Γidxs)
        #cs  = idxs2coef.(cⁱs, Nᴴ) # convert our solutions to coefficient vectors to compare
        Set(cⁱs) ≠ Set(cⁱs′) && throw("Check failed: did not obtain equivalent solutions")
    end       

    # print some stuff, if requested
    verbose && println("   νᵀ = ", νₘᵢₙ, ": ", length(cⁱs), " solutions")

    return cⁱs, νₘᵢₙ, nsᴴ
end

function find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                    verbose=verbose)
    verbose && println("SG ", sgnum)

    sb, Γidxs = compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    nsᴴ = matrix(sb)
    Nⁱʳʳ = size(nsᴴ, 1)
    notΓidxs = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]
    νsᴴ = fillings(sb)
    νᴴₘₐₓ = maximum(νsᴴ)
    νᴴₘᵢₙ = minimum(νsᴴ)   
    
    # Here, the irrep of 1L is regular (Γ₁) and the irrep of 2T is irregular (i.e. has 
    # negative coefficients). As a result, it is impossible to expand 2T's irrep in the
    # Hilbert basis since it has strictly positive elements and coefficients. We can still
    # can try to find an expansion for 1L+ 2T simultaneously.

    #ntidxs_nsᴴ   = find_symmetry_constrained_bases(sb, ms,   Γidxs)
    ntidxs_nsᴴ   = 1:size(nsᴴ, 2)
    ntidxs¹ᴸ_nsᴴ = find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
    _, pick = findmin(νsᴴ[ntidxs¹ᴸ_nsᴴ])
    ntidxs¹ᴸ_nsᴴ = ntidxs¹ᴸ_nsᴴ[pick:pick]

    νsᴸ = @view νsᴴ[ntidxs¹ᴸ_nsᴴ]
    nsᴸ = @view sb[ntidxs¹ᴸ_nsᴴ]
    νᴸ = maximum(νsᴸ)
    println("   νᴴₘᵢₙ = ", νᴴₘᵢₙ, ", νᴸ = ", νᴸ)
    # find _all_ feasible solutions to ms constraints for fixed and minimal νᵗ; can at 
    # include any number of Hilbert bases in general - we have a fast path for less than 
    # 4 bases, and otherwise fall back to PyNormaliz
    ms′ = similar(ms) # buffer
    idxsᴸ_keep = Int[]
    νsᵀ = Int[]
    νᵗ = 3 # target filling (≥3)
    fargs = (ms, νsᴴ, sb, Γidxs) # fixed args to test_expansion_add_if_valid
    while (isempty(νsᵀ) || !(minimum(νsᵀ) == 2 || νᵗ-νᴸ ≥ minimum(νsᵀ))) && 
           νᵗ ≤ 4*νᴴₘₐₓ

        verbose && print("   νᵗ = ", νᵗ, ": ")
        # determine the maximum number of basis terms for an expansion with filling νᵗ
        max_terms = div(νᵗ,νᴴₘᵢₙ, RoundDown)

        # Determine the solutions to c₁νᴴ₁ + c₂νᴴ₂ + ... = νᵗ subject to the ms constraint        
        if max_terms ≤ 4
            # Fast path by manual looping
            # TODO: Limit number of loops to `max_terms`? Maybe use generated functions...
            cⁱs = Vector{Int}[]
            for idxᵢ in eachindex(ntidxs_nsᴴ)
                i = ntidxs_nsᴴ[idxᵢ]
                ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,), νᵗ, fargs...)
                ν ≥ νᵗ && continue # if ν already meets or overflows νᵗ we cannot add more
                for idxⱼ in idxᵢ:length(ntidxs_nsᴴ)
                    j = ntidxs_nsᴴ[idxⱼ]
                    ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,j), νᵗ, fargs...)
                    ν ≥ νᵗ && continue
                    for idxₖ in idxⱼ:length(ntidxs_nsᴴ)
                        k = ntidxs_nsᴴ[idxₖ]
                        ν = test_expansion_add_if_valid!(cⁱs, ms′, (i,j,k), νᵗ, fargs...)
                        ν ≥ νᵗ && continue
                        for idxₗ in idxₖ:length(ntidxs_nsᴴ)
                            l = ntidxs_nsᴴ[idxₗ]
                            test_expansion_add_if_valid!(cⁱs, ms′, (i,j,k,l), νᵗ, fargs...)
                        end
                    end
                end
            end
            println(length(cⁱs), " candidate expansions (from looping)")
        else
            # Fallback to generic solver from Normaliz for higher number of terms
            cⁱs = solve_from_linear_diophantine_eqs(νsᴴ, νᵗ, sb, ms, Γidxs)
            println(length(cⁱs), " candidate expansions (from Normaliz)")
        end

        # Proceed to check combinations of cᴸ and c
        n = similar(first(sb)) # candidate solution buffer
        for (idxᴸ, nᴸ) in enumerate(nsᴸ) # 1L constraints         
                for cⁱ in cⁱs            # 2T+1L constraints
                    sum_symbases!(n, sb, cⁱ) # compute new candidate vector from cⁱ indices
                    # test 1: n(∉Γ)-nᴸ(∉Γ) ≥ 0
                    if all(≥(0), @views n[notΓidxs] .- nᴸ[notΓidxs])
                        # test 2: [n(Γ)-n¹ᴸ⁺²ᵀₚᵢₙ(Γ)] - [nᴸ(Γ)-n¹ᴸₚᵢₙ(Γ)] ≥ 0
                        if all(≥(0), (n[Γidxs] .- ms) .- (nᴸ[Γidxs] .- ms¹ᴸ)) 
                            nᵀ = n .- nᴸ
                            νᵀ = νᵗ - νsᴸ[idxᴸ]
                            if νᵀ∉νsᵀ
                                push!(νsᵀ, νᵀ)
                                println("      Found solution with νᵀ = ", νᵀ)
                            end
                            @goto earlystop
                        end
                    end
                end
            #end
        end
        @label earlystop
        νᵗ += 1
    end
    
    if !isempty(νsᵀ)
        min_νᵀ = minimum((νsᵀ))
    else
        throw("   Could not satisfy requirements")
    end

    return nothing, min_νᵀ, nothing
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

function idxs2coef(cⁱ, Nᴴ)
    c = zeros(Int, Nᴴ)
    for i in cⁱ
        c[i] += 1
    end
    return c
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


function test_expansion_add_if_valid!(cⁱs, ms′, # push to cⁱs; use ms′ as an updating buffer
                                      ijks::NTuple{N,Int}, νᵗ, ms, νsᴴ, sb, Γidxs) where N
    if N == 1
        i, = ijks
        ν = νsᴴ[i]                                  # update filling
        ν ≠ νᵗ && return ν                          # return early if ν overflows νᵗ
        @inbounds @views ms′ .= ms .- sb[i][Γidxs]  # update Γ-constraints
    elseif N == 2
        i,j = ijks
        ν = νsᴴ[i] + νsᴴ[j]
        ν ≠ νᵗ && return ν
        @inbounds @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs]
    elseif N == 3
        i,j,k = ijks
        ν = νsᴴ[i] + νsᴴ[j] + νsᴴ[k]
        ν ≠ νᵗ && return ν
        @inbounds @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs] .- sb[k][Γidxs]
    elseif N == 4
        i,j,k,l = ijks
        ν = νsᴴ[i] + νsᴴ[j] + νsᴴ[k] + νsᴴ[l]
        ν ≠ νᵗ && return ν
        @inbounds @views ms′ .= ms .- sb[i][Γidxs] .- sb[j][Γidxs] .- sb[k][Γidxs] .- sb[l][Γidxs]
    else
        throw("Unexpected combination of more than three Hilbert bases")
    end

    # check if nᴴᵢ+nᴴⱼ+nᴴₖ has filling νᵗ and fulfil symmetry constraints from `ms`
    if all(≤(0), ms′) # check if nᴴᵢ+nᴴⱼ fulfill `ms` constraints
        add_solution!(cⁱs, ijks) # push a solution "i+j" to storage `cⁱs`
    end

    return ν # return filling associated with `ijks` combination
end

# TODO: Complete this and test?
function find_solutions_recursive(cⁱs, ms′, idxs, νᵗ, ms, νsᴴ, nsᴴ, Γidxs, ntidxs, 
                                  startidx, stopidx, depth, maxdepth)
    if depth == maxdepth
        return nothing
    end
    for idxᵢ in startidx:stopidx
        i = ntidxs[idxᵢ]
        ν = test_expansion_add_if_valid!(cⁱs, ms′, (idxs...,i), νᵗ, ms, νsᴴ, nsᴴ, Γidxs)
        if ν ≥ νᵗ 
            continue
        else
            # try to add more Hilbert basis vectors
            find_solutions_recursive(cⁱs, ms′, (idxs...,i), νᵗ, ms, νsᴴ, nsᴴ, Γidxs, ntidxs,
                                     idxᵢ, stopidx, depth+1, maxdepth)
        end
    end
    return nothing
end

"""
    solve_from_linear_diophantine_eqs(νsᴴ, νᵗ, nsᴴ, ms, Γidxs)

Solves the same problem as `find_minimum_bandreps_regular2T(..)`, but about 20-250× slower,
but with a simpler surface implementation.

Returns the coefficient vector `cs` rather than indices into a coefficient vector. The
associated symmetry vector is thus `n = nsᴴ*c`.
"""
function solve_from_linear_diophantine_eqs(νsᴴ, νᵗ, sb, ms, Γidxs)
    cⁱs = Vector{Int}[]
    #nsᴴ_Γ = matrix(sb)[Γidxs,:]
    cⁱs_νᵗ = filling_constrained_nsᴴ_expansion(νsᴴ, νᵗ, verbose=false)
    n = similar(first(sb))
    for cⁱ in cⁱs_νᵗ
        sum_symbases!(n, sb, cⁱ) # set n = nsᴴ_Γ*c
        if all((@view n[Γidxs]) .≥ ms)
            push!(cⁱs, cⁱ)
        end
    end

    return cⁱs
end

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

function isvalid_solution(idxs::Vector{Int}, sb::SymBasis, ms::Vector{Int}, νᵗ::Int, Γidxs)
    n = sum_symbases(sb, idxs)
    return all(n[Γidxs] .≥ ms) && n[end] == νᵗ
end


"""
    filling_constrained_nsᴴ_expansion(νsᴴ::AbstractVector{<:Int}, νᵗ::Int)

Find all non-negative integer solutions ``{cᵢ}`` to the linear Diophantine equation

> ``c₁ν₁ᴴ + c₂v₂ᴴ + ... =`` `νᵗ`

with `νsᴴ` ``= [ν₁ᴴ, ν₂ᴴ, ...]`` denoting the fillings associated with a Hilbert basis.

Solutions are returned as a `::Vector{Vector{Int}}`. Uses PyNormaliz to solve the integral
polytope defined by the above inhomogeneous equation.

Optionally prints number of solutions, if the kwarg `verbose::Bool=false` is set to `true`.
"""
function filling_constrained_nsᴴ_expansion(νsᴴ::AbstractVector{Int}, νᵗ::Int;
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

function find_zero_freq_gamma_transverse_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_transverse_symvals(lg)

    return find_representation(symvals, lgirs, nothing)
end
function find_zero_freq_gamma_longitudinal_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_longitudinal_symvals(lg)

    return find_representation(symvals, lgirs, nothing)
end
function find_zero_freq_gamma_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_symvals(lg)

    return find_representation(symvals, lgirs)
end
# convenience accessors
for f in (:find_zero_freq_gamma_transverse_representation, 
          :find_zero_freq_gamma_longitudinal_representation,
          :find_zero_freq_gamma_representation)
    @eval function $f(sgnum::Integer; timereversal::Bool=false)
        lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
        timereversal && (lgirs = realify(lgirs))

        return $f(lgirs)
    end
end

# ------------------------------------------------------------------------------------------
# (sum of) symmetry eigenvalues of ω=0 branches
function zero_freq_gamma_transverse_symvals(ops::AbstractVector{SymOperation{3}})
    symvalsᵀ = Vector{ComplexF64}(undef, length(ops))

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        
        if rotval != -1     # not inversion
            n = abs(rotval) # rotation order 
            isimproper = signbit(rotval)

            if !isimproper  # ← "Proper" rotation
                # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
                # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
                θ = 2π/n
                symvalsᵀ[i] = 2cos(θ) # eⁱᶿ + e⁻ⁱᶿ = 2cos(θ)
            else            # ← Roto-inversion or mirror
                # This is not generally possible to infer for transverse plane waves alone. 
                # E.g., by direct example: there are
                # no lines of symmetry from Γ that contain -4₀₀₁ in sg 81; nor any with 
                # -3₀₀₁ in sg 147; nor any with -6₀₀₁ in sg 174. As such, there is no line 
                # of symmetry from Γ along which transverse plane waves could be 
                # symmetry-allowed eigenfunctions for the rotoinversions

                # It _is_ possible for a simple mirror though (i.e., rotation followed by
                # inversion, i.e. -2 === m): the right choice is to pick the symmetry
                # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
                # each high-symmetry k-vector)
                if rotval == -2
                    symvalsᵀ[i] = zero(ComplexF64) # (+1) + (-1)
                elseif rotval ∈ (-3, -4, -6)
                    # SGs with non-mirror rotoinversions 81:88, 111:142, 147:148, 162:167, 
                    # 174:176, 187:194, 200:206, and 215:230
                    θ = 2π/n
                    symvalsᵀ[i] = -2cos(θ) - 2.0 # [2T+1L] - 1L = [(-eⁱᶿ) + (-e⁻ⁱᶿ) + (-1)] - (+1)
                else
                    throw("Unexpected rotation value $rotval")
                end
            end
        else                # ← Inversion
            # This is a bit of a special case because, like rotoinversions, there are no
            # lines of symmetry from Γ that are invariant under inversion. However, if we 
            # recall that, overall, we are considering two transverse and a single
            # longitudinal plane wave branch (≡T²+L), which effectively transform like the
            # three Cartesian vectors at Γ, there should be no difficulty since, regardless
            # of partitioning into T² and L, each symmetry eigenvalue must be -1
            symvalsᵀ[i] = -4.0 # [2T+1L] - 1L = [(-1) + (-1) + (-1)] - (+1)
        end
    end

    return symvalsᵀ
end

function zero_freq_gamma_longitudinal_symvals(ops::AbstractVector{SymOperation{3}})
    symvalsᴸ = ones(ComplexF64, length(ops)) 
    #=
    Vector{ComplexF64}(undef, length(ops))
    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)  
        if rotval != -1     # not inversion
            n = abs(rotval) # rotation order 
            isimproper = signbit(rotval)
            if !isimproper  # ← Ordinary rotation
                symvalsᴸ[i] = one(ComplexF64)
            else            # ← Roto-inversion or mirror
                if rotval == -2
                    symvalsᴸ[i] = one(ComplexF64)
                elseif rotval ∈ (-3, -4, -6)
                    symvalsᴸ[i] = ComplexF64(NaN) # Indeterminate
                else
                    throw("Unexpected rotation value $rotval")
                end
            end
        else                # ← Inversion
            symvalsᴸ[i] = -one(ComplexF64)
        end
    end
    return symvalsᴸ
    =#
end

function zero_freq_gamma_symvals(ops::AbstractVector{SymOperation{3}})
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
# convenience accessors via space/little groups, ensuring primitive basis
for f in (:zero_freq_gamma_transverse_symvals, :zero_freq_gamma_symvals)
    @eval $f(sg::Union{LittleGroup{3}, SpaceGroup{3}}) = $f(operations(primitivize(sg)))
end

# ------------------------------------------------------------------------------------------

if false
    νs = getindex.(find_bandrepresentation_lowest_bands.(1:230, timereversal=true, verbose=false),
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