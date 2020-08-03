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
See also the `*_loop` and `*_normaliz` methods in src/legacy_constrained_expansions.jl
that achieve the same goal by different means. They are retained, unloaded, in the codebase,
despite being less capable or much slower, respectively, in the belief that they might more
provide a simpler illustration of the conceptual approach.
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