# ---------------------------------------------------------------------------------------- #
#                  Methods to test the topology of a given symmetry vector n
#                        and to analyze the properties of sets of EBRs                     #
# ---------------------------------------------------------------------------------------- #

@doc"""
$(TYPEDSIGNATURES)

Enum type enumerating the possible "coarse" topological classifications diagnosable by
symmetry within the topological quantum chemistry / symmtry indicator frameworks.
"""
@enum TopologyKind begin
    TRIVIAL    = 0
    NONTRIVIAL = 1
    FRAGILE    = 2
end

# -----------------------------------------------------------------------------------------
# Trivial/nontrivial solution topology via Smith/BandRepSet/Collection{<:NewBandRep}

@doc """
$(TYPEDSIGNATURES)

Return whether a symmetry vector `n` is a band-combination that is trivial or nontrivial
from a symmetry perspective, i.e. whether it has an integer-coefficient expansion in the
elementary band representation (EBR) basis or not (i.e. a rational-coefficient expansion).
Returns a value from the Enum [`TopologyKind`](@ref) (`TRIVIAL` or `NONTRIVIAL`).

No distinction is made between fragile and trivial symmetry vectors: i.e., a `TRIVIAL`
return value may in fact be a `FRAGILE` state on more careful inspection: such a distinction
can be made by `calc_detailed_topology` from
[SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl)

See also [`symmetry_indicators`](@ref) to obtain the associated symmetry indicators of a
nontrivial symmetry vector.

## Input

The EBR basis can be provided as `::BandRepSet`, `::Matrix{<:Integer}`, or a `Smith`
decomposition.
The length of `n` must equal the EBR basis' number of irreps or the number of irreps plus 1
(i.e. include the band connectivity).

## Implementation

We check whether an integer-coefficient expansion exists via the Smith normal decomposition
of the EBR matrix ``\\mathbf{B} = \\mathbf{S}\\boldsymbol{\\Lambda}\\mathbf{T}``. If

```math
    (\\mathbf{S}^{-1}\\mathbf{n})_j = 0 \\mod \\lambda_j
```

for all ``j = 1, \\ldots, d^{\\text{bs}}`` (``d^{\\text{bs}}`` is the number of non-zero
diagonal elements of ``\\boldsymbol{\\Lambda}``, i.e. the invariant factors of
``\\mathbf{B}``), there exists a solution to ``\\mathbf{B}\\mathbf{c} = \\mathbf{n}`` with
integer coefficients ``c_j \\in \\mathbb{Z}``.

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `iscompatible(n, [...]) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_negative` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function calc_topology(
    n::AbstractVector{<:Integer},
    F::Smith;
    allow_incompatible::Bool=false,
    allow_negative::Bool=false
)

    if !allow_incompatible && !iscompatible(n, F; allow_negative)
        _throw_incompatible_or_negative(n)
    end
    
    S⁻¹, Λ = F.Sinv, F.SNF # Λ = [λ₁, …, λ_{dᵇˢ}, 0, …, 0]
    dᵇˢ = count(!iszero, Λ)

    # n is trivial if (S⁻¹n)ⱼ = 0 mod λⱼ for j = 1, …, dᵇˢ. This is equivalent to checking
    # whether there exists an integer coefficient expansion for `n` in the EBR basis that
    # `F` represents (i.e., whether `Nemo.cansolve(B, n) == true`) but faster.
    # We do the matrix-vector product row-wise to check `mod(S⁻¹[1:dᵇˢ]*n)[i], Λ[i]) = 0`
    # for `i ∈ 1:dᵇˢ` without allocating unnecessarily
    is_trivial = all(1:dᵇˢ) do i
        Λᵢ = Λ[i]
        Λᵢ == 1 && return true # fast path: `mod(x, 1) = 0` for all integer `x`.
        S⁻¹ᵢ = @view S⁻¹[i,:]
        mod(dot(S⁻¹ᵢ, n), Λᵢ) == 0
    end
    return is_trivial ? TRIVIAL : NONTRIVIAL
end

function calc_topology(n::AbstractVector{<:Integer}, B::AbstractMatrix{<:Integer}; kws...)
    length(n)≠size(B, 1) && throw(DimensionMismatch("incompatible dimensions of `n` & `B`"))
    return calc_topology(n, smith(B); kws...)
end

function calc_topology(
    n::AbstractVector{<:Integer},
    brs::Union{Collection{<:NewBandRep}, BandRepSet};
    kws...
)
    B = stack(brs)
    if !includes_connectivity(n, brs)
        return calc_topology(n, (@view B[1:end-1, :]); kws...)
    end
    return calc_topology(n, B; kws...)
end

# -----------------------------------------------------------------------------------------
# Stable topological indices

@doc """
$(TYPEDSIGNATURES)

Return the symmetry indicator indices of a symmetry vector `n`, in the context of a set of
elementary band representations (EBRs) `brs`, provided as a `Collection{<:NewBandRep}`, a
`BandRepSet`, a `Matrix{<:Integer}`, or a `Smith` decomposition thereof.

In detail, the method returns the nontrivial indices ``[\\nu_1, \\ldots, \\nu_n]``
associated with the symmetry indicator group (see, [`indicator_group`](@ref))
``[\\lambda_1, \\ldots, \\lambda_n]`` of the EBR basis.
The indices ``\\nu_i`` are elements of a cyclic group of order ``\\lambda_i``, i.e. 
``\\nu_i ∈ \\mathbb{Z}_{\\lambda_i} = \\{0, 1, \\ldots, \\lambda_i-1\\}``.

See also `calc_topology` to determine whether any symmetry indicator is nonzero (i.e.,
whether the symmetry vector is topologically nontrivial).

## Implementation

The indices are computed using the Smith normal decomposition ``\\mathbf{B} = \\mathbf{S}
\\boldsymbol{\\Lambda}\\mathbf{T}`` of the EBR matrix ``\\mathbf{B}``. 
Specifically, denoting by ``\\mathbf{s}_i^{-1}`` the ``i``th nontrivial row of
``\\mathbf{S}^{-1}``, the symmetry indicator topological indices of a symmetry vector
``\\mathbf{n}`` are computed as ``\\nu_i = \\mathbf{s}_i^{-1}\\mathbf{n}``.[^HCP]

[^HCP]: [H.C. Po, J. Phys. Cond. Matter **32**, 263001 (2020)](https://doi.org/10.1088/1361-648X/ab7adb).

## Keyword arguments
If `n` is not a compatible band structure (i.e., if `iscompatible(n, brs) = false`), an
error is thrown. This behavior can be controlled by two boolean keyword arguments:

- `allow_incompatible` (`false`): if `true`, disables the compatibility check entirely.
- `allow_negative` (`false`): if `true`, allows negative symmetry content, but maintain
  requirement that `n` respects the compatibilty relations in an algebraic sense.
"""
function symmetry_indicators(
    n::AbstractVector{<:Integer},
    F::Smith;
    allow_incompatible::Bool=false,
    allow_negative::Bool=false
)
    if !allow_incompatible && !iscompatible(n, F; allow_negative)
        _throw_incompatible_or_negative(n)
    end

    idxs = findall(x -> x≠0 && x≠1, F.SNF) # find nontrivial factor groups
    Λ    = @view F.SNF[idxs]               # nontrivial invariant factors
    S̃⁻¹  = @view F.Sinv[idxs, :]           # nontrivial rows of S⁻¹

    return mod.(S̃⁻¹*n, Λ)
end
function symmetry_indicators(
    n::AbstractVector{<:Integer},
    B::AbstractMatrix{<:Integer};
    kws...
)
    length(n)≠size(B, 1) && throw(DimensionMismatch("incompatible dimensions of `n` & `B`"))
    return symmetry_indicators(n, smith(B); kws...)
end
function symmetry_indicators(
    n::AbstractVector{<:Integer},
    brs::Union{Collection{<:NewBandRep}, BandRepSet};
    kws...
)
    B = stack(brs)
    if !includes_connectivity(n, brs)
        return symmetry_indicators(n, (@view B[1:end-1, :]); kws...)
    end
    return symmetry_indicators(n, B; kws...)
end

# ---------------------------------------------------------------------------------------- #

@doc """
$(TYPEDSIGNATURES)

Return the symmetry indicator group ``X^{\\text{BS}}`` associated with an input set of band
representations `brs` (or Smith decomposition thereof, `F`), i.e., return the the nontrivial
(i.e., ≠ {0,1}) elementary factors of the Smith normal form of the band representation
matrix. The return value is a `Vector{Int}` containing the nontrivial factors. If no
nontrivial factors exists, the return value is an empty `Vector{Int}`.

See also [`indicator_group_as_string`](@ref) for a formatted string version.

## Understanding

The symmetry indicator group answers the question "what direct product of ``\\mathbb{Z}_n``
groups is the the quotient group ``X^{\\text{BS}} = \\{\\text{BS}\\}/\\{\\text{AI}\\}``
isomorphic to?" (see e.g.,
[Po, Watanabe, & Vishwanath, Nature Commun. **8**, 50 (2017)](https://doi.org/10.1038/s41467-017-00133-2)
for more information).

## Example
```jldoctest
julia> brs = calc_bandreps(2, Val(3));

julia> indicator_group(brs)
4-element Vector{Int64}:
 2
 2
 2
 4
```
"""
function indicator_group(F::Smith)
    Λ = F.SNF
    nontriv_idx = findall(is_not_one_or_zero, Λ)
    return Λ[nontriv_idx]
end
function indicator_group(B::AbstractMatrix{<:Integer})
    F = smith(B, inverse=false)
    return indicator_group(F)
end
function indicator_group(brs::Union{Collection{<:NewBandRep}, BandRepSet})
    return indicator_group(stack(brs))
end
is_not_one_or_zero(x) = !(isone(x) || iszero(x))

@doc """
$(TYPEDSIGNATURES)

Return the symmetry indicator group ``X^{\\text{BS}}`` as a formatted string (i.e., 
as `"Zᵢ×Zⱼ×…"`). See also [`indicator_group`](@ref) for a vector representation.

## Example
```jldoctest
julia> brs = calc_bandreps(2, Val(3));

julia> indicator_group_as_string(brs)
"Z₂×Z₂×Z₂×Z₄"
```
"""
function indicator_group_as_string(nontriv_Λ::AbstractVector{<:Integer})
    if isempty(nontriv_Λ)
        return "Z₁"
    else
        io = IOBuffer()
        for (i, Λᵢ) in enumerate(nontriv_Λ)
            print(io, "Z", subscriptify(string(Λᵢ)))
            i == length(nontriv_Λ) || print(io, "×")
        end
    end
    return String(take!(io))
end
function indicator_group_as_string(
    brs::Union{Collection{<:NewBandRep}, BandRepSet, AbstractMatrix{<:Integer}, Smith}
)
    return indicator_group_as_string(indicator_group(brs))
end

# ---------------------------------------------------------------------------------------- #

@doc """
$(TYPEDSIGNATURES)

Test whether a symmetry vector `n` is a valid band grouping, i.e. whether it fulfils all
compatibility relations in the Brillouin zone and is non-negative. That is, test whether
`n` belong to the set of physical band structures {BS}.

The test compares the symmetry vector `n` to an set of elementary band representations,
provided either as a `BandRepSet`, a `Collection{<:NewBandRep}`, a `Matrix{<:Integer}`,
or a `Smith` decomposition. The irrep sorting of `n` and this set of EBRs must be identical.

## Keyword arguments

- `allow_negative` (`false`): if `true`, allows negative symmetry content. This can be
  relevant if `n` contains negative content that may nevertheless respect the compatibility
  relations in a strictly algebraic sense.

## Implementation

Belonging to {BS} is tested by comparing to a set of elementary band representations
(EBRs). A symmetry vector ``\\mathbf{n}`` is in {BS} if

```math
    \\tilde{\\mathbf{S}}\\tilde{\\mathbf{S}}^{-1}\\mathbf{n} = \\mathbf{n}
```

where ``\\tilde{\\mathbf{S}}`` (``\\tilde{\\mathbf{S}}^{-1}``) denotes the nonsingular
columns (rows) of ``\\mathbf{S}`` (``\\mathbf{S}^{-1}``) in the Smith normal decomposition
of the EBR matrix ``\\mathbf{A} = \\mathbf{S}\\boldsymbol{\\Lambda}\\mathbf{T}``.

## Examples

```julia-repl
julia> brs = calc_bandreps(22, Val(3)); # from Crystalline.jl
julia> n = parse(SymmetryVector, "Z₃, T₃, L₁, Y₃, Γ₃", irreps(brs)) # a compatible vector

# test a compatible symmetry vector
julia> iscompatible(n, brs)
true

# test an invalid symmetry vector
julia> n′ = copy(n);
julia>  n′.multsv[1] .= [1,0,0,0]  # change Z₃ to Z₁; incompatible modification
julia> iscompatible(n′, brs)
false

# test a symmetry vector with negative content
julia> n′′ = brs[1] + brs[2] - brs[3];  # contains negative elements
julia> iscompatible(n′′, brs)
false
julia> iscompatible(n′′, brs; allow_negative=true)
true
```
"""
function iscompatible(
    n::AbstractVector{<:Integer},
    F::Smith;
    allow_negative::Bool=false
)
    allow_negative || all(≥(0), n) || return false # check non-negativity

    # check compatibility relations
    dᵇˢ = count(!iszero, F.SNF)
    S̃   = @view F.S[:,OneTo(dᵇˢ)]     # relevant columns of S only
    S̃⁻¹ = @view F.Sinv[OneTo(dᵇˢ), :] # relevant rows of S⁻¹ only
    return S̃*(S̃⁻¹*n) == n
end
function iscompatible(n::AbstractVector{<:Integer}, B::Matrix{<:Integer}; kws...)
    iscompatible(n, smith(B); kws...)
end
function iscompatible(
    n::AbstractVector{<:Integer}, 
    brs::Union{BandRepSet, Collection{<:NewBandRep}};
    kws...
)
    iscompatible(n, stack(brs); kws...)
end

# ---------------------------------------------------------------------------------------- #
# Utilities/helper functions

function _throw_incompatible_or_negative(n)
    error(DomainError(n, "`n` is not a physically realizable band grouping"))
end

"""
    $(TYPEDSIGNATURES)

Return whether `n` includes the connectivity as an element by comparing with size of `brs`.
"""
function includes_connectivity(
    n::AbstractVector{<:Integer},
    brs::Union{Collection{<:NewBandRep}, BandRepSet}
)
    Nn = length(n)
    Nirr = brs isa Collection{<:NewBandRep} ? length(first(brs))-1 : length(irreplabels(brs))
    if Nn == Nirr+1
        return true
    elseif Nn == Nirr
        return false
    else 
        error(DimensionMismatch("incompatible dimensions of `n` and `brs`"))
    end
end