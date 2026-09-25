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
# Trivial/nontrivial solution topology via Smith/Collection{<:BandRep}

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

The EBR basis can be provided as `::Collection{<:BandRep}`, `::Matrix{<:Integer}`,
or a `Smith`
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
    
    (; Λ, S̃⁻¹) = smith_column_bases(F) # Λ = [λ₁, …, λ_{dᵇˢ}], S̃⁻¹ = S⁻¹[1:dᵇˢ, :]

    # n is trivial if (S⁻¹n)ⱼ = 0 mod λⱼ for j = 1, …, dᵇˢ. This is equivalent to checking
    # whether there exists an integer coefficient expansion for `n` in the EBR basis that
    # `F` represents (i.e., whether `Nemo.cansolve(B, n) == true`) but faster.
    # We do the matrix-vector product row-wise to check `mod((S̃⁻¹*n)[i], Λ[i]) = 0` for
    # `i ∈ 1:dᵇˢ` without allocating unnecessarily
    is_trivial = all(eachindex(Λ)) do i
        Λᵢ = Λ[i]
        Λᵢ == 1 && return true # fast path: `mod(x, 1) = 0` for all integer `x`.
        mod(dot(@view(S̃⁻¹[i,:]), n), Λᵢ) == 0
    end
    return is_trivial ? TRIVIAL : NONTRIVIAL
end

function calc_topology(n::AbstractVector{<:Integer}, B::AbstractMatrix{<:Integer}; kws...)
    length(n)≠size(B, 1) && throw(DimensionMismatch("incompatible dimensions of `n` & `B`"))
    return calc_topology(n, smith(B); kws...)
end

function calc_topology(
    n::AbstractVector{<:Integer},
    brs::Collection{<:BandRep};
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
elementary band representations (EBRs) `brs`, provided as a `Collection{<:BandRep}`, a
`Matrix{<:Integer}`, or a `Smith` decomposition thereof.

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
    brs::Collection{<:BandRep};
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
julia> brs = bandreps(2, Val(3));

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
function indicator_group(brs::AbstractVector{<:AbstractVector{<:Integer}})
    return indicator_group(stack(brs))
end
is_not_one_or_zero(x) = !(isone(x) || iszero(x))

"""
    basisdim(brs::Collection{<:BandRep})   --> Int
    basisdim(B::AbstractMatrix{<:Integer}) --> Int
    basisdim(F::Smith)                     --> Int

Return the dimension of the (linearly independent parts) of a band representation basis.
This is ``d^{\\text{bs}} = d^{\\text{ai}}`` in the notation of [Po, Watanabe, & Vishwanath,
Nature Commun. **8**, 50 (2017)](https://doi.org/10.1038/s41467-017-00133-2), or 
equivalently, the rank of `stack(brs)` over the ring of integers.
This is the number of linearly independent basis vectors that span the expansions of
a band structure viewed as symmetry data.
""" 
basisdim(F::Smith) = count(!iszero, F.SNF) # nonzeros of the Smith normal diagonal matrix
basisdim(B::AbstractMatrix{<:Integer}) = basisdim(smith(B, inverse=false))
basisdim(brs::AbstractVector{<:AbstractVector{<:Integer}}) = basisdim(stack(brs))

"""
    smith_column_bases(brs::Collection{<:BandRep})   --> @NamedTuple{S̃, Λ, S̃⁻¹}
    smith_column_bases(B::AbstractMatrix{<:Integer}) --> @NamedTuple{S̃, Λ, S̃⁻¹}
    smith_column_bases(F::Smith)                     --> @NamedTuple{S̃, Λ, S̃⁻¹}

Return the parts of the Smith normal decomposition of a set of band representations `brs`
(or of its matrix `B`, or of a `Smith` decomposition `F` thereof) that pertain to the column
space of the band representation matrix, i.e. the parts associated with the ``d^{\\text{bs}}
= `` [`basisdim`](@ref)`(brs)` nonzero elementary factors of ``\\boldsymbol{\\Lambda}``:

- `S̃`: the first ``d^{\\text{bs}}`` **columns** of ``\\mathbf{S}``; an integer-coefficient
  basis for all gapped band structures {BS}.
- `Λ`: the first ``d^{\\text{bs}}`` elements of ``\\boldsymbol{\\Lambda}``, i.e. its nonzero
  elementary factors ``\\lambda_1, \\ldots, \\lambda_{d^{\\text{bs}}}``.
- `S̃⁻¹`: the first ``d^{\\text{bs}}`` **rows** of ``\\mathbf{S}^{-1}``; these take a
  symmetry vector `n` to its coefficients in `S̃`, i.e. `S̃⁻¹*n`.

Note that `S̃⁻¹` is a slice of ``\\mathbf{S}^{-1}``, *not* the inverse of the (generally
nonsquare) `S̃`; the two nevertheless satisfy `S̃⁻¹*S̃ == I`. All three returned quantities
are views into `F`, so nothing is allocated.

A basis for the atomic insulators {AI} — the bands induced by localized orbitals at the
Wyckoff positions — is `S̃*Diagonal(Λ)`. A symmetry vector that can be expanded on that
basis with positive integer coefficients is a trivial insulator (i.e., deformable to an
atomic limit); one that cannot is topological, either fragilely (some negative coefficients)
or strongly (fractional coefficients). [`calc_topology`](@ref) distinguishes the strong case
from `Λ` and `S̃⁻¹` alone.

## Implementation

For an n×m integer matrix ``\\mathbf{B}``, the Smith normal form gives integer matrices
``\\mathbf{S}``, ``\\mathrm{diagm}(\\boldsymbol{\\Lambda})`` and ``\\mathbf{T}`` (of size
n×n, n×m and m×m, respectively) with ``\\mathbf{B} =
\\mathbf{S}\\mathrm{diagm}(\\boldsymbol{\\Lambda})\\mathbf{T}``, where
``\\boldsymbol{\\Lambda} = [\\lambda_1, \\ldots, \\lambda_r, 0, \\ldots, 0]`` with
``\\lambda_{j+1}`` divisible by ``\\lambda_j`` and ``r = d^{\\text{bs}} \\leq \\min(n,m)``;
``\\mathbf{S}`` and ``\\mathbf{T}`` have integer-valued inverses.

Applying `S̃⁻¹` to integer symmetry data ``\\mathbf{n}`` gives the integer factors
``q_i C_i`` (``C_i = \\lambda_i`` here) of [Tang, Po, Vishwanath, & Wan, Nature Physics
**15**, 470 (2019)](https://doi.org/10.1038/s41567-019-0418-7).
"""
function smith_column_bases(F::Smith)
    nzidxs = OneTo(basisdim(F))
    return (; S̃   = @view(F.S[:, nzidxs]),     # relevant columns of S only
              Λ   = @view(F.SNF[nzidxs]),      # nonzero elementary factors only
              S̃⁻¹ = @view(F.Sinv[nzidxs, :]))  # relevant rows of S⁻¹ only
end
smith_column_bases(B::AbstractMatrix{<:Integer}) = smith_column_bases(smith(B))
function smith_column_bases(brs::AbstractVector{<:AbstractVector{<:Integer}})
    return smith_column_bases(stack(brs))
end


@doc """
$(TYPEDSIGNATURES)

Return the symmetry indicator group ``X^{\\text{BS}}`` as a formatted string (i.e., 
as `"Zᵢ×Zⱼ×…"`). See also [`indicator_group`](@ref) for a vector representation.

## Example
```jldoctest
julia> brs = bandreps(2, Val(3));

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
    brs::Union{AbstractVector{<:AbstractVector{<:Integer}},
               AbstractMatrix{<:Integer}, Smith}
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
provided either as a `Collection{<:BandRep}`, a `Matrix{<:Integer}`, or a `Smith`
decomposition. The irrep sorting of `n` and this set of EBRs must be identical.

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
julia> brs = bandreps(22, Val(3)); # from Crystalline.jl
julia> n = parse(SymmetryVector, "Z₃, T₃, L₁, Y₃, Γ₃", irreps(brs)) # a compatible vector

# test a compatible symmetry vector
julia> iscompatible(n, brs)
true

# test an invalid symmetry vector
julia> n′ = copy(n);
julia> multiplicities(n′)[1] .= [1,0,0,0]  # change Z₃ to Z₁; incompatible modification
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
    (; S̃, S̃⁻¹) = smith_column_bases(F)
    return S̃*(S̃⁻¹*n) == n
end
function iscompatible(n::AbstractVector{<:Integer}, B::Matrix{<:Integer}; kws...)
    iscompatible(n, smith(B); kws...)
end
function iscompatible(
    n::AbstractVector{<:Integer}, 
    brs::Collection{<:BandRep};
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
    brs::Collection{<:BandRep}
)
    Nn = length(n)
    Nirr = length(first(brs))-1
    if Nn == Nirr+1
        return true
    elseif Nn == Nirr
        return false
    else 
        error(DimensionMismatch("incompatible dimensions of `n` and `brs`"))
    end
end