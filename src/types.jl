# ---------------------------------------------------------------------------------------- #
# SymOperation
# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct SymOperation{D} <: AbstractMatrix{Float64}
    rotation    :: SqSMatrix{D, Float64} # a square stack-allocated matrix
    translation :: SVector{D, Float64}
end
SymOperation(m::AbstractMatrix{<:Real}) = SymOperation{size(m,1)}(float(m))
function SymOperation{D}(m::AbstractMatrix{Float64}) where D
    @boundscheck size(m) == (D, D+1) || throw(DomainError(size(m), "matrix size must be (D, D+1)"))
    # we construct the SqSMatrix's tuple data manually, rather than calling e.g.
    # `convert(SqSMatrix{D, Float64}, @view m[OneTo(D),OneTo(D)])`, just because it is
    # slightly faster that way (maybe the view has a small penalty?)...
    tuple_cols  = ntuple(j -> ntuple(i -> (@inbounds m[i,j]), Val(D)), Val(D))
    rotation    = SqSMatrix{D, Float64}(tuple_cols)
    translation = SVector{D, Float64}(ntuple(j -> (@inbounds m[j, D+1]), Val(D)))
    SymOperation{D}(rotation, translation)
end
function SymOperation(r::Union{SMatrix{D,D,<:Real}, MMatrix{D,D,<:Real}}, 
                      t::Union{SVector{D,<:Real}, MVector{D,<:Real}}=zero(SVector{D,Float64})
                      ) where D
    SymOperation{D}(SqSMatrix{D,Float64}(r), t)
end
SymOperation(t::SVector{D,<:Real}) where D = SymOperation(one(SqSMatrix{D,Float64}), SVector{D,Float64}(t))
SymOperation{D}(t::AbstractVector{<:Real}) where D = SymOperation(one(SqSMatrix{D,Float64}), SVector{D,Float64}(t))
# extracting StaticArray representations of the symmetry operation, amenable to linear algebra
rotation(op::SymOperation{D}) where D = SMatrix(op.rotation)
translation(op::SymOperation{D}) where D = op.translation
matrix(op::SymOperation{D}) where D = 
    SMatrix{D, D+1, Float64, D*(D+1)}((SquareStaticMatrices.flatten(op.rotation)..., 
                                       translation(op)...))

# string constructors
xyzt(op::SymOperation) = matrix2xyzt(matrix(op))
SymOperation{D}(s::AbstractString) where D = SymOperation{D}(xyzt2components(s, Val(D))...)
# type-unstable convenience constructors; avoid for anything non-REPL related, if possible
SymOperation(m::Matrix{<:Real}) = SymOperation{size(m,1)}(float(m))
SymOperation(t::AbstractVector{<:Real}) = SymOperation{length(t)}(t)
SymOperation(s::AbstractString) = SymOperation{count(==(','), s)+1}(s)

# define the AbstractArray interface for SymOperation
@propagate_inbounds getindex(op::SymOperation, i::Int) = matrix(op)[i]
IndexStyle(::Type{<:SymOperation}) = IndexLinear()
size(::SymOperation{D}) where D = (D, D+1)
copy(op::SymOperation) = op # cf. https://github.com/JuliaLang/julia/issues/41918

rotation(m::AbstractMatrix{<:Real}) = @view m[:,1:end-1] # rotational (proper or improper) part of an operation
translation(m::AbstractMatrix{<:Real}) = @view m[:,end]  # translation part of an operation
rotation(m::SMatrix{D,Dp1,<:Real}) where {D,Dp1} = m[:,SOneTo(D)] # needed for type-stability w/ StaticArrays (returns an SMatrix{D,D,...})
translation(m::SMatrix{D,Dp1,<:Real}) where {D,Dp1} = m[:,Dp1]    # not strictly needed for type-stability    (returns an SVector{D,...})

dim(::SymOperation{D}) where D = D
function (==)(op1::SymOperation, op2::SymOperation)
    if dim(op1) == dim(op2) && op1.rotation == op2.rotation && translation(op1) == translation(op2)
        return true
    else
        return false
    end
end
isapprox(op1::SymOperation, op2::SymOperation; kwargs...)= (dim(op1) == dim(op2) && isapprox(matrix(op1), matrix(op2); kwargs...))
unpack(op::SymOperation) = (rotation(op), translation(op))

one(::Type{SymOperation{D}}) where D = SymOperation{D}(one(SqSMatrix{D,Float64}),
                                                       zero(SVector{D,Float64}))
function isone(op::SymOperation{D}) where D
    @inbounds for D₁ in SOneTo(D)
        iszero(op.translation[D₁]) || return false
        for D₂ in SOneTo(D)
            check = D₁ ≠ D₂ ? iszero(op.rotation[D₁,D₂]) : isone(op.rotation[D₁,D₁])
            check || return false
        end
    end
    return true
end

# ---------------------------------------------------------------------------------------- #
# MultTable
# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct MultTable{D} <: AbstractMatrix{SymOperation{D}}
    operations::Vector{SymOperation{D}}
    table::Matrix{Int64} # Cayley table: indexes into `operations`
end
@propagate_inbounds function getindex(mt::MultTable, i::Int)
    mtidx = mt.table[i]
    return mt.operations[mtidx]
end
size(mt::MultTable) = size(mt.table)
IndexStyle(::Type{<:MultTable}) = IndexLinear()


# ---------------------------------------------------------------------------------------- #
# AbstractVec: KVec & RVec
# ---------------------------------------------------------------------------------------- #

# 𝐤-vectors (or, equivalently, 𝐫-vectors) are specified as a pair (k₀, kabc), denoting a 𝐤-vector
#       𝐤 = ∑³ᵢ₌₁ (k₀ᵢ + aᵢα+bᵢβ+cᵢγ)*𝐆ᵢ     (w/ recip. basis vecs. 𝐆ᵢ)
# here the matrix kabc is columns of the vectors (𝐚,𝐛,𝐜) while α,β,γ are free parameters
# ranging over all non-special values (i.e. not coinciding with any high-sym 𝐤).

abstract type AbstractVec{D} end
# A type which must have a vector field `cnst` (subtyping `AbstractVector`, of length `D`)
# and a matrix field `free` (subtyping `AbstractMatrix`; of size `(D,D)`).
# Intended to represent points, lines, planes and volumes in direct (`RVec`) or reciprocal
# space (`KVec`).
constant(v::AbstractVec)  = v.cnst
free(v::AbstractVec)      = v.free
parts(v::AbstractVec)     = (constant(v), free(v))
dim(::AbstractVec{D}) where D = D
isspecial(v::AbstractVec) = iszero(free(v))
"""
$(TYPEDSIGNATURES)

Return a vector whose entries are `true` (`false`) if the free parameters α,β,γ, 
respectively, occur with nonzero (zero) coefficients in `v`.
"""
freeparams(v::AbstractVec)  = map(colⱼ->!iszero(colⱼ), eachcol(free(v)))
"""
$(TYPEDSIGNATURES)

Return total number of free parameters occuring in `v`.
"""
nfreeparams(v::AbstractVec) = count(colⱼ->!iszero(colⱼ), eachcol(free(v)))

function (v::AbstractVec)(αβγ::AbstractVector{<:Real})
    cnst, free = parts(v)
    return cnst + free*αβγ
end
(v::AbstractVec)(αβγ::Vararg{<:Real}) = v([αβγ...])
(v::AbstractVec)(::Nothing) = constant(v)
(v::AbstractVec)()          = v(nothing)


# parsing `AbstractVec`s from string format
function _strip_split(str::AbstractString)
    str = filter(!isspace, strip(str, ['(',')','[',']'])) # tidy up string (remove parens & spaces)
    return split(str,',')
end
function parse_abstractvec(xyz::Vector{<:SubString}, T::Type{<:AbstractVec{D}}) where D
    length(xyz) == D || throw(DimensionMismatch("Dimension D doesn't match input string"))
    cnst = zero(MVector{D, Float64})
    free = zero(MMatrix{D, D, Float64})
    for (i, coord) in enumerate(xyz)
        # --- "free" coordinates, free[i,:] ---
        for (j, matchgroup) in enumerate((('α','u','x'),('β','v','y'),('γ','w','z')))
            pos₂ = findfirst(∈(matchgroup), coord)
            if !isnothing(pos₂)
                free[i,j]  = searchpriornumerals(coord, pos₂)
            end
        end
        
        # --- "fixed"/constant coordinate, cnst[i] ---
        m = match(r"(?:\+|\-)?(?:(?:[0-9]|/|\.)+)(?!(?:[0-9]|\.)*[αuxβvyγwz])", coord)
        # regex matches any digit sequence, possibly including slashes, that is _not_
        # followed by one of the free-part identifiers αuβvγw (this is the '(?!' bit). 
        # If a '+' or '-' exist before the first digit, it is included in the match. 
        # The '(?:' bits in the groups simply makes sure that we don't actually create a
        # capture group, because we only need the match and not the individual captures 
        # (i.e. just a small optimization of the regex).
        # We do not allow arithmetic aside from division here, obviously: any extra numbers 
        # terms are ignored.
        if m===nothing   # no constant terms
            if last(coord) ∈ ('α','u','x','β','v','y','γ','w','z') # free-part only case
                continue # cnst[i] is zero already
            else
                throw(ErrorException("Unexpected parsing error in constant term"))
            end
        else
            cnst[i] = parsefraction(m.match)
        end
    end
    return T(cnst, free)
end

# generate KVec and RVec structs and parsers jointly...
for T in (:KVec, :RVec)
    @eval begin
    struct $T{D} <: AbstractVec{D}
        cnst :: SVector{D,Float64}
        free :: SqSMatrix{D,Float64}
    end
    free(v::$T) = SMatrix(v.free)

    @doc """
        $($T){D}(str::AbstractString) --> $($T){D}
        $($T)(str::AbstractString)    --> $($T)
        $($T)(::AbstractVector, ::AbstractMatrix) --> $($T)
    
    Return a `$($T)` by parsing the string representations `str`, supplied in one of the
    following formats:

    ```jl
    "(\$1,\$2,\$3)"
    "[\$1,\$2,\$3]"
    "\$1,\$2,\$3"
    ```

    where the coordinates `\$1`,`\$2`, and `\$3` are strings that may contain fractions,
    decimal numbers, and "free" parameters {`'α'`,`'β'`,`'γ'`} (or, alternatively and
    equivalently, {`'u'`,`'v'`,`'w'`} or {`'x'`,`'y'`,`'z'`}).
    
    Fractions such as `1/2` and decimal numbers can be parsed: but use of any other special
    operator besides `/` will produce undefined behavior (e.g. do not use `*`).

    ## Example
    ```jldoctest
    julia> $($T)("0.25,α,0")
    [0.25, α, 0.0]
    ```
    """
    function $T{D}(str::AbstractString) where D
        xyz = _strip_split(str)
        return parse_abstractvec(xyz, $T{D})
    end
    function $T(str::AbstractString)
        xyz = _strip_split(str)
        D   = length(xyz)
        return parse_abstractvec(xyz, $T{D})
    end
    function $T(cnst::AbstractVector{<:Real}, 
                free::AbstractMatrix{<:Real}=zero(SqSMatrix{length(cnst), Float64}))
        D = length(cnst)
        @boundscheck D == LinearAlgebra.checksquare(free) || throw(DimensionMismatch("Mismatched argument sizes"))
        $T{D}(SVector{D,Float64}(cnst), SqSMatrix{D,Float64}(free))
    end
    $T(xs::Vararg{<:Real, D}) where D = $T(SVector{D, Float64}(xs))
    $T(xs::NTuple{D, <:Real}) where D = $T(SVector{D, Float64}(xs))
    end
end

# arithmetic with abstract vectors
(-)(v::T) where T<:AbstractVec = T(-constant(v), -free(v))
for op in (:(-), :(+))
    @eval function $op(v1::T, v2::T) where T<:AbstractVec
        cnst1, free1 = parts(v1); cnst2, free2 = parts(v2) 
        return T($op(cnst1, cnst2), $op(free1, free2))
    end
    @eval function $op(v1::T, cnst2::AbstractVector) where T<:AbstractVec
        dim(v1) == length(cnst2) || throw(DimensionMismatch("argument dimensions must be equal"))
        cnst1, free1 = parts(v1)
        return T($op(cnst1, cnst2), free1)
    end
    @eval function $op(cnst1::AbstractVector, v2::T) where T<:AbstractVec
        length(cnst1) == dim(v2) || throw(DimensionMismatch("argument dimensions must be equal"))
        cnst2, free2 = parts(v1)
        return T($op(cnst1, cnst2), $op(free2))
    end
end
function zero(::Type{T}) where T<:AbstractVec{D} where D
    T(zero(SVector{D,Float64}), zero(SqSMatrix{D,Float64}))
end
zero(v::AbstractVec) = zero(typeof(v))

# `isapprox` without considerations of lattice-vectors
function isapprox(v1::AbstractVec, v2::AbstractVec; kwargs...)
    cnst1, free1 = parts(v1); cnst2, free2 = parts(v2)  # ... unpacking
       
    return isapprox(cnst1, cnst2; kwargs...) && isapprox(free1, free2; kwargs...)
end
function (==)(v1::AbstractVec, v2::AbstractVec)
    cnst1, free1 = parts(v1); cnst2, free2 = parts(v2)  # ... unpacking
       
    return cnst1 == cnst2 && free1 == free2
end

"""
    isapprox(kv1::KVec, kv2::KVec[, cntr::Char]; kwargs...) --> Bool
                                                            
Compute approximate equality of two `KVec`s `k1` and `k2` modulo any primitive reciprocal
lattice vectors. To ensure that primitive reciprocal lattice vectors are used, the centering
type `cntr` (see [`Bravais.centering`](@ref)) must be given (the dimensionality is inferred
from `kv1` and `kv2`).
Optionally, keyword arguments (e.g., `atol` and `rtol`) can be provided, to be forwarded
to `Base.isapprox`.

If `cntr` is not provided, the comparison will not account for equivalence by primitive
reciprocal lattice vectors.
"""
function isapprox(kv1::KVec{D}, kv2::KVec{D}, cntr::Char; kwargs...) where D
    k₀1, kabc1 = parts(kv1); k₀2, kabc2 = parts(kv2)  # ... unpacking

    # check if k₀ ≈ k₀′ differ by a _primitive_ 𝐆 vector
    diff = primitivebasismatrix(cntr, Val(D))' * (k₀1 .- k₀2)
    kbool = all(el -> isapprox(el, round(el); kwargs...), diff) 
    # check if kabc1 ≈ kabc2; no need to check for difference by a 
    # 𝐆 vector, since kabc is in interior of BZ
    abcbool = isapprox(kabc1, kabc2; kwargs...)

    return kbool && abcbool
end


# Note that the _coefficients_ of a general 𝐤- or 𝐫-vector transforms differently than the
# reciprocal or direct _basis_, which transforms from non-primed to primed variants. See
# discussion in Bravais.jl /src/transform.jl.
@doc raw"""
    transform(v::AbstractVec, P::AbstractMatrix)  -->  v′::typeof(v)

Return a transformed coordinate vector `v′` from an original coordinate vector `v` using a
basis change matrix `P`.

Note that a basis change matrix ``\mathbf{P}`` transforms direct coordinate vectors (`RVec`)
as ``\mathbf{r}' = \mathbf{P}^{-1}\mathbf{r}`` but transforms reciprocal coordinates
(`KVec`) as ``\mathbf{k}' = \mathbf{P}^{\mathrm{T}}\mathbf{k}`` [^ITA6]

[^ITA6]: M.I. Aroyo, International Tables of Crystallography, Vol. A, 6th edition (2016):
         Secs. 1.5.1.2 and 1.5.2.1.
"""
transform(::AbstractVec, ::AbstractMatrix{<:Real})

function transform(kv::KVec{D}, P::AbstractMatrix{<:Real}) where D
    # P maps an "original" reciprocal-space coefficient vector (k₁ k₂ k₃)ᵀ to a transformed
    # coefficient vector (k₁′ k₂′ k₃′)ᵀ = Pᵀ(k₁ k₂ k₃)ᵀ
    k₀, kabc = parts(kv)
    return KVec{D}(P'*k₀, P'*kabc)
end
function transform(rv::RVec{D}, P::AbstractMatrix{<:Real}) where D
    # P maps an "original" direct-space coefficient vector (r₁ r₂ r₃)ᵀ to a transformed
    # coefficient vector (r₁′ r₂′ r₃′)ᵀ = P⁻¹(r₁ r₂ r₃)ᵀ
    rcnst, rfree = parts(rv)
    return RVec{D}(P\rcnst, P\rfree)
end

@doc raw"""
    primitivize(v::AbstractVec, cntr::Char)  -->  v′::typeof(v)

Transforms a conventional coordinate vector `v` to a standard primitive basis (specified by
the centering type `cntr`), returning the primitive coordinate vector `v′`.

Note that a basis change matrix ``\mathbf{P}`` (as returned e.g. by
[`Bravais.primitivebasismatrix`](@ref)) transforms direct coordinate vectors
([`RVec`](@ref)) as ``\mathbf{r}' = \mathbf{P}^{-1}\mathbf{r}`` but transforms reciprocal
coordinates ([`KVec`](@ref)) as ``\mathbf{k}' = \mathbf{P}^{\text{T}}\mathbf{k}`` [^ITA6].
Recall also the distinction between transforming a basis and the coordinates of a vector.

[^ITA6]: M.I. Aroyo, International Tables of Crystallography, Vol. A, 6th edition (2016):
         Secs. 1.5.1.2 and 1.5.2.1.
"""
function primitivize(v::AbstractVec{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p'
        return v
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(v, P)
    end
end

"""
    conventionalize(v′::AbstractVec, cntr::Char)  -->  v::typeof(v′)

Transforms a primitive coordinate vector `v′` back to a standard conventional basis
(specified by the centering type `cntr`), returning the conventional coordinate vector `v`.

See also [`primitivize`](@ref) and [`transform`](@ref).
"""
function conventionalize(v′::AbstractVec{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p'
        return v′
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(v′, inv(P))
    end
end


# ---------------------------------------------------------------------------------------- #
# AbstractGroup: Generic Group, SpaceGroup, PointGroup, LittleGroup
# ---------------------------------------------------------------------------------------- #

abstract type AbstractGroup{D} <: AbstractVector{SymOperation{D}} end
num(g::AbstractGroup) = g.num
operations(g::AbstractGroup) = g.operations
dim(::AbstractGroup{D}) where D = D
order(g::AbstractGroup) = length(g)

# define the AbstractArray interface for AbstractGroup
@propagate_inbounds getindex(g::AbstractGroup, i::Int) = operations(g)[i]    # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
@propagate_inbounds setindex!(g::AbstractGroup, op::SymOperation, i::Int) = (operations(g)[i] .= op)
size(g::AbstractGroup) = size(operations(g))
IndexStyle(::Type{<:AbstractGroup}) = IndexLinear()

# --- Generic group ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct GenericGroup{D} <: AbstractGroup{D}
    operations::Vector{SymOperation{D}}
end
num(::GenericGroup) = 0
label(::GenericGroup) = ""

# --- Space group ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct SpaceGroup{D} <: AbstractGroup{D}
    num::Int64
    operations::Vector{SymOperation{D}}
end
label(sg::SpaceGroup) = iuc(sg)

# --- Point group ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct PointGroup{D} <: AbstractGroup{D}
    num::Int64
    label::String
    operations::Vector{SymOperation{D}}
end
label(pg::PointGroup) = pg.label
iuc(pg::PointGroup) = label(pg)

# --- Little group ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct LittleGroup{D} <: AbstractGroup{D}
    num::Int64
    kv::KVec{D}
    klab::String
    operations::Vector{SymOperation{D}}
end
LittleGroup(num::Int64, kv::KVec, klab::String, ops::AbstractVector{SymOperation{D}}) where D = LittleGroup{D}(num, kv, klab, ops)
LittleGroup(num::Int64, kv::KVec, ops::AbstractVector{SymOperation{D}}) where D = LittleGroup{D}(num, kv, "", ops)
kvec(lg::LittleGroup) = lg.kv
klabel(lg::LittleGroup) = lg.klab
label(lg::LittleGroup)  = iuc(num(lg), dim(lg))*" at "*klabel(lg)*" = "*string(kvec(lg))

# ---------------------------------------------------------------------------------------- #
# Reality <: Enum{Int8}:    1 (≡ real), -1 (≡ pseudoreal), 0 (≡ complex)
# ---------------------------------------------------------------------------------------- #
 
"""
    Reality <: Enum{Int8}

Enum type with instances

```  
REAL = 1
PSEUDOREAL = -1
COMPLEX = 0
```

The return value of [`reality(::AbstractIrrep)`](@ref) and [`calc_reality`](@ref) is an
instance of `Reality`. The reality type of an irrep is relevant for constructing "physically
real" irreps (co-reps) via [`realify`](@ref).
"""
@enum Reality::Int8 begin
    REAL       = 1
    PSEUDOREAL = -1
    COMPLEX    = 0
end

# ---------------------------------------------------------------------------------------- #
# AbstractIrrep: PGIrrep, LGIrrep
# ---------------------------------------------------------------------------------------- #

""" 
    AbstractIrrep{D}

Abstract supertype for irreps of dimensionality `D`: must have fields `cdml`, `matrices`,
`g` (underlying group), `reality` (and possibly `translations`). May overload a function
`irreps` that returns the associated irrep matrices; if not, will simply be `matrices`.
"""
abstract type AbstractIrrep{D} end
(ir::AbstractIrrep)(αβγ=nothing) = ir.matrices
group(ir::AbstractIrrep, αβγ=nothing) = ir.g
label(ir::AbstractIrrep) = ir.cdml
matrices(ir::AbstractIrrep) = ir.matrices    
reality(ir::AbstractIrrep) = ir.reality
translations(ir::T) where T<:AbstractIrrep = hasfield(T, :translations) ? ir.translations : nothing
characters(ir::AbstractIrrep, αβγ::Union{AbstractVector{<:Real},Nothing}=nothing) = tr.(ir(αβγ))
irdim(ir::AbstractIrrep)  = size(first(matrices(ir)),1)
klabel(ir::AbstractIrrep) = klabel(label(ir))
order(ir::AbstractIrrep)  = order(group(ir))
operations(ir::AbstractIrrep) = operations(group(ir))
num(ir::AbstractIrrep) = num(group(ir))
dim(::AbstractIrrep{D}) where D = D
function klabel(cdml::String)
    idx = findfirst(c->isdigit(c) || issubdigit(c) || c=='ˢ', cdml) # look for regular digit or subscript digit
    previdx = idx !== nothing ? prevind(cdml, idx) : lastindex(cdml)
    return cdml[firstindex(cdml):previdx]
end
"""
    $TYPEDSIGNATURES --> Bool

Return whether the provided irrep has been made "physically real" (i.e. is a corep) so that
it differs from the underlying irrep (i.e. whether the irrep and "derived" corep differ).

For an irrep that has not been passed to [`realify`](@ref), this is always `false`.
For an irrep produced by `realify`, this can be either `false` or `true`: if the reality
type is `REAL` it is `false`; if the reality type is `PSEUDOREAL` or `COMPLEX` it is `true`.
"""
iscorep(ir::AbstractIrrep) = ir.iscorep

# --- Point group irreps ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct PGIrrep{D} <: AbstractIrrep{D}
    cdml::String
    g::PointGroup{D}
    matrices::Vector{Matrix{ComplexF64}}
    reality::Reality
    iscorep::Bool
end
function PGIrrep{D}(cdml::String, pg::PointGroup{D}, matrices::Vector{Matrix{ComplexF64}},
                    reality::Reality) where D 
    PGIrrep{D}(cdml, pg, matrices, reality, false)
end

# printing
function prettyprint_irrep_matrix(io::IO, pgir::PGIrrep, i::Integer, prefix::AbstractString)
    P = pgir.matrices[i]
    prettyprint_scalar_or_matrix(io, P, prefix, false)
end

# --- Little group irreps ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct LGIrrep{D} <: AbstractIrrep{D}
    cdml::String # CDML label of irrep (including k-point label)
    g::LittleGroup{D} # contains sgnum, kvec, klab, and operations that define the little group
    matrices::Vector{Matrix{ComplexF64}}
    translations::Vector{Vector{Float64}}
    reality::Reality
    iscorep::Bool
end
function LGIrrep{D}(cdml::String, lg::LittleGroup{D}, 
                    matrices::Vector{Matrix{ComplexF64}}, 
                    translations::Vector{Vector{Float64}},
                    reality::Reality) where D
    return LGIrrep{D}(cdml, lg, matrices, translations, reality, false)
end
function LGIrrep{D}(cdml::String, lg::LittleGroup{D}, 
                    matrices::Vector{Matrix{ComplexF64}}, 
                    translations_sentinel::Nothing, # sentinel value for all-zero translations
                    reality::Reality) where D
    translations = [zeros(Float64,D) for _=OneTo(order(lg))]
    return LGIrrep{D}(cdml, lg, matrices, translations, reality)
end
kvec(lgir::LGIrrep)  = kvec(group(lgir))
isspecial(lgir::LGIrrep)  = isspecial(kvec(lgir))
issymmorph(lgir::LGIrrep) = issymmorph(group(lgir))
kstar(lgir::LGIrrep) = kstar(spacegroup(num(lgir), dim(lgir)), 
                             kvec(lgir), centering(num(lgir), dim(lgir)))
function (lgir::LGIrrep)(αβγ::Union{Vector{<:Real},Nothing}=nothing)
    P = lgir.matrices
    τ = lgir.translations
    if !iszero(τ)
        k = kvec(lgir)(αβγ)
        P = deepcopy(P) # needs deepcopy rather than a copy due to nesting; otherwise we overwrite..!
        for (i,τ′) in enumerate(τ)
            if !iszero(τ′) && !iszero(k)
                P[i] .*= cis(2π*dot(k,τ′))  # note cis(x) = exp(ix)
                # NOTE/TODO/FIXME:
                # This follows the convention in Eq. (11.37) of Inui as well as the Bilbao
                # server, i.e. has Dᵏ({I|𝐭}) = exp(i𝐤⋅𝐭); but disagrees with several other
                # references (e.g. Herring 1937a and Kovalev's book; and even Bilbao's
                # own _publications_?!).
                # In these other references one take Dᵏ({I|𝐭}) = exp(-i𝐤⋅𝐭), while Inui takes
                # Dᵏ({I|𝐭}) = exp(i𝐤⋅𝐭) [cf. (11.36)]. The former choice, i.e. Dᵏ({I|𝐭}) =
                # exp(-i𝐤⋅𝐭) actually appears more natural, since we usually have symmetry 
                # operations acting _inversely_ on functions of spatial coordinates and
                # Bloch phases exp(i𝐤⋅𝐫).
                # Importantly, the exp(i𝐤⋅τ) is also the convention adopted by Stokes et al.
                # in Eq. (1) of Acta Cryst. A69, 388 (2013), i.e. in ISOTROPY (also
                # expliciated at https://stokes.byu.edu/iso/irtableshelp.php), so, overall,
                # this is probably the sanest choice for this dataset.
                # This weird state of affairs was also noted explicitly by Chen Fang in
                # https://doi.org/10.1088/1674-1056/28/8/087102 (near Eqs. (11-12)).
                #
                # If we wanted swap the sign here, we'd likely have to swap t₀ in the check
                # for ray-representations in `check_multtable_vs_ir(::MultTable, ::LGIrrep)`
                # to account for this difference. It is not enough just to swap the sign
                # - I checked (⇒ 172 failures in test/multtable.jl) - you would have 
                # to account for the fact that it would be -β⁻¹τ that appears in the 
                # inverse operation, not just τ. Same applies here, if you want to 
                # adopt the other convention, it should probably not just be a swap 
                # to -τ, but to -β⁻¹τ. Probably best to stick with Inui's definition.
            end
        end
    end
    # FIXME: Attempt to flip phase convention. Does not pass tests.
    #=
    lg = group(lgir)
    if !issymmorph(lg)
        k = kvec(lgir)(αβγ)
        for (i,op) in enumerate(lg)
            P[i] .* cis(-4π*dot(k, translation(op)))
        end
    end
    =#

    return P
end

"""
    israyrep(lgir::LGIrrep, αβγ=nothing) -> (::Bool, ::Matrix)

Computes whether a given little group irrep `ir` is a ray representation 
by computing the coefficients αᵢⱼ in DᵢDⱼ=αᵢⱼDₖ; if any αᵢⱼ differ 
from unity, we consider the little group irrep a ray representation
(as opposed to the simpler "vector" representations where DᵢDⱼ=Dₖ).
The function returns a boolean (true => ray representation) and the
coefficient matrix αᵢⱼ.
"""
function israyrep(lgir::LGIrrep, αβγ::Union{Nothing,Vector{Float64}}=nothing) 
    k = kvec(lgir)(αβγ)
    lg = group(lgir) # indexing into/iterating over `lg` yields the LittleGroup's operations
    Nₒₚ = length(lg)
    α = Matrix{ComplexF64}(undef, Nₒₚ, Nₒₚ)
    # TODO: Verify that this is OK; not sure if we can just use the primitive basis 
    #       here, given the tricks we then perform subsequently?
    mt = MultTable(primitivize(lg)) 
    for (row, oprow) in enumerate(lg)
        for (col, opcol) in enumerate(lg)
            t₀ = translation(oprow) + rotation(oprow)*translation(opcol) - translation(lg[mt.table[row,col]])
            ϕ  = 2π*dot(k,t₀) # include factor of 2π here due to normalized bases
            α[row,col] = cis(ϕ)
        end
    end
    return any(x->norm(x-1.0)>DEFAULT_ATOL, α), α
end


# ---------------------------------------------------------------------------------------- #
# CharacterTable
# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct CharacterTable{D} <: AbstractMatrix{ComplexF64}
    ops::Vector{SymOperation{D}}
    irlabs::Vector{String}
    chartable::Matrix{ComplexF64} # Stored as irreps-along-columns & operations-along-rows
    # TODO: for LGIrreps, it might be nice to keep this more versatile and include the 
    #       translations and kvec as well; then we could print a result that doesn't  
    #       specialize on a given αβγ choice (see also CharacterTable(::LGirrep))
    tag::String
end
function CharacterTable{D}(ops::AbstractVector{SymOperation{D}},
            irlabs::AbstractVector{String},
            chartable::AbstractMatrix{<:Real}) where D
    return CharacterTable{D}(ops, irlabs, chartable, "")
end

operations(ct::CharacterTable) = ct.ops
labels(ct::CharacterTable) = ct.irlabs
characters(ct::CharacterTable) = ct.chartable # TODO: maybe remove this? Redudant cf. `getindex`
tag(ct::CharacterTable) = ct.tag

@propagate_inbounds getindex(ct::CharacterTable, i::Int) = characters(ct)[i]
size(ct::CharacterTable) = size(ct.chartable)
IndexStyle(::Type{<:CharacterTable}) = IndexLinear()

"""
    CharacterTable(irs::AbstractVector{<:AbstractIrrep}, αβγ=nothing)

Returns a `CharacterTable` associated with vector of `AbstractIrrep`s `irs`. 

Optionally, an `αβγ::AbstractVector{<:Real}` variable can be passed to evaluate the irrep
(and associated characters) with concrete free parameters (e.g., for `LGIrrep`s, a concrete
k-vector sampled from a "line-irrep"). Defaults to `nothing`, indicating it being either 
irrelevant (e.g., for `PGIrrep`s) or all free parameters implicitly set to zero.
"""
function CharacterTable(irs::AbstractVector{<:AbstractIrrep{D}},
                        αβγ::Union{AbstractVector{<:Real}, Nothing}=nothing) where D
    table = Array{ComplexF64}(undef, order(first(irs)), length(irs))
    for (i,col) in enumerate(eachcol(table))
        col .= characters(irs[i], αβγ)
    end
    g = group(first(irs))
    tag = "#"*string(num(g))*" ("*label(g)*")"
    return CharacterTable{D}(operations(first(irs)), label.(irs), table, tag)
end

# ---------------------------------------------------------------------------------------- #
# BandRep & BandRepSet (band representations)
# ---------------------------------------------------------------------------------------- #

# --- BandRep ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct BandRep <: AbstractVector{Int}
    wyckpos::String  # Wyckoff position that induces the BR (TODO: type as `::WyckPos` instead of `::String`)
    sitesym::String  # Site-symmetry point group of Wyckoff pos (IUC notation)
    label::String    # Symbol ρ↑G, with ρ denoting the irrep of the site-symmetry group
    dim::Int         # Dimension (i.e. # of bands) in band rep
    decomposable::Bool  # Whether a given bandrep can be decomposed further
    spinful::Bool       # Whether a given bandrep involves spinful irreps ("\bar"'ed irreps)
    irvec::Vector{Int}  # Vector that references irlabs of a parent BandRepSet; nonzero
                           # entries correspond to an element in the band representation
    irlabs::Vector{String} # A reference to the labels; same as in the parent BandRepSet
end
wyck(BR::BandRep)    = BR.wyckpos
sitesym(BR::BandRep) = BR.sitesym
label(BR::BandRep)   = BR.label
irreplabels(BR::BandRep) = BR.irlabs

"""
    dim(BR::BandRep) --> Int

Return the number of bands included in the provided `BandRep`.

If the bands are "nondetachable" (i.e. if `BR.decomposable = false`), this is equal to a
band connectivity μ.
"""
dim(BR::BandRep)     = BR.dim

# define the AbstractArray interface for BandRep
size(BR::BandRep) = (size(BR.irvec)[1] + 1,) # number of irreps sampled by BandRep + 1 (filling)
@propagate_inbounds function getindex(BR::BandRep, i::Int)
    return i == length(BR.irvec)+1 ? dim(BR) : BR.irvec[i]
end
IndexStyle(::Type{<:BandRep}) = IndexLinear()
function iterate(BR::BandRep, i=1)
    # work-around performance issue noted in https://discourse.julialang.org/t/iteration-getindex-performance-of-abstractarray-wrapper-types/53729
    if i == length(BR)
        return dim(BR), i+1
    else
        return iterate(BR.irvec, i) # also handles `nothing` when iteration is done
    end
end

# --- BandRepSet ---
"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct BandRepSet <: AbstractVector{BandRep}
    sgnum::Int              # space group number, sequential
    bandreps::Vector{BandRep}
    kvs::Vector{<:KVec}     # Vector of 𝐤-points # TODO: Make parametric?
    klabs::Vector{String}   # Vector of associated 𝐤-labels (in CDML notation)
    irlabs::Vector{String}  # Vector of (sorted) CDML irrep labels at _all_ 𝐤-points
    allpaths::Bool          # Whether all paths (true) or only maximal 𝐤-points (false) are included
    spinful::Bool           # Whether the band rep set includes (true) or excludes (false) spinful irreps
    timeinvar::Bool         # Whether the band rep set assumes time-reversal symmetry (true) or not (false) 
end
num(BRS::BandRepSet)         = BRS.sgnum
klabels(BRS::BandRepSet)     = BRS.klabs
kvecs(BRS::BandRepSet)       = BRS.kvs
hasnonmax(BRS::BandRepSet)   = BRS.allpaths
irreplabels(BRS::BandRepSet) = BRS.irlabs
isspinful(BRS::BandRepSet)   = BRS.spinful
istimeinvar(BRS::BandRepSet) = BRS.timeinvar
reps(BRS::BandRepSet)        = BRS.bandreps

# define the AbstractArray interface for BandRepSet
size(BRS::BandRepSet) = (length(reps(BRS)),) # number of distinct band representations
@propagate_inbounds getindex(BRS::BandRepSet, i::Int) = reps(BRS)[i]
IndexStyle(::Type{<:BandRepSet}) = IndexLinear()

"""
    matrix(BRS::BandRepSet; includedim::Bool=true)

Return a matrix representation of `BRS::BandRepSet`, with band representations as columns 
and irreps over rows.

By default, the last row will give the "filling" of each `BandRep` (or, more precisely,
number of included bands per `BandRep`, i.e. `dim.(BRS)`. To toggle this off, set the
keyword argument `includedim` to `false` (default is `includedim = true`).
"""
function matrix(BRS::BandRepSet; includedim::Bool=true)
    Nⁱʳʳ, Nᵉᵇʳ = length(first(BRS)), length(BRS)
    M = Matrix{Int}(undef, Nⁱʳʳ - !includedim, Nᵉᵇʳ)
    @inbounds for (j, BR) in enumerate(BRS)
        for (i, v) in enumerate(BR.irvec)
            M[i,j] = v
        end
        if includedim
            M[Nⁱʳʳ,j] = dim(BR)
        end
    end

    return M
end 
@deprecate matrix(BRS::BandRepSet, includedim::Bool) matrix(BRS::BandRepSet; includedim=includedim)