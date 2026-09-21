# --- SU2 ---
"""
$(TYPEDEF)$(TYPEDFIELDS)

An element of SU(2), stored by the two complex parameters `a` and `b` of

``U = \\begin{pmatrix} a & b \\\\ -b^* & a^* \\end{pmatrix}``

with ``|a|^2 + |b|^2 = 1``. This is the spin-½ part of a double group operation: the two
SU(2) elements of a spatial operation differ by an overall sign, `u` and `-u`, and that
sign is what distinguishes an operation from its ``\\bar{E}``-barred partner.
"""
@struct_hash_equal struct SU2 <: AbstractMatrix{ComplexF64}
    a :: ComplexF64
    b :: ComplexF64
    @inline function SU2(a::ComplexF64, b::ComplexF64)
        # `det(U) = |a|² + |b|²`, so this is the full SU(2) constraint for the stored form
        @boundscheck if !isapprox(abs2(a) + abs2(b), 1, atol=DEFAULT_ATOL)
            throw(DomainError((a, b),
                "SU(2) elements must be normalized, i.e. have det(U) = |a|²+|b|² = 1"))
        end
        return new(a, b)
    end
end
SU2(a::Number, b::Number) = SU2(ComplexF64(a), ComplexF64(b))
function SU2(U::AbstractMatrix{<:Number})
    @boundscheck begin
        size(U) == (2,2) || throw(DomainError(size(U), "matrix size must be (2,2)"))
        # only `U[1,1]` and `U[1,2]` are stored, so check that the rest is consistent
        (isapprox(U[2,1], -conj(U[1,2]); atol=DEFAULT_ATOL) &&
         isapprox(U[2,2],  conj(U[1,1]); atol=DEFAULT_ATOL)) || throw(DomainError(U,
            "matrix must be of the SU(2) form [a b; -b* a*]"))
    end
    return SU2(U[1,1], U[1,2])
end

matrix(u::SU2) = SMatrix{2,2,ComplexF64}(u.a, -conj(u.b), u.b, conj(u.a))

# ::: AbstractArray interface :::
Base.size(::SU2) = (2, 2)
Base.IndexStyle(::Type{SU2}) = IndexCartesian()
@propagate_inbounds Base.getindex(u::SU2, i::Int, j::Int) = matrix(u)[i, j]

# `(a, b)` multiply as the corresponding matrices do; written out to avoid building them.
# A product of SU(2) elements is normalized already, so skip the constructor's check.
(*)(u₁::SU2, u₂::SU2) =
    @inbounds SU2(u₁.a*u₂.a - u₁.b*conj(u₂.b), u₁.a*u₂.b + u₁.b*conj(u₂.a))
one(::Type{SU2}) = SU2(one(ComplexF64), zero(ComplexF64))
one(::SU2) = one(SU2)
inv(u::SU2) = SU2(conj(u.a), -u.b)
(-)(u::SU2) = SU2(-u.a, -u.b)
function Base.isapprox(u₁::SU2, u₂::SU2; atol::Real=DEFAULT_ATOL, kws...)
    return isapprox(u₁.a, u₂.a; atol, kws...) && isapprox(u₁.b, u₂.b; atol, kws...)
end

# --- DSymOperation ---
"""
$(TYPEDEF)$(TYPEDFIELDS)

A double group operation: a spatial operation `op` together with its SU(2) element `su2`
(see [`su2`](@ref)), which is how it acts on spin-½ degrees of freedom.

The SU(2) element is carried rather than recomputed, because it is not determined by the
rotation part alone (see [`su2`](@ref)) — and because carrying it is what keeps
[`compose`](@ref) closed.

Whether the operation is ``\\bar{E}``-barred is *not* stored: it follows from `su2` alone,
in any setting, so [`isbarred`](@ref) computes it on demand rather than every
[`compose`](@ref) maintaining it.
"""
@struct_hash_equal struct DSymOperation{D} <: AbstractOperation{D}
    op  :: SymOperation{D}
    su2 :: SU2
end
SymOperation{D}(dop::DSymOperation{D}) where D = dop.op
SymOperation(dop::DSymOperation) = dop.op

"""
    su2(dop::DSymOperation) --> SU2

Return the SU(2) part of the double group operation `dop`.
"""
su2(dop::DSymOperation) = dop.su2

function compose(dop₁::DSymOperation{D}, dop₂::DSymOperation{D}, modτ::Bool=true) where D
    return DSymOperation{D}(compose(dop₁.op, dop₂.op, modτ), dop₁.su2 * dop₂.su2)
end
(*)(dop₁::DSymOperation{D}, dop₂::DSymOperation{D}) where D = compose(dop₁, dop₂)

inv(dop::DSymOperation{D}) where D = DSymOperation{D}(inv(dop.op), inv(dop.su2))

one(::Type{DSymOperation{D}}) where D = DSymOperation{D}(one(SymOperation{D}), one(SU2))
# a pure lattice translation, which acts trivially on spin
_translation_operation(::Type{DSymOperation{D}}, t) where D =
    DSymOperation{D}(SymOperation{D}(t), one(SU2))
one(dop::DSymOperation) = one(typeof(dop))
# the SU(2) parameters are irrational for most operations, so unlike `isone(::SymOperation)`
# the check below must be approximate
# Ē carries `su2 = -𝟙`, which is not `≈ one(SU2)`, so no barring check is needed here
isone(dop::DSymOperation) = isone(dop.op) && isapprox(dop.su2, one(SU2))

function Base.isapprox(
    dop₁::DSymOperation{D},
    dop₂::DSymOperation{D},
    vs...;
    kws...
) where D
    return isapprox(dop₁.su2, dop₂.su2) && isapprox(dop₁.op, dop₂.op, vs...; kws...)
end

xyzt(dop::DSymOperation) = xyzt(dop.op) # does not feature the SU(2) part!
function seitz(dop::DSymOperation)
    s = seitz(dop.op)
    isbarred(dop) || return s
    # `seitz` omits the braces when the translation part vanishes
    return startswith(s, '{') ? "{ᵈ" * SubString(s, 2) : "ᵈ" * s
end

# --- change of lattice basis ---
# A change of lattice basis keeps the Cartesian frame fixed, so the SU(2) element is
# unchanged. (A rotation of the Cartesian frame by `V` would instead act as `U → VUV†`.)
function transform(
    dop::DSymOperation{D},
    P::AbstractMatrix{<:Real},
    p::Union{AbstractVector{<:Real}, Nothing}=nothing,
    modw::Bool=true
) where D
    return DSymOperation{D}(transform(dop.op, P, p, modw), dop.su2)
end
function primitivize(dop::DSymOperation{D}, cntr::Char, modw::Bool=true) where D
    return DSymOperation{D}(primitivize(dop.op, cntr, modw), dop.su2)
end
function conventionalize(dop::DSymOperation{D}, cntr::Char, modw::Bool=true) where D
    return DSymOperation{D}(conventionalize(dop.op, cntr, modw), dop.su2)
end
