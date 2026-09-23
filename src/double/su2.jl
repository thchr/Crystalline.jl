# --- SU(2) elements of spatial operations ---

"""
    su2(op::SymOperation{3}, sgnum::Integer)   --> SU2
    su2(op::SymOperation{3}, hexagonal::Bool)  --> SU2

Return the SU(2) element of the spatial operation `op`: the 2×2 matrix by which `op` acts
on spin-½ degrees of freedom.

The crystal system must be supplied, either as a space group number or directly as whether
it is hexagonal or trigonal. It is needed because the SU(2) element is not fixed by the
rotation part alone: a rotation matrix in fractional coordinates does not determine the
Cartesian rotation axis, and four rotation parts are assigned differently in the hexagonal
and trigonal settings than elsewhere.

The assignment follows Altmann & Herzig, *Point-Group Theory Tables* (1994), as used by the
Bilbao Crystallographic Server. For a rotation by `φ` about the Cartesian unit axis `𝐧`,

``U = \\exp[-i(φ/2)\\,𝐧⋅𝛔] = \\cos(φ/2)𝟙 - i\\sin(φ/2)\\,(𝐧⋅𝛔)``

with ``𝛔 = (σ_x, σ_y, σ_z)`` the Pauli matrices.

An improper operation acts through its proper part alone, spin being axial; inversion
therefore maps to the identity.

The values are tabulated rather than evaluated from that expression: the expression fixes
`U` for every operation with ``φ ≠ π``, but two-fold rotations and mirrors have ``φ = π``,
where `U` and `-U` describe the same spatial operation and the choice between them is
convention.

## Cartesian frame
`op` must be given in a conventional setting, whose basis is taken in a fixed orientation
relative to a Cartesian frame: this fixes the physical operation, and hence its SU(2)
element. The tabulated elements follow the Bilbao Crystallographic Server and Altmann &
Herzig, whose frame, relative to the conventional basis `(𝐚, 𝐛, 𝐜)`, is:
- **all but hexagonal and trigonal lattices**: `𝐚 ∥ x`, `𝐛` in the `xy`-plane; the same as
  Crystalline's (e.g., as used by [`directbasis`](@ref) and [`crystal`](@ref));
- **hexagonal and trigonal lattices**: `x ∥ 𝐚 + 2𝐛`, `y ∥ 𝐚`, and `z ∥ -𝐜`. This follows
  from Altmann & Herzig's rules placing the two classes of binary axes perpendicular to the
  principal axis along `x` and `y`. It differs from Crystalline's frame (`𝐚 ∥ x`,
  `𝐜 ∥ z`) by a two-fold rotation about `(1, 1, 0)` (or, equivalently for every rotation
  of these lattices, about `(1, -1, 0)`).
We keep the tabulated elements unchanged, so that they agree with those of the Bilbao
Crystallographic Server and of Altmann & Herzig. The frame does not affect any
group-theoretical result (irreps, characters, band representations); it matters only when
SU(2) elements are combined with quantities that refer to Cartesian axes (e.g., spin
operators taken as physical directions, orbitals that transform under Cartesian rotations,
or external fields), which must then use the same frame.

## Changes of basis and of setting
The element stays with the physical operation under [`transform`](@ref) (and so
[`primitivize`](@ref) and [`conventionalize`](@ref)): the operation is rewritten in another
basis, but not changed. When `transform` is instead used for a change of *setting* (e.g.,
with [`conjugacy_relations`](@ref)), the new basis is generally not in the orientation
above, and the SU(2) elements need not equal those tabulated for the new setting (e.g.,
`2₀₀₁` written as `2₀₁₀` keeps `-iσz`, where the table has `-iσy` for `2₀₁₀`). To obtain the
double group operations of the new setting, transform the spatial operations and attach
the SU(2) elements tabulated for that setting, i.e., `su2(op′, sgnum′)`.
"""
su2(op::SymOperation{3}, sgnum::Integer) = su2(op, _ishexagonal(sgnum))
# whether space group `sgnum`, or group `g`, uses the hexagonal Cartesian frame (see
# `SU2_BY_ROTATION_HEX`); point groups 16-27 are the trigonal and hexagonal ones
_ishexagonal(sgnum::Integer) = crystalsystem(sgnum, 3) ∈ ("hexagonal", "trigonal")
_ishexagonal(g::Union{SpaceGroup{3}, LittleGroup{3}, SiteGroup{3}}) = _ishexagonal(num(g))
_ishexagonal(pg::PointGroup{3}) = 16 ≤ num(pg) ≤ 27
function su2(op::SymOperation{3}, hexagonal::Bool)
    k = _rotation_key(op)
    if hexagonal
        u = get(SU2_BY_ROTATION_HEX, k, nothing)
        u === nothing || return u
    end
    u = get(SU2_BY_ROTATION, k, nothing)
    u === nothing && throw(DomainError(rotation(op),
        "no tabulated SU(2) element for this rotation part; the operation must be given " *
        "in a conventional setting of a crystallographic group"))
    return u
end

function _rotation_key(op::SymOperation{3})
    W = rotation(op)
    all(w -> abs(w - round(w)) ≤ DEFAULT_ATOL, W) || throw(DomainError(W,
        "rotation part is not integer-valued; SU(2) elements are tabulated only for " *
        "operations in a conventional (fractional-coordinate) setting"))
    return NTuple{9,Int}(round.(Int, W))
end

"""
    _binary_axis(u::SU2) --> SVector{3,Float64}

Return the oriented rotation axis `𝐧` of a two-fold rotation or mirror whose SU(2) element
is `u = -i(𝐧⋅𝛔)`. Since `-𝐧` gives `-u`, the orientation of the axis picks out one of the
two SU(2) elements of the operation. Altmann calls `𝐧` the operation's *pole*.
"""
_binary_axis(u::SU2) = SVector{3,Float64}(-imag(u.b), -real(u.b), -imag(u.a))

# The two-fold axis directions that occur in the conventional settings, with each lattice in
# the Cartesian frame of the tables above, and oriented as Altmann orients them. Nothing is
# special about the directions themselves; the orientations are convention, and are tabulated
# rather than derived because Altmann fixes them as a set (by requiring the matrices to
# represent the double group of `D₂`) and not by any rule on the individual axis.
const SU2_BINARY_AXES = let ns = SVector{3,Float64}[]
    for tbl in (SU2_BY_ROTATION, SU2_BY_ROTATION_HEX), u in values(tbl)
        if abs(real(u.a)) ≤ DEFAULT_ATOL
            n = _binary_axis(u)
            any(m -> abs(dot(n, m)) > 1 - 1e-8, ns) || push!(ns, n)
        end
    end
    ns
end

"""
    isbarred(u::SU2)            --> Bool
    isbarred(op::DSymOperation) --> Bool

Return whether `u` is the barred one of the two SU(2) elements that share a spatial
operation.

Every spatial operation has two SU(2) elements, `U` and `-U`, since a 2π rotation acts as
`-𝟙` on spin; the barred element is the one obtained by composing with that 2π rotation.
Which of the two is *called* barred is a convention. We follow Altmann's, in which the
rotation angle is taken in ``[-π, π]``, so an unbarred element always has
`real(u.a) = cos(φ/2) ≥ 0`. The sign of `real(u.a)` therefore settles every operation with
``φ ≠ π``.

Two-fold rotations and mirrors have `real(u.a) = 0` for both elements, and are settled by
their oriented "axis" (see `_binary_axis`) instead, against `SU2_BINARY_AXES`: the
directions that two-fold axes take in the conventional settings of the space groups, with
each lattice in the Cartesian frame used by the Bilbao tables, and oriented as Altmann
orients them. Any other direction lies outside Altmann's tables; there, we take as unbarred
the element whose axis has a positive final non-zero component.

Reads only the SU(2) element, never the rotation part, and is hence unaffected by a change
of lattice basis (e.g., to a primitive setting).
"""
function isbarred(u::SU2)
    r = real(u.a)
    abs(r) > DEFAULT_ATOL && return r < 0
    n = _binary_axis(u)
    for m in SU2_BINARY_AXES                    # unit vectors, along distinct axes
        d = dot(n, m)
        d >  1 - 1e-8 && return false
        d < -1 + 1e-8 && return true
    end
    return n[findlast(x -> abs(x) > DEFAULT_ATOL, n)] < 0   # an untabulated direction
end
isbarred(op::DSymOperation) = isbarred(op.su2)
