# --- Double group types ---
# Each mirrors its spinless counterpart field-for-field, differing only in operation type.
# They are separate types, rather than a type parameter on the spinless ones, so that
# `SymOperation`, `LittleGroup`, `LGIrrep` and the rest keep their current definitions and
# previously serialized data still loads.
#
# A double group stores all `2|G|` operations, the barred ones included, so that it is a
# group in its own right: multiplication tables, conjugacy classes and the Herring
# criterion then need no special casing.

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct DSpaceGroup{D} <: AbstractGroup{D, DSymOperation{D}}
    num        :: Int
    operations :: Vector{DSymOperation{D}}
end
label(sg::DSpaceGroup) = iuc(num(sg), dim(sg))

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct DPointGroup{D} <: AbstractGroup{D, DSymOperation{D}}
    num        :: Int
    label      :: String
    operations :: Vector{DSymOperation{D}}
end
label(pg::DPointGroup) = pg.label
iuc(pg::DPointGroup) = label(pg)
centering(::DPointGroup) = nothing

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct DLittleGroup{D} <: AbstractGroup{D, DSymOperation{D}}
    num        :: Int
    kv         :: KVec{D}
    klab       :: String
    operations :: Vector{DSymOperation{D}}
end
Base.position(lg::DLittleGroup) = lg.kv
klabel(lg::DLittleGroup) = lg.klab
label(lg::DLittleGroup) = iuc(num(lg), dim(lg))

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct DSiteGroup{D} <: AbstractGroup{D, DSymOperation{D}}
    num        :: Int
    wp         :: WyckoffPosition{D}
    operations :: Vector{DSymOperation{D}}
    cosets     :: Vector{SymOperation{D}}
end
Base.position(g::DSiteGroup) = g.wp
label(g::DSiteGroup) = iuc(num(g), dim(g))
centering(::DSiteGroup) = nothing

# --- construction ---
"""
    doubled_operations(ops, hexagonal::Bool) --> Vector{DSymOperation{3}}

Attach the SU(2) element (see [`su2`](@ref)) to each operation of `ops`, returning the
`2|G|` operations of the associated double group: the operations themselves first, then
their ``\\bar{E}``-barred partners, in the same order.
"""
function doubled_operations(ops::AbstractVector{SymOperation{3}}, hexagonal::Bool)
    n = length(ops)
    dops = Vector{DSymOperation{3}}(undef, 2n)
    for (i, op) in enumerate(ops)
        u = su2(op, hexagonal)
        dops[i]   = DSymOperation{3}(op,  u)
        dops[i+n] = DSymOperation{3}(op, -u)
    end
    return dops
end

@noinline _only_3d(D) =
    throw(DomainError(D, "double groups are currently only supported in 3D"))

"""
    spacegroup(sgnum::Integer, Dᵛ::Val{D}, spinfulᵛ::Val{S})
                                        --> SpaceGroup{D} or DSpaceGroup{D}

Return the space group `sgnum`, as its double (spinful) group if `spinfulᵛ` is `Val(true)`.

Double groups are currently supported in 3D only.
"""
function spacegroup(sgnum::Integer, Dᵛ::Val{D}, ::Val{true}) where D
    D == 3 || _only_3d(D)
    sg = spacegroup(sgnum, Dᵛ)
    hex = crystalsystem(sgnum, 3) ∈ ("hexagonal", "trigonal")
    return DSpaceGroup{D}(sgnum, doubled_operations(operations(sg), hex))
end
spacegroup(sgnum::Integer, Dᵛ::Val, ::Val{false}) = spacegroup(sgnum, Dᵛ)
function spacegroup(sgnum::Integer, D::Integer, spinful::Bool) # type-unstable convenience
    return spacegroup(sgnum, Val(D), Val(spinful))
end
