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
    doubled_operations(g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}})
                                                            --> Vector{DSymOperation{3}}

Attach the SU(2) element (see [`su2`](@ref)) to each operation of `g`, returning the `2|G|`
operations of the associated double group: the operations themselves first, then their
``\\bar{E}``-barred partners, in the same order.
"""
function doubled_operations(g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}})
    hexagonal = _ishexagonal(g)
    n = length(g)
    dops = Vector{DSymOperation{3}}(undef, 2n)
    for (i, op) in enumerate(g)
        u = su2(op, hexagonal)
        dops[i]   = DSymOperation{3}(op,  u)
        dops[i+n] = DSymOperation{3}(op, -u)
    end
    return dops
end

@noinline _only_3d(D) = throw(DomainError(D, "double groups are currently only supported in 3D"))

centering(g::Union{DSpaceGroup{D}, DLittleGroup{D}}) where D = centering(num(g), D)

"""
    doublegroup(g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}})
                                --> DSpaceGroup{3}, DLittleGroup{3}, or DPointGroup{3}

Return the double group of the space, little, or point group `g` (see
[`doubled_operations`](@ref)).
"""
doublegroup(sg::SpaceGroup{3}) = DSpaceGroup{3}(num(sg), doubled_operations(sg))
function doublegroup(lg::LittleGroup{3})
    return DLittleGroup{3}(num(lg), position(lg), klabel(lg), doubled_operations(lg))
end
doublegroup(pg::PointGroup{3}) = DPointGroup{3}(num(pg), label(pg), doubled_operations(pg))

# --- change of lattice basis ---
# Little groups hold no centring copies, so no operations become equivalent (unlike for
# space groups, cf. `reduce_ops`)
function primitivize(lg::DLittleGroup{D}, modw::Bool=true) where D
    cntr = centering(lg)
    kv′  = primitivize(position(lg), cntr)
    ops′ = primitivize.(operations(lg), cntr, modw)
    return DLittleGroup{D}(num(lg), kv′, klabel(lg), ops′)
end
