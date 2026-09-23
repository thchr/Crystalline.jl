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
struct DSpaceGroup{D} <: AbstractSpaceGroup{D, DSymOperation{D}}
    num        :: Int
    operations :: Vector{DSymOperation{D}}
end
label(sg::DSpaceGroup) = iuc(num(sg), dim(sg))

"""
$(TYPEDEF)$(TYPEDFIELDS)
"""
struct DPointGroup{D} <: AbstractPointGroup{D, DSymOperation{D}}
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
struct DLittleGroup{D} <: AbstractLittleGroup{D, DSymOperation{D}}
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
struct DSiteGroup{D} <: AbstractSiteGroup{D, DSymOperation{D}}
    num        :: Int
    wp         :: WyckoffPosition{D}
    operations :: Vector{DSymOperation{D}}
    cosets     :: Vector{DSymOperation{D}}
end
Base.position(g::DSiteGroup) = g.wp
label(g::DSiteGroup) = iuc(num(g), dim(g))
centering(::DSiteGroup) = nothing
cosets(g::DSiteGroup) = g.cosets

positionlabel(g::DLittleGroup) = klabel(g)
positionlabel(g::DSiteGroup) = label(position(g))

# --- construction ---
"""
    doubled_operations(g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}, SiteGroup{3}})
                                                            --> Vector{DSymOperation{3}}

Attach the SU(2) element (see [`su2`](@ref)) to each operation of `g`, returning the `2|G|`
operations of the associated double group: the operations themselves first, then their
``\\bar{E}``-barred partners, in the same order.
"""
function doubled_operations(
    g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}, SiteGroup{3}}
)
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

"""
    doublegroup(g::Union{SpaceGroup{3}, LittleGroup{3}, PointGroup{3}, SiteGroup{3}})
                    --> DSpaceGroup{3}, DLittleGroup{3}, DPointGroup{3}, or DSiteGroup{3}

Return the double group of the space, little, point, or site symmetry group `g` (see
[`doubled_operations`](@ref)).

The coset representatives of a site symmetry group are lifted to the double group with
their unbarred SU(2) elements.
"""
doublegroup(sg::SpaceGroup{3}) = DSpaceGroup{3}(num(sg), doubled_operations(sg))
function doublegroup(lg::LittleGroup{3})
    return DLittleGroup{3}(num(lg), position(lg), klabel(lg), doubled_operations(lg))
end
doublegroup(pg::PointGroup{3}) = DPointGroup{3}(num(pg), label(pg), doubled_operations(pg))
function doublegroup(siteg::SiteGroup{3})
    hexagonal = _ishexagonal(siteg)
    cosets′ = [DSymOperation{3}(op, su2(op, hexagonal)) for op in cosets(siteg)]
    return DSiteGroup{3}(num(siteg), position(siteg), doubled_operations(siteg), cosets′)
end

# --- matching two double groups ---
"""
    _lift_signs(ops::AbstractVector{DSymOperation{3}},
                ops′::AbstractVector{DSymOperation{3}})  -->  Vector{Int}

Given the unbarred operations `ops` and `ops′` of two double groups whose spatial parts
correspond one-to-one (`ops[i] ↔ ops′[i]`, as a group isomorphism), return signs `s` such
that mapping `ops[i]` to `ops′[i]` if `s[i] = 1`, or to its barred partner if `s[i] = -1`,
is an isomorphism of the double groups.

Matching the spatial parts does not settle the signs: for two-fold rotations and mirrors,
which of the two SU(2) elements is unbarred is a convention (see [`isbarred`](@ref)), and a
change of setting need not respect it. Several choices of signs are valid in general; they
differ by a one-dimensional ±1 representation of the group. We return the valid choice with
the fewest `-1`s, i.e., the one that maps unbarred operations to unbarred operations as
often as possible (on ties, the first found).

For all site symmetry groups of the 230 space groups, the tied choices give the same
irrep characters.
"""
function _lift_signs(
    ops::AbstractVector{DSymOperation{3}},
    ops′::AbstractVector{DSymOperation{3}}
)
    n = length(ops)
    # products: `ops[i]*ops[j]` is `c[i,j]` times `ops[k[i,j]]` (unbarred), and likewise
    # `c′` for `ops′`; the spatial parts compose alike, so `k` is shared
    k, c = _unbarred_products(ops)
    k′, c′ = _unbarred_products(ops′)
    k == k′ || error("the spatial parts of the two groups do not correspond")

    # a valid `s` must satisfy `s[i] s[j] c′[i,j] = c[i,j] s[k[i,j]]`; so `s` is fixed by
    # its values on a set of generators. NB: index 1 must be the identity (as it is in every
    # group built by Crystalline), both here and in `_generated`
    gens = Int[]
    for i in 1:n
        i ∈ _generated(gens, k) || push!(gens, i)
    end
    best = Int[]
    for bits in 0:(2^length(gens) - 1)
        s = zeros(Int, n)
        s[1] = 1
        for (b, g) in enumerate(gens)
            s[g] = isodd(bits >> (b-1)) ? -1 : 1
        end
        todo = [1; gens]
        while !isempty(todo)
            i = pop!(todo)
            for g in gens
                j = k[i, g]
                if iszero(s[j])
                    s[j] = s[i] * s[g] * c′[i, g] * c[i, g]
                    push!(todo, j)
                end
            end
        end
        all(s[k[i,j]] == s[i] * s[j] * c′[i,j] * c[i,j] for i in 1:n, j in 1:n) || continue
        if isempty(best) || count(==(-1), s) < count(==(-1), best)
            best = s
        end
    end
    isempty(best) && error("the two double groups are not isomorphic under the given \
                            correspondence")
    return best
end

function _unbarred_products(ops::AbstractVector{DSymOperation{3}})
    n = length(ops)
    k = Matrix{Int}(undef, n, n)
    c = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        W = rotation(ops[i]) * rotation(ops[j])
        k[i,j] = something(findfirst(op -> rotation(op) ≈ W, ops))
        u = ops[i].su2 * ops[j].su2
        c[i,j] = u ≈ ops[k[i,j]].su2 ? 1 : -1
    end
    return k, c
end

function _generated(gens::AbstractVector{Int}, k::AbstractMatrix{Int})
    S = Set([1])
    todo = [1]
    while !isempty(todo)
        i = pop!(todo)
        for g in gens
            j = k[i, g]
            j ∈ S || (push!(S, j); push!(todo, j))
        end
    end
    return S
end
