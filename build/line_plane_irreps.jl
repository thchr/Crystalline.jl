"""
    LinePlaneIrreps

Derive the little group irreps of the plane (2D) and line (1D) groups from those of a parent
space group, and write them to `data/irreps/lgs/{2,1}d/`, in the layout of the 3D data (see
`write_dsg_irreps.jl`). Works for single- and double-valued irreps alike (via `spinful`).

## Call tree

 write_lgirreps(Val(D); spinful)               write all D-dimensional irreps to disk
 └─ derive_lgirreps_from_spacegroup(sgnum, Val(D); spinful)  ↓
    │                                          all irreps of one group, by k-label
    ├─ PARENTS[D][sgnum]                       the parent space group, and the embedding
    ├─ embedding_signs                         (spinful) SU(2) sign of each `embed(g)`
    └─ derive_lgirreps_from_spacegroup(lg, parent, …)  ↓
       │                                       the irreps of one little group `lg`
       ├─ find_source                          a tabulated k-manifold of the parent
       │  └─ locate                            whether a k-manifold contains another
       └─ remap_to_kstar                       (if needed) move it to the k-manifold
                                               that contains `embed(k)`

## Parent space groups

A plane or line group acts as the layer or rod group with no out-of-plane or out-of-line
symmetry (see `SU2`), and so as a space group in which the out-of-plane or out-of-line
directions are left fixed: its *parent*. The parent is that of `PLANE2SPACE_NUMS` for plane
groups, and of `FRIEZE2SPACE_NUMS[LINE2FRIEZE_NUMS[n]]` for line groups.

An operation `g = {W|w}` of the D-dimensional plane/line group then corresponds to the
parent's operation `embed(g)`, whose rotation part is `W ⊕ I(3-D)` (block-diagonal) and
translation part is `w ⊕ 0(3-D)`, up to a permutation of the parent's basis vectors
(see `Parent`). Likewise, a k-vector `k` corresponds to `embed(k) = k ⊕ 0(3-D)`.

## The irreps at a k-manifold

The little groups and their k-labels are those already stored in
`data/irreps/lgs/{2,1}d/littlegroups_data.jld2`. For each, `find_source` picks a k-manifold
tabulated for the parent, with the same little group (under `embed`), that _contains_
`embed(k)`: e.g., a line in 2D may lie in a plane in 3D. Often, the tabulated k-manifold
contains `embed(k)` only up to a reciprocal lattice vector, or only `h·embed(k)` for some
operation `h` of the parent (i.e., another arm of its star); the tabulated irreps are then
first moved to a k-manifold that contains `embed(k)` itself by `remap_to_kstar`.

The parent's irreps are stored as matrices `𝐏` and translations `𝛕`, with `D(g) = 𝐏 exp(2πi
𝐤⋅𝛕)` at every `𝐤` of the k-manifold, and so in particular at every `𝐤` of a k-manifold it
contains. The irrep of `g` is therefore given by the `𝐏` and `𝛕` of `embed(g)`, with `𝛕`
restricted to the plane or line (`embed(k)` has no component outside it).

For double-valued irreps, `embed(g)` must also be given an SU(2) element, which is fixed by
that of `g` only up to a sign for two-fold rotations and mirrors. The sign is chosen by
`_lift_signs`, mapping unbarred operations to unbarred operations as often as possible (in
practice, always).

Each irrep takes its label from the parent irrep, with the k-label replaced by its own:
e.g., Σ₃ˢ of p1m1 (3) derives from F₃ˢ of Pm (6). Likewise its reality, which is checked by
the Herring criterion.
"""
module LinePlaneIrreps

using Crystalline, JLD2, LinearAlgebra, StaticArrays
using Crystalline: PLANE2SPACE_NUMS, LINE2FRIEZE_NUMS, FRIEZE2SPACE_NUMS, AbstractLGIrrep,
                   reduce_ops, _ishexagonal, _lift_signs, _to_value, _doubled_matrices,
                   _doubled_translations, can_intersect, _check_valid_sgnum_and_dim,
                   _collect_lgirsd_sorted, nfreeparams, free

export write_lgirreps, derive_lgirreps_from_spacegroup, PARENTS

# ---------------------------------------------------------------------------------------- #
# Writing to disk

"""
    write_lgirreps(::Val{D}; spinful::Val) --> String

Derive the irreps of all D-dimensional space groups (see `derive_lgirreps_from_spacegroup`)
and write them to `data/irreps/lgs/{D}d/irreps_data_spinful.jld2` (if `spinful = Val(true)`)
or `data/irreps/lgs/{D}d/irreps_data.jld2` (if `spinful = Val(false)`). Only the operations
of the ordinary little group are stored, as for the 3D data.
"""
function write_lgirreps(::Val{D}; spinful::Val) where D
    D in (1, 2) || throw(ArgumentError("D must be 1 or 2"))
    S = _to_value(spinful)
    outpath = joinpath(dirname(@__DIR__), "data", "irreps", "lgs", "$(D)d",
                       S ? "irreps_data_spinful.jld2" : "irreps_data.jld2")
    # Crystalline keeps the data files open for reading, which blocks overwriting them
    jldfile = (S ? Crystalline.DLGIRREPS_JLDFILES : Crystalline.LGIRREPS_JLDFILES)[D]
    isassigned(jldfile) && close(jldfile[]) # close the open file, so we can write a new one

    JLD2.jldopen(outpath, "w") do f
        for sgnum in 1:MAX_SGNUM[D]
            lgirsd = derive_lgirreps_from_spacegroup(sgnum, Val(D); spinful)
            lgirsv = _collect_lgirsd_sorted(lgirsd)
            matrices_list     = Vector{Vector{Matrix{ComplexF64}}}[]
            translations_list = Vector{Union{Nothing, Vector{Vector{Float64}}}}[]
            for lgirs in lgirsv
                n = S ? order(group(lgirs)) ÷ 2 : order(group(lgirs)) # unbarred operations
                matrices = [lgir.matrices[1:n] for lgir in lgirs]
                translations = map(lgirs) do lgir
                    τs = lgir.translations[1:n]
                    all(iszero, τs) ? nothing : Vector{Float64}.(τs)
                end
                push!(matrices_list, matrices)
                push!(translations_list, translations)
            end
            f["$(sgnum)/klab_list"]         = klabel.(lgirsv)
            f["$(sgnum)/matrices_list"]     = matrices_list
            f["$(sgnum)/translations_list"] = translations_list
            f["$(sgnum)/realities_list"]    = [Int8.(reality.(lgirs)) for lgirs in lgirsv]
            f["$(sgnum)/cdml_list"]         = [label.(lgirs) for lgirs in lgirsv]
        end
    end
    println("wrote $(outpath)  ($(round(filesize(outpath)/2^10; digits=1)) KiB)")
    return outpath
end

# ---------------------------------------------------------------------------------------- #
# Parent space groups

"""
    Parent{D}

The parent space group `num` of a D-dimensional space group, together with the embedding of
the latter into the former: the `i`th basis vector of the D-dimensional plane/line group is
the `axes[i]`th conventional basis vector of the parent space group. The remaining basis
vector(s) of the parent span the out-of-plane (or out-of-line) directions, which its
operations leave fixed.

## Example
The parent of plane group p1m1 (3) is Pm (6); the mirror plane of Pm is normal to `𝐛`, as
`Parent{2}(6, (2, 3))`: the basis vectors `𝐚` and `𝐛` of p1m1 are `𝐛` and `𝐜` of Pm,
and `𝐚` of Pm is out-of-plane. Accordingly, the mirror `{-x,y}` of p1m1 is embedded as the
mirror `{x,-y,z}` of Pm, and the k-point `(k₁, k₂)` as `(0, k₁, k₂)`.
"""
struct Parent{D}
    num  :: Int
    axes :: NTuple{D, Int}
end

# A plane group lies in the 𝐚𝐛-plane of its parent (`axes = (1, 2)`), except if the parent
# is monoclinic, with unique axis 𝐛. `MONOCLINIC_PLANE_AXES` lists these exceptions:
# keys are plane group numbers, and values are the `axes` of their `Parent`, i.e., which of
# the parent's basis vectors 𝐚, 𝐛, 𝐜 (as 1, 2, 3) the plane group's 𝐚 and 𝐛 map to
const MONOCLINIC_PLANE_AXES = Dict( # plane sgnum => parent.axes
    2 => (1, 3),  # p2   ⇒ P2: in the 𝐚𝐜-plane, normal to the two-fold axis 𝐛
    3 => (2, 3),  # p1m1 ⇒ Pm: in the 𝐛𝐜-plane, normal to the mirror plane
    4 => (2, 3),  # p1g1 ⇒ Pc: in the 𝐛𝐜-plane, normal to the glide plane
    5 => (2, 1)   # c1m1 ⇒ Cm: in the 𝐛𝐚-plane, whose centring is that of Cm
)

"""
    PARENTS :: Tuple{Vector{Parent{1}}, Vector{Parent{2}}}

`PARENTS[D][n]` is the `Parent` of the D-dimensional space group `n`, i.e., of line group
`n` (`D = 1`) or plane group `n` (`D = 2`). A line group lies along the `𝐛` axis of its
parent, normal to the mirror plane of Pm.
"""
const PARENTS = (
    # 1D: always the 𝐛-axis of the parent
    [Parent{1}(FRIEZE2SPACE_NUMS[LINE2FRIEZE_NUMS[n]], (2,)) for n in 1:MAX_SGNUM[1]],
    # 2D: usually the 𝐚𝐛-plane of the parent, except for monoclinic parents
    [Parent{2}(PLANE2SPACE_NUMS[n], get(MONOCLINIC_PLANE_AXES, n, (1, 2)))
     for n in 1:MAX_SGNUM[2]]
)

# `embed(g)` and `embed(k)`, i.e., `W ⊕ I(3-D)`, `w ⊕ 0(3-D)`, and `k ⊕ 0(3-D)`, with the
# D-dimensional components placed at `parent.axes` (the same placement applies to direct and
# reciprocal coordinates, since the parent's basis is merely permuted)
function embed(op::SymOperation{D}, parent::Parent{D}) where D
    W, w = rotation(op), translation(op)
    W³, w³ = MMatrix{3,3,Float64}(I), zero(MVector{3,Float64})
    for i in 1:D
        w³[parent.axes[i]] = w[i]
        for j in 1:D
            W³[parent.axes[i], parent.axes[j]] = W[i,j]
        end
    end
    return SymOperation{3}(W³, w³)
end
function embed(kv::KVec{D}, parent::Parent{D}) where D
    k₀, kabc = parts(kv)
    k₀³, kabc³ = zero(MVector{3,Float64}), zero(MMatrix{3,3,Float64})
    for i in 1:D
        k₀³[parent.axes[i]] = k₀[i]
        for j in 1:D # the free parameters `αβγ` keep their order
            kabc³[parent.axes[i], j] = kabc[i,j]
        end
    end
    return KVec{3}(k₀³, kabc³)
end

# ---------------------------------------------------------------------------------------- #
# Deriving the irreps

"""
    derive_lgirreps_from_spacegroup(sgnum, ::Val{D}; spinful::Val)
                    --> Dict{String, Collection{LGIrrep{D}}} or
                        Dict{String, Collection{DLGIrrep{D}}}

Derive the single- or double-valued little group irreps of the D-dimensional plane or line
group with conventional number `sgnum` from those of its parent space group (see `PARENTS`),
for each little group of `littlegroups(sgnum, Val(D))`.
"""
function derive_lgirreps_from_spacegroup(sgnum::Integer, ::Val{D}; spinful::Val) where D
    D in (1, 2) || throw(ArgumentError("D must be 1 or 2"))
    _check_valid_sgnum_and_dim(sgnum, D)
    parent = PARENTS[D][sgnum]::Parent{D}
    lgirsd³ = lgirreps(parent.num, Val(3); spinful)
    lgs³ = littlegroups(parent.num, Val(3)) # spatial operations only (no SU2), to match k-manifolds
    s₀ = _to_value(spinful) ? embedding_signs(sgnum, Val(D), parent) : nothing
    lgs = littlegroups(sgnum, Val(D); spinful)
    return Dict(klab => derive_lgirreps_from_spacegroup(lg, parent, lgs³, lgirsd³, s₀) for (klab, lg) in lgs)
end

"""
    derive_lgirreps_from_spacegroup(lg, parent, lgs³, lgirsd³, s₀) --> Collection

Derive the irreps of the little group `lg` from the parent's irreps `lgirsd³` (with little
groups `lgs³`). `s₀` holds the `signbit`s of `embedding_signs` if `lg` is a double group,
and is `nothing` otherwise.
"""
function derive_lgirreps_from_spacegroup(
    lg::AbstractLittleGroup{D},
    parent::Parent{D},
    lgs³::AbstractDict{String, <:AbstractLittleGroup{3}},
    lgirsd³::AbstractDict{String, <:Collection{<:AbstractLGIrrep{3}}},
    s₀::Union{Nothing, AbstractDict},
) where D
    spinful = isspinful(lg)
    kv = position(lg)
    hexagonal³ = _ishexagonal(parent.num, Val(3))

    # find parent's irreps, on a k-manifold that contains `embed(kv)` (via `find_source`)
    lg′ = if spinful # `lg′` contains *only* the *spatial* (not SU2) operations of `lg`
        ops′ = [dop.op for dop in operations(lg)[1:order(lg)÷2]]
        LittleGroup{D}(num(lg), kv, klabel(lg), ops′)
    else
        lg
    end
    klab³, kv′, h = find_source(lg′, parent, lgs³)
    lgirs³ = lgirsd³[klab³]
    if !isapprox(kv′, position(first(lgirs³)), nothing, #=modw=# false)
        # target `kv′` is only *equivalent* to the k-manifold of `position(lgirs³)` but not
        # equal; and is instead merely related to `position(lgirs³)` by the operation `h⁻¹`
        # (possibly modulo a reciprocal lattice vector): remap `lgirs³` to `kv′`. We pass
        # `h⁻¹` as the only coset representative, since `find_source` constructed `kv′` as
        # the image of `position(lgirs³)` under `h⁻¹`: the search for a representative
        # that maps one onto the other then cannot fail, nor pick a different one
        h⁻¹ = spinful ? DSymOperation(inv(h), hexagonal³) : inv(h)
        lgirs³ = remap_to_kstar(lgirs³, kv′, [h⁻¹])
    end

    # the index of each `embed(g)` among the parent's operations, with the SU(2) element
    # chosen by `embedding_signs` if spinful
    ops³ = operations(first(lgirs³))
    idxs = map(lg′) do op
        op³ = embed(op, parent)
        if spinful
            u = SU2(op³, hexagonal³)
            op³ = DSymOperation{3}(op³, s₀[rotation(op)] ? -u : u)
        end
        @something(
            findfirst(op′ -> isapprox(op′, op³, nothing, #=modw=# false), ops³),
            error("⋕$(num(lg)), $(klabel(lg)): $(seitz(op³)) is not an operation of the \
                   parent's little group at $(klab³)")
        )
    end

    # the space group without centring copies, as needed by `calc_reality`
    sgops = reduce_ops(spacegroup(num(lg), Val(D); spinful=Val(spinful)), true)

    # construct the irreps of `lg` from its parents'
    lgirs = map(lgirs³) do lgir³
        Ps = lgir³.matrices[idxs]
        τs = [[τ[a] for a in parent.axes] for τ in lgir³.translations[idxs]]

        cdml = klabel(lg) * chopprefix(label(lgir³), klabel(lgir³))
        lgir = if spinful
            # `Ps` and `τs` only cover the unbarred half of operations of `lg` (namely those
            # in `lg′`): so start by getting the missing half
            Ps = _doubled_matrices(Ps)
            τs = _doubled_translations(τs)
            DLGIrrep{D}(cdml, lg, Ps, τs, reality(lgir³))
        else
            LGIrrep{D}(cdml, lg, Ps, τs, reality(lgir³))
        end
        if calc_reality(lgir, sgops) ≠ reality(lgir) # sanity-check
            error("reality of $cdml (⋕$(num(lg))) differs from that of $(label(lgir³))")
        end
        lgir
    end
    return Collection(lgirs)
end

"""
    find_source(lg, parent, lgs³) --> (klab³, kv′, h)

Find a k-manifold `klab³` among the parent's little groups `lgs³` whose irreps are those of
`lg`. Returns, together with `klab³`, a k-manifold `kv′` that contains `embed(position(lg))`
exactly, and an operation `h` of the parent, such that `kv′` is equivalent to the k-manifold
of `klab³` by `h⁻¹`, up to a reciprocal lattice vector (i.e., such that the irreps of
`klab³` can be moved to `kv′` by `remap_to_kstar`).

Prefers `h` to be the identity, and then the lowest-dimensional k-manifold.
"""
function find_source(lg::LittleGroup{D}, parent::Parent{D}, lgs³::AbstractDict) where D
    cntr³ = centering(parent.num, 3)
    kv_embed = embed(position(lg), parent)
    ops_embed = [embed(op, parent) for op in lg]
    best = nothing
    sgops³ = reduce_ops(spacegroup(parent.num, Val(3)), true)
    for h in sgops³
        for (klab³, lg³) in lgs³
            kv³ = position(lg³) # a tabulated k-manifold, to compare with `h·kv_embed`
            # the little groups must have the same order …
            length(lg³) == length(lg) || continue
            # … and `kv³` must contain `h·kv_embed`, up to a reciprocal lattice vector `G`
            # (otherwise, its irreps are irreps at other k-points) …
            containment = locate(h * kv_embed, kv³, cntr³)
            containment.bool || continue
            # … and the little group of `kv³` must be that of `kv_embed`, conjugated by
            # `h`: a k-manifold that contains another can have a smaller little group, even
            # if, as here, both groups have the same order
            all(ops_embed) do op
                op′ = compose(h, compose(op, inv(h), false), false)
                any(op³ -> isapprox(op³, op′, cntr³, #=modw=# true), lg³)
            end || continue

            # score the match of `h` and `lg³`: prefer `h` to be the identity; and otherwise
            # prefer the lowest-dimensional k-manifold
            score = (!isone(h), nfreeparams(kv³))
            if isnothing(best) || score < first(best)
                k₀³, kabc³ = parts(kv³)
                kv′ = inv(h) * KVec{3}(k₀³ + containment.G, kabc³)
                best = (score, (klab³, kv′, h))
            end
        end
    end
    isnothing(best) && error("no k-manifold of ⋕$(parent.num) matches $(klabel(lg))")
    return last(best)
end

"""
    locate(kv, kv′, cntr) --> (; bool, G)

Return whether the k-manifold `kv` is contained in the k-manifold `kv′`, modulo a reciprocal
lattice vector `G` of the `cntr`-centred lattice, and if so (`bool = true`), that vector, in
the conventional basis: i.e., such that every point of `kv` is a point of `kv′`, plus `G`.
"""
function locate(kv::KVec{3}, kv′::KVec{3}, cntr::Char)
    # `kv` is contained in `kv′` iff they intersect (checked by `can_intersect` below) *and*
    # if every free direction of `kv` is also a free direction of `kv′`: i.e., the
    # directions spanned by the columns of `free(kv)` must lie in the span of those of
    # `free(kv′)`, so that appending the former to the latter does not increase its rank.
    # E.g., the lines `kv = [α, 0, 0]` and `kv′ = [0, α, 0]` intersect (at `Γ`), but neither
    # contains the other, and `rank(hcat(free(kv′), free(kv))) = 2 ≠ rank(free(kv′)) = 1`.
    if rank(hcat(free(kv′), free(kv))) ≠ rank(free(kv′))
        return (; bool=false, G=zero(SVector{3,Int})) # not contained
    end
    # compare in a primitive basis, so that only the reciprocal lattice vectors of the
    # centred lattice count (see `can_intersect`)
    intersection = can_intersect(primitivize(kv, cntr), primitivize(kv′, cntr))
    intersection.bool || return (; bool=false, G=zero(SVector{3,Int})) # not contained

    # `G` is `intersection.L`, but in the conventional rather than the primitive basis
    G = kv(intersection.αβγ) - kv′(intersection.αβγ′)
    all(x -> abs(x - round(x)) < 1e-10, G) || error("unexpected non-integer G = $G")
    return (; bool=true, G=round.(Int, G)) # contained (mod G)
end

"""
    embedding_signs(sgnum, ::Val{D}, parent) --> Dict

For each rotation part `W` of the `D`-dimensional space group `sgnum`, the signbit `s`
(`true`: negative sign; `false`: positive sign) with which the (unbarred) SU(2) element of
an operation `g` with rotation part `W` corresponds to that of `embed(g)` in the parent
(see `_lift_signs`). The signs do not depend on the k-point.
"""
function embedding_signs(sgnum::Integer, ::Val{D}, parent::Parent{D}) where D
    pgops = pointgroup(spacegroup(sgnum, Val(D)))
    hexagonal, hexagonal³ = _ishexagonal(sgnum, Val(D)), _ishexagonal(parent.num, Val(3))
    dops = [DSymOperation(op, hexagonal) for op in pgops]
    dops³ = [DSymOperation(embed(op, parent), hexagonal³) for op in pgops]
    signs = _lift_signs(dops, dops³)
    return Dict(rotation(op) => signbit(s) for (op, s) in zip(pgops, signs))
end

end # module LinePlaneIrreps
