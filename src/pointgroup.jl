# ===== CONSTANTS =====

# We include several axes settings; as a result, there are more than 32 point groups 
# under the 3D case, because some variations are "setting-degenerate" (but are needed
# to properly match all space group settings)
const PG_NUM2IUC = (
    (["1"], ["m"]),                                           # 1D
    (["1"], ["2"], ["m"], ["mm2"], ["4"], ["4mm"], ["3"],     # 2D
     ["3m1", "31m"],       # C₃ᵥ setting variations
     ["6"], ["6mm"]),
    (["1"], ["-1"], ["2"], ["m"], ["2/m"], ["222"], ["mm2"],  # 3D
     ["mmm"], ["4"], ["-4"], ["4/m"], ["422"], ["4mm"], 
     ["-42m", "-4m2"],     # D₂d setting variations
     ["4/mmm"], ["3"], ["-3"],
     ["312", "321"],       # D₃ setting variations  (hexagonal axes)
     ["3m1", "31m"],       # C₃ᵥ setting variations (hexagonal axes)
     ["-31m", "-3m1"],     # D₃d setting variations (hexagonal axes)
     ["6"], ["-6"], ["6/m"], ["622"], ["6mm"], 
     ["-62m", "-6m2"],     # D₃ₕ setting variations
     ["6/mmm"], ["23"], ["m-3"], ["432"], ["-43m"], ["m-3m"])
)
# a flat tuple-listing of all the iuc labels in PG_NUM2IUC; sliced across dimensions
const PG_IUCs = map(x->collect(Iterators.flatten(x)), PG_NUM2IUC)
# a tuple of ImmutableDicts, giving maps from iuc label to point group number
const PG_IUC2NUM = tuple([ImmutableDict([lab=>findfirst(∋(lab), PG_NUM2IUC[D])
                           for lab in PG_IUCs[D]]...) for D in (1,2,3)]...)
# The IUC notation for point groups can be mapped to the Schoenflies notation, but the 
# mapping is not one-to-one but rather one-to-many; e.g. 3m1 and 31m maps to C₃ᵥ but 
# correspond to different axis orientations. 
# When there is a choice of either hexagonal vs. rhombohedral or unique axes b vs unique
# axes a/c we choose hexagonal and unique axes b, respectively.
const PG_IUC2SCHOENFLIES = ImmutableDict(
    "1"     => "C₁",   "-1"   => "Cᵢ",
    "2"     => "C₂",   "m"    => "Cₛ",   "2/m"  => "C₂ₕ",  # unique axes b setting
    "222"   => "D₂",   "mm2"  => "C₂ᵥ",  "mmm"  => "D₂ₕ",  "4"    => "C₄",
    "-4"    => "S₄",   "4/m"  => "C₄ₕ",  "422"  => "D₄",   "4mm"  => "C₄ᵥ",
    "-42m"  => "D₂d",  "-4m2" => "D₂d",  # D₂d setting variations
    "4/mmm" => "D₄ₕ",  "3"    => "C₃",   "-3"   => "C₃ᵢ",
    "312"   => "D₃",   "321"  => "D₃",   # D₃ setting variations  (hexagonal axes)
    "3m1"   => "C₃ᵥ",  "31m"  => "C₃ᵥ",  # C₃ᵥ setting variations (hexagonal axes)
    "-31m"  => "D₃d",  "-3m1" => "D₃d",  # D₃d setting variations (hexagonal axes)
    "6"     => "C₆",   "-6"   => "C₃ₕ",  "6/m"  => "C₆ₕ",  "622"  => "D₆",
    "6mm"   => "C₆ᵥ",
    "-62m"  => "D₃ₕ",  "-6m2" => "D₃ₕ",  # D₃ₕ setting variations
    "6/mmm" => "D₆ₕ",  "23"   => "T",
    "m-3"   => "Tₕ",   "432"  => "O",    "-43m" => "Td",   "m-3m" => "Oₕ"
)

# ===== METHODS =====

# --- Notation ---
function pointgroup_iuc2num(iuclab::String, D::Integer)
    pgnum = get(PG_IUC2NUM[D], iuclab, nothing)
    if pgnum === nothing
        throw(DomainError(iuclab, "invalid point group IUC label"))
    else
        return pgnum
    end
end
schoenflies(pg::PointGroup) = PG_IUC2SCHOENFLIES[iuc(pg)]

@inline function pointgroup_num2iuc(pgnum::Integer, Dᵛ::Val{D}, setting::Integer) where D
    @boundscheck D ∉ (1,2,3) && _throw_invalid_dim(D)
    @boundscheck 1 ≤ pgnum ≤ length(PG_NUM2IUC[D]) || throw(DomainError(pgnum, "invalid pgnum; out of bounds of Crystalline.PG_NUM2IUC"))
    iucs = @inbounds PG_NUM2IUC[D][pgnum]
    @boundscheck 1 ≤ setting ≤ length(iucs) || throw(DomainError(setting, "invalid setting; out of bounds of Crystalline.PG_NUM2IUC[pgnum]"))
    return @inbounds iucs[setting]
end

# --- POINT GROUPS VS SPACE & LITTLE GROUPS ---
function find_parent_pointgroup(g::AbstractGroup{D}) where D
    # Note: this method will only find parent point groups with the same setting (i.e. 
    #       basis) as `g` *and* with the same operator sorting. From a more general
    #       perspective, one might be interested in finding any isomorphic parent point
    #       group (that is achieved by `find_isomorphic_parent_pointgroup`).
    xyzt_pgops = sort!(xyzt.(pointgroup(g)))

    @inbounds for iuclab in PG_IUCs[D]
        P = pointgroup(iuclab, D)
        if sort!(xyzt.(P)) == xyzt_pgops # the sorting/xyzt isn't strictly needed; belts & buckles...
            return P
        end
    end

    return nothing
end

# --- POINT GROUP IRREPS ---
_unmangle_pgiuclab(iuclab) = replace(iuclab, "/"=>"_slash_")

# loads 3D point group data from a .jld2 file (by default, the one in `PGIRREPS_JLDFILE`)
function _load_pgirreps_data(iuclab::String, jldfile::JLD2.JLDFile=PGIRREPS_JLDFILE[])
    jldgroup = jldfile[_unmangle_pgiuclab(iuclab)]
    matrices::Vector{Vector{Matrix{ComplexF64}}} = jldgroup["matrices"]
    realities::Vector{Int8}                      = jldgroup["realities"]
    cdmls::Vector{String}                        = jldgroup["cdmls"]

    return matrices, realities, cdmls
end

# 3D
"""
    pgirreps(iuclab::String, ::Val{D}=Val(3), ::Val{S}=Val(false); mulliken::Bool=false)
    pgirreps(iuclab::String, D::Integer, spinful::Bool=false; mulliken::Bool=false)

Return the (crystallographic) point group irreps of the IUC label `iuclab` of dimension `D`
as a `Vector{PGIrrep{D}}`.

If `S` (or `spinful`) is `true`, the double-valued irreps of the double point group are
returned instead, as `DPGIrrep{D}`s (currently available in 3D only). Their labels are the
CDML labels with an appended `ˢ` (e.g., `"Γ₆ˢ"`); Mulliken labels are not yet available for
them. The single-valued irreps of a double group coincide with those of the ordinary group
and are not included.

See `Crystalline.PG_IUC2NUM[D]` for possible IUC labels in dimension `D`.

## Notation

The irrep labelling follows the conventions of CDML [^1] [which occasionally differ from
those in e.g. Bradley and Cracknell, *The Mathematical Theory of Symmetry in Solids* 
(1972)].

To use Mulliken ("spectroscopist") irrep labels instead, set the keyword argument
`mulliken = true` (default, `false`). See also [`mulliken`](@ref).

## Data sources

The data is sourced from the Bilbao Crystallographic Server [^2]. If you are using this 
functionality in an explicit fashion, please cite the original reference [^3].

## References
[^1]: Cracknell, Davies, Miller, & Love, Kronecher Product Tables 1 (1979).

[^2]: Bilbao Crystallographic Database's
      [Representations PG program](https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=spg).

[^3]: Elcoro et al., 
      [J. of Appl. Cryst. **50**, 1457 (2017)](https://doi.org/10.1107/S1600576717011712)
"""
function pgirreps(
    iuclab::String,
    Dᵛ::Val{3}=Val(3),
    spinfulᵛ::Val{S}=Val(false);
    mulliken::Bool=false
) where S
    pg = pointgroup(iuclab, Dᵛ, spinfulᵛ) # operations

    matrices, realities, cdmls = _load_pgirreps_data(iuclab, _pgirreps_jldfile(spinfulᵛ))
    S && mulliken && error("Mulliken labels are not yet available for double-valued irreps")
    pgirlabs = !mulliken ? cdmls : _mulliken.(Ref(iuclab), cdmls, false)
    
    return Collection(_pgirrep.(pgirlabs, Ref(pg), matrices, Reality.(realities)))
end
_pgirrep(cdml, pg::PointGroup{D}, P, reality) where D = PGIrrep{D}(cdml, pg, P, reality)

_pgirreps_jldfile(::Val{false}) = PGIRREPS_JLDFILE[]
function _pgirreps_jldfile(#=Val{S}=# ::Val{true})
    isassigned(DPGIRREPS_JLDFILE) || error("spinful point group irrep data file not found")
    return DPGIRREPS_JLDFILE[]
end
# 2D
function pgirreps(iuclab::String, ::Val{2}, ::Val{false}=Val(false); mulliken::Bool=false)
    pg = pointgroup(iuclab, Val(2)) # operations

    # Because the operator sorting and setting is identical* between the shared point groups
    # of 2D and 3D, we can just do a whole-sale transfer of shared irreps from 3D to 2D.
    # (*) Actually, "2" and "m" have different settings in 2D and 3D; but they just have two
    #     operators and irreps each, so the setting difference doesn't matter.
    #     That the settings and sorting indeed agree between 2D and 3D is tested in 
    #     scripts/compare_pgops_3dvs2d.jl
    matrices, realities, cdmls = _load_pgirreps_data(iuclab)
    pgirlabs = !mulliken ? cdmls : _mulliken.(Ref(iuclab), cdmls, false)
    
    return Collection(PGIrrep{2}.(pgirlabs, Ref(pg), matrices, Reality.(realities)))
end
# 1D
function pgirreps(iuclab::String, ::Val{1}, ::Val{false}=Val(false); mulliken::Bool=false)
    pg = pointgroup(iuclab, Val(1))
    # Situation in 1D is sufficiently simple that we don't need to bother with loading from 
    # a disk; just branch on one of the two possibilities
    if iuclab == "1"
        matrices = [[fill(one(ComplexF64), 1, 1)]]
        cdmls    = ["Γ₁"]
    elseif iuclab == "m"
        matrices = [[fill(one(ComplexF64), 1, 1),  fill(one(ComplexF64), 1, 1)], # even
                    [fill(one(ComplexF64), 1, 1), -fill(one(ComplexF64), 1, 1)]] # odd
        cdmls    = ["Γ₁", "Γ₂"]
    else
        throw(DomainError(iuclab, "invalid 1D point group IUC label"))
    end
    pgirlabs = !mulliken ? cdmls : _mulliken.(Ref(iuclab), cdmls, false)
    
    return Collection(PGIrrep{1}.(pgirlabs, Ref(pg), matrices, REAL))
end
# invalid dimensions, and spinful irreps in 1D or 2D
function pgirreps(iuclab::String, ::Val{D}, ::Val=Val(false); kws...) where D
    return D ∈ (1,2) ? _only_3d(D) : _throw_invalid_dim(D)
end
function pgirreps(iuclab::String, D::Integer, spinful::Bool=false; kws...)
    return pgirreps(iuclab, Val(D), Val(spinful); kws...)
end
function pgirreps(
    pgnum::Integer,
    Dᵛ::Val{D}=Val(3),
    spinfulᵛ::Val=Val(false);
    setting::Integer=1,
    kws...
) where D
    iuc = pointgroup_num2iuc(pgnum, Dᵛ, setting)
    return pgirreps(iuc, Dᵛ, spinfulᵛ; kws...)
end
function pgirreps(pgnum::Integer, D::Integer, spinful::Bool=false; kws...)
    return pgirreps(pgnum, Val(D), Val(spinful); kws...) :: Collection{<:AbstractPGIrrep}
end

function ⊕(pgir1::PGIrrep{D}, pgir2::PGIrrep{D}) where D
    if label(group(pgir1)) ≠ label(group(pgir2)) || order(pgir1) ≠ order(pgir2)
        error("The direct sum of two PGIrreps requires identical point groups")
    end
    
    pgirlab = label(pgir1)*"⊕"*label(pgir2)
    g   = group(pgir1)
    T   = eltype(eltype(pgir1.matrices))
    z12 = zeros(T, irdim(pgir1), irdim(pgir2))
    z21 = zeros(T, irdim(pgir2), irdim(pgir1))
    matrices = [[m1 z12; z21 m2] for (m1, m2) in zip(pgir1.matrices, pgir2.matrices)]
    reality = UNDEF
    iscorep = pgir1.iscorep || pgir2.iscorep

    return PGIrrep{D}(pgirlab, g, matrices, reality, iscorep)
end

# ---------------------------------------------------------------------------------------- #

"""
    find_isomorphic_parent_pointgroup(g::AbstractVector{SymOperation{D}}) 
                                                    --> PointGroup{D}, Vector{Int}, Bool

Given a group `g` (or a collection of operators, defining a group), identifies a "parent"
point group that is isomorphic to `g`.

Three variables are returned:
- `pg`: the identified "parent" point group, with operators sorted to match the sorting
  of `g`'s operators.
- `Iᵖ²ᵍ`: a permutation vector which transforms the standard sorting of point group
  operations (as returned by [`pointgroup(::String)`](@ref)) to the operator sorting of `g`.
- `equal`: a boolean, identifying whether the point group parts of `g` operations are
  identical (`true`) or merely isomorphic to the point group operations in `g`.
  In practice, this indicates whether `pg` and `g` are in the same setting or not.

## Implementation
The identification is made partly on the basis of comparison of operators (this is is
sufficient for the `equal = true` case) and partly on the basis of comparison of 
multiplication tables (`equal = false` case); the latter can be combinatorially slow if
the sorting of operators is unlucky (i.e., if permutation between sortings in `g` and `pg`
differ by many pairwise permutations).

Beyond mere isomorphisms of multiplication tables, the search also guarantees that all
rotation orders are shared between `pg` and `g`; similarly, the rotation senses (e.g., 4⁺ &
4⁻ have opposite rotation senses or directions) are shared. This disambiguates point groups 
that are intrinsically isomorphic to eachother, e.g. "m" and "-1", but which still differ in
their spatial interpretation.

## Properties
The following properties hold for `g`, `pg`, and `Iᵖ²ᵍ`:
```jl
pg, Iᵖ²ᵍ, equal = find_isomorphic_parent_pointgroup(g)
@assert MultTable(pg) == MultTable(pointgroup(g))
pg′ = pointgroup(label(pg), dim(pg)) # "standard" sorting
@assert pg′[Iᵖ²ᵍ] == pg
```
If `equal = true`, the following also holds:
```jl
pointgroup(g) == operations(pg) == operations(pg′)[Iᵖ²ᵍ]
```

## Example
```jl
sgnum = 141
wp    = wyckoffs(sgnum, Val(3))[end] # 4a Wyckoff position
sg    = spacegroup(sgnum, Val(3))
siteg = sitegroup(sg, wp)
pg, Iᵖ²ᵍ, equal = find_isomorphic_parent_pointgroup(siteg)
```
"""
function find_isomorphic_parent_pointgroup(g::AbstractVector{SymOperation{D}}) where D
    Nᵍ = length(g)
    Nᵍ == 1 && isone(first(g)) && return pointgroup("1", Val(D)), [1], true

    # a site group is always similar to a point group (i.e. a group w/o translations); it
    # can be transformed to that setting by transforming each element g∈`siteg` as t⁻¹∘g∘t
    # where t denotes translation by the site group's wyckoff position. In practice, we can
    # just "drop the translation" parts to get the same result though; this is faster
    g′ = pointgroup(g) # get rid of any repeated point group operations, if they exist
    ordersᵍ = rotation_order_and_sense.(g′)

    # TODO: For the isomorphism search below, it is _NOT_ correct to assume we can always
    #       group up the operations of `g′` by their _sense_ as well. In fact, the best we
    #       can do is to group them by their conjugacy classes (see `classes(g′)`).
    #       For the "standard-input" point groups (as tabulated by `pointgroup`), it happens
    #       to work out alright to add the wrong "sense-grouping" assumption - but it
    #       quickly goes awry. E.g., if we "primitivize" a set of such point groups, through
    #       e.g., the I centering (which has a negative determinant, crucially), the
    #       assumption can lead us sufficiently wrong that we can't find the correct
    #       isomorphic group.
    #       On the other hand, the rotation+sense grouping is usually more effective, since
    #       it produces more distinct groups that just using the conjugacy classes would.

    # first sorting step: by rotation order
    Iᵍ = sortperm(ordersᵍ)
    permute!(ordersᵍ, Iᵍ) 
    permute!(g′, Iᵍ)

    # identify rotation orders as indexing groups
    I_groups = _find_equal_groups_in_sorted(ordersᵍ)

    # --------------------------------------------------------------------------------------
    # first, we check if we literally have the exact same group operations; in that case
    # we can check for group equality rather than isomorphism, and get combinatorial speedup
    # since we don't have to attempt a search over permutations of multiplication tables
    Iᵖ      = Vector{Int}(undef, Nᵍ)
    ordersᵖ = Vector{Tuple{Int, Int}}(undef, Nᵍ)
    for iuclab in PG_IUCs[D]
        pg = _fast_path_necessary_checks!(iuclab, ordersᵖ, Iᵖ, Nᵍ, ordersᵍ, Val(D))
        isnothing(pg) && continue

        # check if there is a permutation of operators that makes `g` and `pg` equal
        bool, permᵖ²ᵍ = _has_equal_pg_operations(g′, pg, I_groups, nothing)
        bool || continue

        # we have a match: compute "unwinding" compound permutation, s.t. `g ~ pg[Iᵖ²ᵍ]`
        Iᵖ²ᵍ = invpermute!(permute!(Iᵖ, permᵖ²ᵍ), Iᵍ) # Iᵖ[permᵖ²ᵍ][invperm(Iᵍ)], but faster
        invpermute!(permute!(pg.operations, permᵖ²ᵍ), Iᵍ) # sort `pg` to match `g`

        return pg, Iᵖ²ᵍ, true
    end

    # --------------------------------------------------------------------------------------
    # if there is a rotation order with only one element, which is not simply the identity,
    # then we can compare that element directly with an associated point group element,
    # even if we are in different settings: whatever the mapping is between them, it will be
    # uniquely defined (up to a sign); so if this the case, we can do the "exact check"
    # also for groups that are related by an orthogonal transformation
    unique_op_idx_into_I_groups = findfirst(I_groups) do I_group
        length(I_group) == 1 || return false
        isone(g′[I_group[1]]) && return false
        return true
    end
    if !isnothing(unique_op_idx_into_I_groups)
        unique_op_idx = only(I_groups[something(unique_op_idx_into_I_groups)]) # into `g′`
        unique_opᵍ = g′[unique_op_idx]
        unique_Wᵍ  = rotation(unique_opᵍ)
        for iuclab in PG_IUCs[D]
            pg = _fast_path_necessary_checks!(iuclab, ordersᵖ, Iᵖ, Nᵍ, ordersᵍ, Val(D))
            isnothing(pg) && continue
    
            # compute the unique transformation between the two group's "unique" elements
            # s.t. `unique_Wᵖ = transformation * unique_Wᵍ * transformation⁻¹`
            unique_Wᵖ = rotation(pg[unique_op_idx])
            transformation = find_similarity_transform(unique_Wᵍ, unique_Wᵖ)
            # check if there is a permutation of operators that makes `g` and `pg` equal
            bool, permᵖ²ᵍ = _has_equal_pg_operations(g′, pg, I_groups, transformation)
            bool || continue
            
            # we have a match: compute "unwinding" compound permutation, s.t. `g ~ pg[Iᵖ²ᵍ]`
            Iᵖ²ᵍ = invpermute!(permute!(Iᵖ, permᵖ²ᵍ), Iᵍ) # = Iᵖ[permᵖ²ᵍ][invperm(Iᵍ)]
            invpermute!(permute!(pg.operations, permᵖ²ᵍ), Iᵍ) # sort `pg` to match `g`
            
            return pg, Iᵖ²ᵍ, false
        end
    end

    # --------------------------------------------------------------------------------------
    # if we arrived here, the group is not equal to any group in our "storage", so we must
    # make do with finding an isomorphic group; this is potentially much more demanding
    # because we have to check all (meaningful) permutations of the multiplication tables;
    # for some groups, this is practically impossible (even with the tricks we play to
    # reduce the combinatorial explosion into a product of smaller combinatorial problems).
    # fortunately, for the space groups we are interested in, the earlier equality check
    # takes care of bad situations like that.
    # we deliberately split this into a separate loop (even though it incurs a bit of
    # repetive work) because we'd prefer to return an identical group rather than an
    # isomorphic group, if that is at all possible
    mtᵍ = MultTable(g′).table
    for iuclab in PG_IUCs[D]
        # NB: in principle, this could probably be sped up a lot by using subgroup relations
        #     to "group up" meaningful portions (in addition to grouping by rotation order)
        pg = _fast_path_necessary_checks!(iuclab, ordersᵖ, Iᵖ, Nᵍ, ordersᵍ, Val(D))
        isnothing(pg) && continue

        # create `pg`'s multiplication table, "order-sorted"
        mtᵖ = MultTable(pg).table
        
        # rigorous, slow, combinatorial search
        P = rigorous_isomorphism_search(I_groups, mtᵍ, mtᵖ)

        if P !== nothing # found isomorphic point group!
            # compute "combined/unwound indexing" Iᵖ²ᵍ of Iᵍ, Iᵖ, and P, where Iᵖ²ᵍ has
            # produces g ~ pg[Iᵖ²ᵍ] in the sense that they have the same character table
            # and same rotation orders. We can get that by noting that we already have
            # g[Iᵍ] ~ pg[Iᵖ[P]], and then `invperm` gets us the rest of the way.
            Iᵖ²ᵍ = invpermute!(permute!(Iᵖ, P), Iᵍ) # = Iᵖ[P][invperm(Iᵍ)], but faster
            invpermute!(permute!(pg.operations, P), Iᵍ) # sort pg to match g
            return pg, Iᵖ²ᵍ, false
        end
    end

    error("could not find an isomorphic parent point group")
end

"""
    _find_equal_groups_in_sorted(v::AbstractVector) --> Vector{UnitRange}

Returns indices into groups of equal values in `v`. Input `v` *must* be sorted so that
identical values are adjacent.

## Example
```jl
_find_equal_groups_in_sorted([-4,1,1,3,3,3]) # returns [1:1, 2:3, 4:6]
```
"""
function _find_equal_groups_in_sorted(v::AbstractVector)
    vᵘ = unique(v)
    Nᵘ, N = length(vᵘ), length(v)
    idx_groups = Vector{UnitRange{Int}}(undef, Nᵘ)
    start = stop = 1
    for (i,o) in enumerate(vᵘ)
        while stop < N && v[stop+1] == o
            stop += 1
        end
        idx_groups[i] = start:stop
        stop += 1
        start = stop
    end
    return idx_groups
end

# checks two necessary (but not sufficient) conditions for a reference (point) group ("g";
# specified by its order `Nᵍ` and sorted rotation-order-&-sense vector `ordersᵍ`) to be 
# isomorphic to the conventional point group ("p") with label `iuclab` and dimension `D`; 
# the rotation-order-&-sense values of `pg` are written into `ordersᵖ` which is mutated
# (provided the first check succeeds); similarly, the sorting-permutation is written into
# `Iᵖ` in mutating fashion
# if the conditions are not met, we can bail early to save time, returning `nothing` - but
# if they are met, we return the rotation-order-&-sense sorted point group operations `pg`
# associated with `iuclab` as well as a permutation vector `Iᵖ` that implements the sorting
function _fast_path_necessary_checks!(
        iuclab::AbstractString,
        ordersᵖ::Vector{Tuple{Int, Int}}, # ← mutated: must be `similar` to `ordersᵍ`
        Iᵖ::Vector{Int},                  # ← mutated: must have same length as `ordersᵍ`
        Nᵍ::Integer,
        ordersᵍ::Vector{Tuple{Int, Int}}, 
        Dᵛ::Val{D} where D
        )
    
    # fast-path check: number of group elements (order) must agree
    PG_ORDERs[iuclab] == Nᵍ || return nothing

    # go ahead and load the point group operations for closer inspection
    pg = pointgroup(iuclab, Dᵛ)

    # compare in rotation-order-&-sense-sorting + check if they agree between "g" and "p"
    ordersᵖ .= rotation_order_and_sense.(pg)
    sortperm!(Iᵖ, ordersᵖ)
    permute!(ordersᵖ, Iᵖ)
    ordersᵍ == ordersᵖ || return nothing

    # checks passed; return point group rotation-order-&-sense-sorted
    permute!(pg.operations, Iᵖ)
    return pg
end

# check if the group `g′` has the same operations as the point group with label `iuclab`,
# possibly with in a different sorting
function _has_equal_pg_operations(
        g′::AbstractVector{SymOperation{D}},
        pg::AbstractVector{SymOperation{D}},
        I_groups::AbstractVector{<:AbstractVector{<:Integer}},
        transformation::Union{Nothing, AbstractMatrix{<:Real}}
        # ↑ assumes mapping `Ref(pg) = transformation .* Ref(g′) .* inv(transformation)`
    ) where D

    # check if the group operations literally are equal, but (maybe) just shuffled
    permᵖ²ᵍ = Vector{Int}(undef, sum(length, I_groups))
    for I_group in I_groups
        for i in I_group
            opᵢ = if isnothing(transformation)
                g′[i]
            else
                transform(g′[i], transformation) # `transformation * g′[i] * transformation⁻¹`
            end
            idx = findfirst(≈(opᵢ), (@view pg[I_group]))
            idx === nothing && return false, permᵖ²ᵍ # ... not equal
            permᵖ²ᵍ[i] = idx + I_group[1] - 1
        end
    end
    
    return true, permᵖ²ᵍ
end

"""
    permute_multtable(mt::Matrix{Int}, P::AbstractVector{Int}) -> Matrix{Int}

Returns a multiplication table derived from `mt` but under a reordering permutation `P` of
the operator indices in `mt`. In practice, this corresponds to a permutation of the rows
and columns of `mt` as well as a subsequent relabelling of elements via `P`.

## Example
```jl
julia> pg  = pointgroup("-4m2", Val(3))
julia> mt  = MultTable(pg)
julia> P   = [2,3,1,4,5,8,6,7] # permutation of operator indices
julia> mt′ = permute_multtable(mt.table, P)
julia> mt′ == MultTable(pg[P]).table
```
"""
function permute_multtable(
    mt::Matrix{Int},
    P::AbstractVector{Int},
    # work-arrays, for reducing allocations
    mt′::Matrix{Int} = similar(mt, Int, (length(P), length(P))), # output matrix
    invP::AbstractVector{Int} = Vector{Int}(undef, length(P))
)
    # The implementation below is equivalent to (but much faster than)
    #       `replace!(mt[P, P], (P .=> Base.OneTo(length(P)))...)`
    N = length(P)
    @assert extrema(P) == (1, N)
    size(mt, 1) ≥ N && size(mt, 2) ≥ N || throw(DimensionMismatch("`mt` must be a square matrix whose size exceeds the length of P"))
    size(mt′) == (N, N) || throw(DimensionMismatch("`mt′` must be a square matrix whose size equals the length of P"))
    length(invP) == N || throw(DimensionMismatch("`invP` must have the same length as P"))
    for i in 1:N # compute the inverse permutation vector for value reassignment
        @inbounds invP[P[i]] = i
    end
    # permute row/cols of multiplication table & assign permuted values also
    for i in 1:N, j in 1:N
        v = @inbounds mt[P[i], P[j]] # get the old value at the permuted indices
        if v > N # if `v` exceeds `N`, we keep whatever value it has
            @inbounds mt′[i, j] = v
        else
            @inbounds mt′[i, j] = invP[v]
        end
    end
    return mt′
end

# this is a rigorous search; it checks every possible permutation of each "order group".
# it is much faster (sometimes by ~14 orders of magnitude) than a naive permutation search
# which doesn't divide into "order groups" first, but can still be intractably slow (e.g.
# for site groups of space group 216).
function rigorous_isomorphism_search(I_orders::AbstractVector{<:AbstractVector{Int}}, 
            mtᵍ::Matrix{Int}, mtᵖ::Matrix{Int}, Jᵛ::Val{J}=Val(length(I_orders))) where J
    # we're using `j` as the rotation-order-group index below
    J == length(I_orders) || throw(DimensionMismatch("static parameter and length must match"))
    Ns = ntuple(j->factorial(length(I_orders[j])), Jᵛ)

    # conceptually, we now want to loop over all possible products of permutations of the
    # elements of `I_orders` (each a vector); effectively, this is a kroenecker product of 
    # permutations. In total, there are `prod(Iⱼ -> factorial(length(Iⱼ), I_orders)` such
    # permutations; this is much less than a brute-force search over permutations that
    # doesn't exploit groupings (with `factorial(sum(length, I_orders))`), but it can still
    # be prohibively slow, even overwhelmingly so. To improve matters more, we arrange the
    # search in such a way that we check the permutations in "blocks" of equal order. This
    # is because there is no point in exploring permutations of `I_orders[j′]` for `j′ > j`
    # if `j`th permutation of `I_orders[j]` wasn't matching up with `mtᵍ`; exploiting this
    # can give enormous scaling improvements since it approximately turns the `prod(...`
    # above into a (complicated) `sum(...`. The implementation is most naturally done by
    # recursively going into the blocks, and then backtracking if the latest block
    # permutation was bad: this is done in `_search_permuted_blocks_recursively!` below
    Nops = sum(length, I_orders)
    P = collect(1:Nops) # permutation vector to be mutated
    P, success = _search_permuted_blocks_recursively!(P, 1, Ns, I_orders, mtᵍ, mtᵖ)
    if success
        return P # found a valid isomorphism (the permutation) between mtᵍ and mtᵖ
    end
    # TODO: We could probably still improve matters a bit by sorting `I_orders` (and the
    #       associated operations in `mtᵍ` and `mtᵖ`) in such a way that blocks with fewest
    #       permutations are checked first; for now, it's still just by rotation order

    return nothing # failed to find an isomorphism (because there wasn't any!)
end
# TODO: We could probably do much better than this by adapting the below to compute a
#       transformation matrix for each "locked-in" block (probably only need one nontrivial
#       block?) and then use that matrix to transform all the _operations_ in e.g., `g` to
#       `p`: if the operations in the are identical (up to sorting) after the transformation
#       then we have demonstrated that they are isomorphic (and we can also immediately) do
#       find the sorting permutation between them.

function _search_permuted_blocks_recursively!(
    P::AbstractVector{Int}, # mutated permutation vector
    j::Int, # current block index
    Ns::NTuple{N, Int}, # permutation lengths for each block
    I_orders::AbstractVector{<:AbstractVector{Int}}, # order groups
    mtᵍ::Matrix{Int}, # multiplication table of group g
    mtᵖ::Matrix{Int}, # multiplication table of group p
) where N
    Iⱼ = I_orders[j]
    idxs = Base.OneTo(maximum(Iⱼ))
    mtᵍⱼ = @view mtᵍ[idxs, idxs]
    Nⱼ = Ns[j]
    mtᵖⱼ′ = Matrix{Int}(undef, length(idxs), length(idxs)) # preallocate
    invPⱼ = Vector{Int}(undef, length(idxs))
    for nⱼ in 1:Nⱼ
        # permute `P[Iⱼ]` in-place, to the `nⱼ`th permutation of the current block: the line
        # below is equivalent to `pⱼ = nthperm(Iⱼ, nⱼ); P[Iⱼ] .= pⱼ`, but it doesn't
        # allocate anything
        nthperm!(copyto!((@view P[Iⱼ]), Iⱼ), nⱼ) # update with the current block permutation
        mtᵖⱼ′ = permute_multtable(mtᵖ, (@view P[idxs]), mtᵖⱼ′, invPⱼ) # permute the multiplication table
        _is_good_perm(mtᵍⱼ, mtᵖⱼ′) || continue # if the current block doesn't match, skip
        if j == length(I_orders) # if we are at the last block, we have a match
            return P, true # return the permutation vector that matches mtᵍ and mtᵖ
        else # if not, recurse into the next block
            P, success = _search_permuted_blocks_recursively!(P, j+1, Ns, I_orders, mtᵍ, mtᵖ)
            success && return P, true # the recursion returned a valid permutation, return
        end
    end
    return P, false # if no good permutation was found, return `false` for `success`
end

function _is_good_perm(mtᵍⱼ::AbstractMatrix{Int}, mtᵖⱼ′::AbstractMatrix{Int})
    # check if the multiplication tables match under the current permutation (but we cannot
    # at this point generally require `mtᵍ == mtᵖ` because they may still differ any 
    # "yet-to-be-permuted" situations; i.e., any index above `size(mtᵍ, 1)` is allowed to
    # differ)
    N = size(mtᵍⱼ, 1)
    for (vᵍ, vᵖ′) in zip(mtᵍⱼ, mtᵖⱼ′)
        (vᵍ > N && vᵖ′ > N) && continue # "not-yet-permuted" indices are allowed to differ
        vᵍ ≠ vᵖ′ && return  false
    end
    return true # multiplication tables match up to current permutation indices
    # Equiv. to:      x = copy(mtᵍⱼ); x[x .> N] .= 0
    # (but faster)    y = copy(mtᵖⱼ′); y[y .> N] .= 0
    #                 return x == y
end

function rotation_order_and_sense(op::SymOperation{D}) where D
    r = rotation_order(op)
    D == 1 && return (r, 0) # no "rotation sense" in 1D; `0` as sentinel

    o = abs(r)
    if o ≤ 2 # no rotation sense for identity, mirrors, 2-fold rotation, or inversion
        return (r, 0) # use `0` as sentinel
    end
    
    # henceforth, `o > 2`, and it is possible to define a ±1 (left/right) sense relative to
    # an axis; the calculation below is borrowed from `seitz` (TODO: consolidate)
    W = rotation(op)
    detW = det(W)
    if D == 2 # augment to be 3×3 for ease of implementation
        W = @inbounds SMatrix{3,3,Float64,9}( # build by column (= [W zeros(2); 0 0 1])
                W[1], W[2], 0.0, W[3], W[4], 0.0, 0.0, 0.0, 1.0)
    end
    
    # --- rotation axis ---
    u = if D == 2
        SVector{3,Int}(0, 0, 1)
    else # D == 3
        rotation_axis_3d(rotation(op), detW, o)
    end
    
    # --- rotation sense ---
    # ±-rotation sense is determined from sign of det(𝐙) where  𝐙 ≡ [𝐮|𝐱|det(𝐖)𝐖𝐱]
    # with 𝐱 an arbitrary vector that is not parallel to 𝐮. [ITA6 vol. A, p. 16, sec.
    # 1.2.2.4(1)(c)]
    x = rand(-1:1, SVector{3, Int})
    while iszero(x×u) # check that generated 𝐱 is not parallel to 𝐮 (if it is, 𝐱×𝐮 = 0)
        x = rand(-1:1, SVector{3, Int})
        iszero(u) && error("rotation axis has zero norm; input is likely invalid")
    end
    Z = hcat(u, x, detW*(W*x))
    sense_float = sign(det(Z))
    sense = round(Int, sense_float)
    sense_float ≈ sense || error("rotation sense is not an integer; unexpected failure")
    
    return (r, sense)
end

#       find_similarity_transform(X, Y)  -->  Matrix{<:Real}
# Given two similar (diagonalizable) matrices, `X` and `Y`, find an orthogonal
# transformation matrix `P` s.t. `Y = P⁻¹XP`. Errors if `X` and `Y` are not similar.
# The idea works like this: denote by `Λ` the (identical, by necessity of being similar)
# eigenvalue matrix of `X` and `Y` and by `Vˣ` and `Vʸ` the associated (invertible)
# eigenvector matrices. Then `X = VˣΛ(Vˣ)⁻¹` and `Y = VʸΛ(Vʸ)⁻¹` and equivalently
# `(Vʸ)⁻¹YVʸ = Λ = (Vˣ)⁻¹XVˣ`; left- & right-multiplying this relation by `Vʸ` & `(Vʸ)⁻¹`, 
# respectively, gives `Vʸ(Vʸ)⁻¹YVʸ(Vʸ)⁻¹ = Y = Vʸ(Vˣ)⁻¹XVˣ(Vʸ)⁻¹ = [Vʸ(Vˣ)⁻¹]X[Vˣ(Vʸ)⁻¹]`,
# which, finally, by comparison with `Y = P⁻¹XP` lets us identify `P = Vʸ(Vˣ)⁻¹`.
# If `X` and `Y` are real symmetric matrices, we have `(Vˣʸ)⁻¹ = (Vˣʸ)ᵀ` (i.e., orthogonal
# eigenvector matrices), in which case the transformation is also orthogonal (`Pᵀ=P⁻¹`).
function find_similarity_transform(X, Y)
    eˣ  = eigen(X); eʸ = eigen(Y)
    if !(eˣ.values ≈ eʸ.values)
        error("matrices `X` and `Y` are not similar; a basis change between `X` and `Y` does not exist")
    end
    Vˣ = eˣ.vectors; Vʸ = eʸ.vectors   
    P = Vˣ*inv(Vʸ)

    if P isa Matrix{<:Real}
        # oftentimes, the eigendecomposition will simply be real already
        return P
    else
        # we don't want to be working with complex matrices, and physically, this should never
        # occur, so we check here and error if we have non-negligible imaginary entries
        if sum(abs∘imag, P) > DEFAULT_ATOL
            error("non-neglible imaginary parts in transformation matrix")
        end
        return real(P)
    end
end