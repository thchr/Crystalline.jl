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

# loads 3D point group data from the .jld2 file opened in `PGIRREPS_JLDFILE`
function _load_pgirreps_data(iuclab::String)
    jldgroup = PGIRREPS_JLDFILE[][_unmangle_pgiuclab(iuclab)] 
    matrices::Vector{Vector{Matrix{ComplexF64}}} = jldgroup["matrices"]
    realities::Vector{Int8}                      = jldgroup["realities"]
    cdmls::Vector{String}                        = jldgroup["cdmls"]

    return matrices, realities, cdmls
end

# 3D
"""
    pgirreps(iuclab::String, ::Val{D}=Val(3); mullikken::Bool=false) where D ∈ (1,2,3)
    pgirreps(iuclab::String, D; mullikken::Bool=false)

Return the (crystallographic) point group irreps of the IUC label `iuclab` of dimension `D`
as a `Vector{PGIrrep{D}}`.

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
function pgirreps(iuclab::String, ::Val{3}=Val(3); mulliken::Bool=false)
    pg = pointgroup(iuclab, Val(3)) # operations

    matrices, realities, cdmls = _load_pgirreps_data(iuclab)
    pgirlabs = !mulliken ? cdmls : _mulliken.(Ref(iuclab), cdmls, false)
    
    return IrrepCollection(PGIrrep{3}.(pgirlabs, Ref(pg), matrices, Reality.(realities)))
end
# 2D
function pgirreps(iuclab::String, ::Val{2}; mulliken::Bool=false)
    pg = pointgroup(iuclab, Val(2)) # operations

    # Because the operator sorting and setting is identical* between the shared point groups
    # of 2D and 3D, we can just do a whole-sale transfer of shared irreps from 3D to 2D.
    # (*) Actually, "2" and "m" have different settings in 2D and 3D; but they just have two
    #     operators and irreps each, so the setting difference doesn't matter.
    #     That the settings and sorting indeed agree between 2D and 3D is tested in 
    #     scripts/compare_pgops_3dvs2d.jl
    matrices, realities, cdmls = _load_pgirreps_data(iuclab)
    pgirlabs = !mulliken ? cdmls : _mulliken.(Ref(iuclab), cdmls, false)
    
    return IrrepCollection(PGIrrep{2}.(pgirlabs, Ref(pg), matrices, Reality.(realities)))
end
# 1D
function pgirreps(iuclab::String, ::Val{1}; mulliken::Bool=false)
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
    
    return IrrepCollection(PGIrrep{1}.(pgirlabs, Ref(pg), matrices, REAL))
end
pgirreps(iuclab::String, ::Val{D}; kws...) where D = _throw_invalid_dim(D) # if D ∉ (1,2,3)
pgirreps(iuclab::String, D::Integer; kws...) = pgirreps(iuclab, Val(D); kws...)
function pgirreps(pgnum::Integer, Dᵛ::Val{D}=Val(3);
                  setting::Integer=1, kws...) where D
    iuc = pointgroup_num2iuc(pgnum, Dᵛ, setting)
    return pgirreps(iuc, Dᵛ; kws...)
end
function pgirreps(pgnum::Integer, D::Integer; kws...)
    return pgirreps(pgnum, Val(D); kws...) :: IrrepCollection{<:PGIrrep}
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
rotation orders are shared between `pg` and `g`. This disambiguates point groups that are
intrinsically isomorphic to eachother, e.g. "m" and "-1", but which still differ in their
spatial interpretation.

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
    ordersᵍ = rotation_order.(g′)

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
    ordersᵖ = Vector{Int}(undef, Nᵍ)
    for iuclab in PG_IUCs[D]
        PG_ORDERs[iuclab] == Nᵍ || continue
        pg = pointgroup(iuclab, Val(D))

        # sorting by rotation order & check if rotation orders agree
        ordersᵖ .= rotation_order.(pg)
        sortperm!(Iᵖ, ordersᵖ)
        permute!(ordersᵖ, Iᵖ)
        ordersᵍ == ordersᵖ || continue

        # check if the group operations literally are equal, but (maybe) just shuffled
        permute!(pg.operations, Iᵖ)
        P = Vector{Int}(undef, sum(length, I_groups))
        equal = true
        for I_group in I_groups
            for i in I_group
                idx = findfirst(==(g′[i]), (@view pg[I_group]))
                idx === nothing && (equal = false; break) # ... not equal
                P[i] = idx + I_group[1] - 1
            end
            equal || break
        end
        if equal
            # "unwind" compound permutation, s.t. g ~ pg[Iᵖ²ᵍ]
            Iᵖ²ᵍ = invpermute!(permute!(Iᵖ, P), Iᵍ) # = Iᵖ[P][invperm(Iᵍ)], but faster
            invpermute!(permute!(pg.operations, P), Iᵍ) # sort pg to match g
            return pg, Iᵖ²ᵍ, true
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
        Crystalline.PG_ORDERs[iuclab] == Nᵍ || continue
        pg = pointgroup(iuclab, Val(D))

        # sorting by rotation order & check if rotation orders agree
        ordersᵖ .= rotation_order.(pg)
        sortperm!(Iᵖ, ordersᵖ)
        permute!(ordersᵖ, Iᵖ)
        ordersᵍ == ordersᵖ || continue

        # create pg's multiplication table, "order-sorted"
        permute!(pg.operations, Iᵖ)
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

"""
    $TYPEDSIGNATURES --> Matrix{Int}

Returns a multiplication table derived from `mt` but under a reordering permutation `P` of
the operator indices in `mt`. In practice, this corresponds to a permutation of the rows
and columns of `mt` as well as a subsequent relabelling of elements via `P`.

## Example
```jl
pg  = pointgroup("-4m2", Val(3))
mt  = MultTable(pg)
P   = [2,3,1,4,5,8,6,7] # permutation of operator indices
mt′ = permute_multtable(mt.table, P)
mt′ == MultTable(pg[P]).table
```
"""
function permute_multtable(mt::Matrix{Int}, P::AbstractVector{Int})
    replace!(mt[P, P], (P .=> Base.OneTo(length(P)))...)
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
    P  = Vector{Int}(undef, sum(length, I_orders))
    # loop over all possible products of permutations of the elements of I_orders (each a 
    # vector); effectively, this is a kroenecker product of permutations. In total, there
    # are `prod(Iⱼ -> factorial(length(Iⱼ), I_orders)` such permutations; this is much less
    # than a brute-force search over permutations that doesn't exploit groupings (with
    # `factorial(sum(length, I_orders))`) - but it can still be overwhelming occassionally
    ns_prev = zeros(Int, J)
    for ns in CartesianIndices(Ns)
        # TODO: This should be done in "blocks" of equal order: i.e., there is no point in
        #       exploring permutations of I_orders[j+1] etc. if the current permutation of
        #       I_orders[j] isn't matching up with mtᵍ - this could give very substantial
        #       scaling improvements since it approximately turns a `prod(...` above into a
        #       `sum(...`
        for j in Base.OneTo(J)
            nⱼ = ns[j]
            if nⱼ != ns_prev[j]
                Iⱼ = I_orders[j]
                P[Iⱼ] .= nthperm(Iⱼ, nⱼ)
                ns_prev[j] = nⱼ
            end
        end
        mtᵖ′ = permute_multtable(mtᵖ, P)
        if mtᵍ == mtᵖ′
            return P # found a valid isomorphism (the permutation) between mtᵍ and mtᵖ
        end
    end
    return nothing # failed to find an isomorphism (because there wasn't any!)
end