# ---------------------------------------------------------------------------------------- #
# CONSTRUCTORS/GETTERS FOR WYCKPOS

"""
$(TYPEDSIGNATURES)

Return the Wyckoff positions of space group `sgnum` in dimension `D` as a 
`Vector{WyckoffPosition{D}`.

The positions are given in the conventional basis setting, following the conventions of the
Bilbao Crystallographic Server (from which the underlying data is sourced [^1]).

## Example
```jldoctest
julia> wps = wyckoffs(16, 2)
4-element Vector{WyckoffPosition{2}}:
 6d: [α, β]
 3c: [1/2, 0]
 2b: [1/3, 2/3]
 1a: [0, 0]
```

## References
[^1]: Aroyo, *et al.*,
      [Z. Kristallogr. **221**, 15-27 (2006)](https://doi.org/0.1524/zkri.2006.221.1.15)
"""
function wyckoffs(sgnum::Integer, ::Val{D}=Val(3)) where D
    strarr = open(joinpath(DATA_DIR, "wyckpos/$(D)d/"*string(sgnum)*".csv")) do io
        DelimitedFiles.readdlm(io, '|', String, '\n')
    end::Matrix{String}
    mults   = parse.(Int, @view strarr[:,1])::Vector{Int}
    letters = only.(@view strarr[:,2])::Vector{Char}
    rvs     = RVec{D}.((@view strarr[:,3]))::Vector{RVec{D}}

    return WyckoffPosition{D}.(mults, letters, rvs)
end
wyckoffs(sgnum::Integer, D::Integer) = wyckoffs(sgnum, Val(D))

# ---------------------------------------------------------------------------------------- #
# METHODS

function compose(op::SymOperation{D}, wp::WyckoffPosition{D}) where D
    WyckoffPosition{D}(multiplicity(wp), wp.letter, compose(op, parent(wp)))
end
(*)(op::SymOperation{D}, wp::WyckoffPosition{D}) where D = compose(op, wp)


"""
$(TYPEDSIGNATURES)

Return the site symmetry group `g::SiteGroup` for a Wyckoff position `wp` in space group
`sg` (or with space group number `sgnum`; in this case, the dimensionality is inferred from
`wp`).

`g` is a group of operations that are isomorphic to the those listed in `sg` (in the sense
that they might differ by lattice vectors) and that leave the Wyckoff position `wp`
invariant, such that `all(op -> wp == op*wp, g) == true`.

The returned `SiteGroup` also contains the coset representatives of the Wyckoff position
(that are again isomorphic to those featured in `sg`), accessible via [`cosets`](@ref),
which e.g. generate the orbit of the Wyckoff position (see [`orbit(::SiteGroup)`](@ref)) and
define a left-coset decomposition of `sg` jointly with the elements in `g`.

## Example
```jldoctest sitegroup
julia> sgnum = 16;

julia> D = 2;

julia> wp = wyckoffs(sgnum, D)[3] # pick a Wyckoff position
2b: [1/3, 2/3]

julia> sg = spacegroup(sgnum, D);

julia> g  = sitegroup(sg, wp)
SiteGroup{2} ⋕16 (p6) at 2b = [1/3, 2/3] with 3 operations:
 1
 {3⁺|1,1}
 {3⁻|0,1}
```

The group structure of a `SiteGroup` can be inspected with `MultTable`:
```jldoctest sitegroup
julia> MultTable(g)
3×3 MultTable{SymOperation{2}}:
──────────┬──────────────────────────────
          │        1  {3⁺|1,1}  {3⁻|0,1} 
──────────┼──────────────────────────────
        1 │        1  {3⁺|1,1}  {3⁻|0,1} 
 {3⁺|1,1} │ {3⁺|1,1}  {3⁻|0,1}         1
 {3⁻|0,1} │ {3⁻|0,1}         1  {3⁺|1,1}
──────────┴──────────────────────────────
```

The original space group can be reconstructed from a left-coset decomposition, using the
operations and cosets contained in a `SiteGroup`:
```jldoctest sitegroup
julia> ops = [opʰ*opᵍ for opʰ in cosets(g) for opᵍ in g];

julia> Set(sg) == Set(ops)
true
```

## Terminology

Mathematically, the site symmetry group is a *stabilizer group* for a Wyckoff position,
in the same sense that the little group of **k** is a stabilizer group for a **k**-point.
"""
function sitegroup(sg::SpaceGroup{D}, wp::WyckoffPosition{D}) where D
    Nsg  = order(sg)
    Ncoset = multiplicity(wp)
    Nsite, check = divrem(Nsg, Ncoset)
    if check != 0
        throw(DomainError((Ncoset, Nsg), "Wyckoff multiplicity must divide space group order"))
    end

    siteops  = Vector{SymOperation{D}}(undef, Nsite)
    cosets   = Vector{SymOperation{D}}(undef, Ncoset)
    orbitrvs = Vector{RVec{D}}(undef, Ncoset)
    rv       = parent(wp)
    
    # both cosets and site symmetry group contains the identity operation, and the orbit 
    # automatically contains rv; add them outside loop
    siteops[1]  = cosets[1] = one(SymOperation{D})
    orbitrvs[1] = rv

    isite = icoset = 1
    for op in sg
        icoset == Ncoset && isite == Nsite && break # stop if all representatives found
        isone(op) && continue # already added identity outside loop; avoid adding twice

        W, w = unpack(op) # rotation and translation
        rv′  = op * rv
        Δ    = rv′ - rv
        Δcnst, Δfree = parts(Δ)

        # Check whether difference between rv and rv′ is a lattice vector: if so, `op` is 
        # isomorphic to a site symmetry operation; if not, to a coset operation.
        # We check this in the original lattice basis, i.e. do not force conversion to a
        # primitive basis. This is consistent with e.g. Bilbao and makes good sense.
        # The caller is of course free to do this themselves (via their choice of basis for
        # the specified `sg` and `wp`).
        if ( # tolerance'd equiv. of `all(isinteger, Δcnst) && all(iszero, Δfree)`
             all(x->isapprox(x, round(x), atol=DEFAULT_ATOL), Δcnst) && 
             all(x->abs(x)≤(DEFAULT_ATOL), Δfree) )             # ⇒ site symmetry operation

            w′ = w - Δcnst
            siteops[isite+=1] = SymOperation{D}(W, w′)

        else                                                    # ⇒ coset operation
            icoset == Ncoset && continue # we only need `Ncoset` representatives in total

            # reduce generated Wyckoff representative to coordinate range q′ᵢ∈[0,1)
            rv′′ = RVec(reduce_translation_to_unitrange(constant(rv′)), free(rv′))
            if any(rv->isapprox(rv, rv′′, nothing, false), (@view orbitrvs[OneTo(icoset)]))
                # ⇒ already included a coset op that maps to this rv′′; don't include twice
                continue
            end
            
            # shift operation so generated `WyckoffPosition`s has coordinates q′ᵢ∈[0,1)
            w′ = w - (constant(rv′) - constant(rv′′))

            # add coset operation and new Wyckoff position to orbit
            cosets[icoset+=1] = SymOperation{D}(W, w′)
            orbitrvs[icoset]  = rv′′
            
            # TODO: For certain Wyckoff positions in space groups 151:154 (8 in total), the
            #       above calculation of w′ lead to very small (~1e-16) negative values for
            #       the new representatives generated from these coset operations, i.e.
            #       `orbit(g)` can contain Wyckoff positions with negative coordinates,
            #       differing from zero by ~1e-16.
        end        
    end
    return SiteGroup{D}(num(sg), wp, siteops, cosets)
end
function sitegroup(sgnum::Integer, wp::WyckoffPosition{D}) where D
    sg = spacegroup(sgnum, Val(D))
    return sitegroup(sg, wp)
end

# `MulTable`s of `SiteGroup`s should be calculated with `modτ = false` always
function MultTable(g::SiteGroup)
    MultTable(operations(g); modτ=false)
end

"""
    orbit(g::SiteGroup)  -->  Vector{WyckoffPosition}

Compute the orbit of the Wyckoff position associated with the site symmetry group `g`.

## Extended help

The orbit of a Wyckoff position ``\\mathbf{r}`` in a space group ``G`` is defined as the
set of inequivalent points in the unit cell that can be obtained by applying the elements of
``G`` to ``\\mathbf{r}``.
Equivalently, every element of the orbit of ``\\mathbf{r}`` can be written as the
composition of a coset representative of the Wyckoff position's site group in ``G`` with
``\\mathbf{r}``.
"""
function orbit(g::SiteGroup)
    rv′s = cosets(g) .* Ref(position(g))
end

"""
    findmaximal(sitegs::AbstractVector{<:SiteGroup})

Given a vector of `SiteGroup`s associated with the Wyckoff positions of a space group,
return those `SiteGroup`s that are associated with a maximal Wyckoff positions.

Results are returned as a `view` into the input vector (i.e. as an 
`AbstractVector{<:SiteGroup}`). The associated Wyckoff positions can be retrieved via
[`position`](@ref).

## Definition
A Wyckoff position is maximal if its site symmetry group has higher order than the site
symmetry groups of any "nearby" Wyckoff positions (i.e. Wyckoff positions that can be 
connected, i.e. made equivalent, through parameter variation to the considered Wyckoff
position).

## Example
```jldoctest
julia> sgnum = 5;

julia> D = 2;

julia> wps = wyckoffs(sgnum, Val(D));

julia> sg  = spacegroup(sgnum, Val(D));


julia> sitegs = sitegroup.(Ref(sg), wps)
2-element Vector{SiteGroup{2}}:
 [1] (4b: [α, β])
 [1, m₁₀] (2a: [0, β])

julia> only(findmaximal(sitegs))
SiteGroup{2} ⋕5 (c1m1) at 2a = [0, β] with 2 operations:
 1
 m₁₀
```
"""
function findmaximal(sitegs::AbstractVector{SiteGroup{D}}) where D
    maximal = Int[]
    for (idx, g) in enumerate(sitegs)
        wp = position(g)
        v  = parent(wp)

        # if `wp` is "special" (i.e. has no "free" parameters), then it must
        # be a maximal wyckoff position, since if there are other lines
        # passing through it, they must have lower symmetry (otherwise, only
        # the line would have been included in the listings of wyckoff positions)
        isspecial(v) && (push!(maximal, idx); continue)

        # if `wp` has "free" parameters, it can still be a maximal wyckoff 
        # position, provided that it doesn't go through any other special 
        # points or intersect any other wyckoff lines/planes of higher symmetry
        has_higher_sym_nearby = false
        for (idx′, g′) in enumerate(sitegs)
            idx′ == idx && continue # don't check against self
            order(g′) < order(g) && continue # `g′` must have higher symmetry to "supersede" `g`
            
            wp′_orbit = orbit(g′) # must check for all orbits of wp′ in general
            for wp′′ in wp′_orbit
                v′ = parent(wp′′)
                if can_intersect(v, v′).bool # `wp′` can "intersect" `wp` & is higher order
                    has_higher_sym_nearby = true
                    break
                end
            end
            has_higher_sym_nearby && break
        end

        # check if there were "nearby" points that had higher symmetry, if not, `wp` is a
        # maximal wyckoff position
        has_higher_sym_nearby || push!(maximal, idx)
    end
    return @view sitegs[maximal]
end