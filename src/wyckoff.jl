# ---------------------------------------------------------------------------------------- #
# STRUCTS

# Wyckoff positions
struct WyckPos{D} <: AbstractVec{D}
    mult   :: Int
    letter :: Char
    qv     :: RVec{D} # associated with a single representative
end
parent(wp::WyckPos)   = wp.qv
free(wp::WyckPos)     = free(parent(wp))
constant(wp::WyckPos) = constant(parent(wp))

multiplicity(wp::WyckPos) = wp.mult
label(wp::WyckPos) = string(multiplicity(wp), wp.letter)
function transform(wp::WyckPos, P::AbstractMatrix{<:Real})
    return typeof(wp)(wp.mult, wp.letter, transform(parent(wp), P))
end

function show(io::IO, ::MIME"text/plain", wp::WyckPos)
    print(io, wp.mult, wp.letter, ": ")
    show(io, MIME"text/plain"(), parent(wp))
end

# Site symmetry groups
struct SiteGroup{D} <: AbstractGroup{D}
    num::Int
    wp::WyckPos{D}
    operations::Vector{SymOperation{D}}
    cosets::Vector{SymOperation{D}}
end
label(g::SiteGroup) = iuc(num(g), dim(g))*" at "*string(wyck(g))

"""
$(TYPEDSIGNATURES)

Return the cosets of a `SiteGroup` `g`.

The cosets generate the orbit of the Wyckoff position `wyck(g)` (see 
[`orbit(::SiteGroup, ::WyckPos)`](@ref)) and furnish a left-coset decomposition of the
underlying space group, jointly with the operations in `g` itself.
"""
cosets(g::SiteGroup) = g.cosets

"""
$(TYPEDSIGNATURES)

Return the Wyckoff position associated with a `SiteGroup`.
"""
wyck(g::SiteGroup) = g.wp

function summary(io::IO, g::SiteGroup)
    print(io, typeof(g), " #", num(g), " at ", label(wyck(g)), " = ")
    show(io, MIME"text/plain"(), parent(wyck(g)))
    print(io, " with ", length(g), " operations")
end

# ---------------------------------------------------------------------------------------- #
# CONSTRUCTORS/GETTERS FOR WYCKPOS

"""
$(TYPEDSIGNATURES)

Return the Wyckoff positions of space group `sgnum` in dimension `D` as a 
`Vector{WyckPos{D}`.

The positions are given in the conventional basis setting, following the conventions of the
Bilbao Crystallographic Server (from which the underlying data is sourced [^1]).

## Example
```jldoctest
julia> wps = get_wycks(16, 2)
4-element Vector{WyckPos{2}}:
 6d: [α, β]
 3c: [0.5, 0.0]
 2b: [0.3333333333333333, 0.6666666666666666]
 1a: [0.0, 0.0]
```

## References
[^1]: Aroyo, *et al.*,
      [Z. Kristallogr. **221**, 15-27 (2006)](https://doi.org/0.1524/zkri.2006.221.1.15)
"""
function get_wycks(sgnum::Integer, ::Val{D}=Val(3)) where D
    strarr = open(joinpath(DATA_DIR, "wyckpos/$(D)d/"*string(sgnum)*".csv")) do io
        DelimitedFiles.readdlm(io, '|', String, '\n')
    end
    mults   = parse.(Int, @view strarr[:,1])
    letters = only.(@view strarr[:,2])
    qvs     = RVec{D}.(@view strarr[:,3])

    return WyckPos{D}.(mults, letters, qvs)
end
get_wycks(sgnum::Integer, D::Integer) = get_wycks(sgnum, Val(D))

# ---------------------------------------------------------------------------------------- #
# METHODS

"""
    compose(op::SymOperation, qv::RVec) --> RVec

Return the composition of `op` ``= \\{W|w\\}`` with a real-space vector `qv::RVec`.

The operation is taken to act directly, i.e. returns ``\\{W|w\\}```qv` ``= W```qv```+w``
rather than ``\\{W|w\\}^{-1}```qv` ``= W^{-1}```qv` ``- W^{-1}w``, which can instead be
obtained from `compose(inv(op), qv)`.

Can also be called via the multipliication operator, i.e. `op * qv = compose(op, qv)`.
"""
function compose(op::SymOperation{D}, qv::RVec{D}) where D
    cnst, free = parts(qv)
    W, w       = unpack(op)

    cnst′ = W*cnst + w
    free′ = W*free

    return RVec{D}(cnst′, free′)
end
function compose(op::SymOperation{D}, wp::WyckPos{D}) where D
    WyckPos{D}(multiplicity(wp), wp.letter, compose(op, parent(wp)))
end

(*)(op::SymOperation{D}, qv::RVec{D}) where D = compose(op, qv)
(*)(op::SymOperation{D}, wp::WyckPos{D}) where D = compose(op, wp)


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
which e.g. generate the orbit of the Wyckoff position (see
[`orbit(::SiteGroup, ::WyckPos)`](@ref)) and define a left-coset decomposition of `sg`
jointly with the elements in `g`.

## Example
```jldoctest sitegroup
julia> sgnum = 16;

julia> D = 2;

julia> wp = get_wycks(sgnum, D)[3] # pick a Wyckoff position
2b: [0.3333333333333333, 0.6666666666666666]

julia> sg = spacegroup(sgnum, D);

julia> g  = SiteGroup(sg, wp)
SiteGroup{2} #16 at 2b = [0.333333, 0.666667] with 3 operations:
 1
 {3⁺|1,1}
 {3⁻|0,1}
```

The group structure of a `SiteGroup` can be inspected with `MultTable`:
```jldoctest sitegroup
julia> MultTable(g)
3×3 MultTable{2}:
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
function SiteGroup(sg::SpaceGroup{D}, wp::WyckPos{D}) where D
    Nsg  = order(sg)
    Ncoset = multiplicity(wp)
    Nsite, check = divrem(Nsg, Ncoset)
    if check != 0
        throw(DomainError((Ncoset, Nsg), "Wyckoff multiplicity must divide space group order"))
    end

    siteops  = Vector{SymOperation{D}}(undef, Nsite)
    cosets   = Vector{SymOperation{D}}(undef, Ncoset)
    orbitqvs = Vector{RVec{D}}(undef, Ncoset)
    qv       = parent(wp)
    
    # both cosets and site symmetry group contains the identity operation, and the orbit 
    # automatically contains qv; add them outside loop
    siteops[1]  = cosets[1] = one(SymOperation{D})
    orbitqvs[1] = qv

    isite = icoset = 1
    for op in sg
        icoset == Ncoset && isite == Nsite && break # stop if all necessary representatives found
        isone(op) && continue # already added identity outside loop; avoid adding twice

        W, w = unpack(op) # rotation and translation
        qv′  = op * qv
        Δ    = qv′ - qv
        Δcnst, Δfree = parts(Δ)

        # Check whether difference between qv and qv′ is a lattice vector: if so, `op` is 
        # isomorphic to a site symmetry operation; if not, to a coset operation.
        # We check this in the original lattice basis, i.e. do not force  conversion to a
        # primitive basis. This is consistent with e.g. Bilbao and makes good sense.
        # The caller is of course free to do this themselves (via their choice of basis for
        # the specified `sg` and `wp`).
        if ( # tolerance'd equiv. of `all(isinteger, Δcnst) && all(iszero, Δfree))`
             all(x->isapprox(x, round(x), atol=DEFAULT_ATOL), Δcnst) && 
             all(x->abs(x)≤(DEFAULT_ATOL), Δfree) )             # ⇒ site symmetry operation

            w′ = w - Δcnst
            siteops[isite+=1] = SymOperation{D}(W, w′)

        else                                                    # ⇒ coset operation
            icoset == Ncoset && continue # we only need `Ncoset` representatives in total

            # reduce generated Wyckoff representative to coordinate range q′ᵢ∈[0,1)
            qv′′ = RVec(reduce_translation_to_unitrange(constant(qv′)), free(qv′))
            if any(≈(qv′′), (@view orbitqvs[OneTo(icoset)]))
                # ⇒ already included a coset op that maps to this qv′′; don't include twice
                continue
            end
            
            # shift operation so generated `WyckPos`s has coordinates q′ᵢ∈[0,1)
            w′ = w - (constant(qv′) - constant(qv′′))

            # add coset operation and new Wyckoff position to orbit
            cosets[icoset+=1] = SymOperation{D}(W, w′)
            orbitqvs[icoset]  = qv′′
            
            # TODO: For certain Wyckoff positions in space groups 151:154 (8 in total), the
            #       above calculation of w′ lead to very small (~1e-16) negative values for
            #       the new representatives generated from these coset operations, i.e.
            #       `orbit(g, wp)` can contain Wyckoff positions with negative coordinates,
            #       differing from zero by ~1e-16.
        end        
    end
    return SiteGroup{D}(num(sg), wp, siteops, cosets)
end
function SiteGroup(sgnum::Integer, wp::WyckPos{D}) where D
    sg = spacegroup(sgnum, Val(D))
    return SiteGroup(sg, wp)
end

# `MulTable`s of `SiteGroup`s should be calculated with `modτ = false` always
function MultTable(g::SiteGroup)
    MultTable(operations(g); modτ=false)
end

function orbit(g::SiteGroup{D}, wp::WyckPos{D}) where D
    qv′s = cosets(g) .* Ref(wp)  # TODO: Remove this method; superfluous input `wp`
end
function orbit(g::SiteGroup{D}) where D
    qv′s = cosets(g) .* Ref(wyck(g))
end



"""
    findmaximal(sitegs::AbstractVector{<:SiteGroup})

Given a vector of `SiteGroup`s associated with the Wyckoff positions of a space group,
return those `SiteGroup`s that are associated with a maximal Wyckoff positions.

Results are returned as a `view` into the input vector (i.e. as an 
`AbstractVector{<:SiteGroup}`). The associated Wyckoff positions can subsequently be
retrieved via [`wyck`](@ref).

## Definition
A Wyckoff position is maximal if its site symmetry group has higher order than the site
symmetry groups of any "nearby" Wyckoff positions (i.e. Wyckoff positions that can be 
connected, i.e. made equivalent, through parameter variation to the considered Wyckoff
position).

## Example
```jldoctest
julia> sgnum = 5;

julia> D = 2;

julia> wps = get_wycks(sgnum, Val(D));

julia> sg  = spacegroup(sgnum, Val(D));


julia> sitegs = SiteGroup.(Ref(sg), wps)
2-element Vector{SiteGroup{2}}:
 SiteGroup{2}[1]
 SiteGroup{2}[1, m₁₀]

julia> only(findmaximal(sitegs))
SiteGroup{2} #5 at 2a = [0.0, β] with 2 operations:
 1
 m₁₀
```
"""
function findmaximal(sitegs::AbstractVector{SiteGroup{D}}) where D
    maximal = Int[]
    for (idx, g) in enumerate(sitegs)
        wp = wyck(g)
        v  = parent(wp)
        N  = order(g)

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
                if _can_intersect(v, v′) # `wp′` can "intersect" `wp` and is higher order
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


function _can_intersect(v::AbstractVec{D}, v′::AbstractVec{D};
                        atol::Real=DEFAULT_ATOL) where D
    # check if solution exists to [A] v′ = v(αβγ) or [B] v′(αβγ′) = v(αβγ) by solving
    # a least squares problem and then checking if it is a strict solution. Details:
    #   Let v(αβγ) = v₀ + V*αβγ and v′(αβγ′) = v₀′ + V′*αβγ′
    #   [A] v₀′ = v₀ + V*αβγ             ⇔  V*αβγ = v₀′-v₀
    #   [B] v₀′ + V′*αβγ′ = v₀ + V*αβγ   ⇔  V*αβγ - V′*αβγ′ = v₀′-v₀  
    #                                    ⇔  hcat(V,-V′)*vcat(αβγ,αβγ′) = v₀′-v₀
    # these equations can always be solved in the least squares sense using the
    # pseudoinverse; we can then subsequently check if the residual of that solution is in
    # fact zero, in which can the least squares solution is a "proper" solution, signaling
    # that `v` and `v′` can intersect (at the found values of `αβγ` and `αβγ′`)
    Δcnst = constant(v′) - constant(v)
    Δfree = if isspecial(v′) # `v′` is special; `v` is not
        free(v)
    else                     # neither `v′` nor `v` are special
        hcat(free(v), -free(v′))
    end
    Δfree⁻¹ = pinv(Δfree)

    # tedious detail:
    # to be safe, we have to check for equivalence between `v` and `v′` while accounting
    # for the fact that they could differ by a lattice vector; in practice, for the wyckoff
    # listings that we have have in 3D, this seems to only make a difference in a single 
    # case (SG 130, wyckoff position 8f) - but there the distinction is actually needed
    for V in Iterators.product(ntuple(_->(0.0,-1.0,1.0), Val(D))...) # loop over nearest lattice vecs
        Δcnst_plus_V = Δcnst + SVector{D,Float64}(V)
        αβγ = Δfree⁻¹*Δcnst_plus_V   # either D-dim `αβγ` or 2D-dim `hcat(αβγ, αβγ′)`
        Δ = Δcnst_plus_V - Δfree*αβγ # residual of least squares solve
        norm(Δ) < atol && return true
    end

    return false
end