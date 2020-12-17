# ---------------------------------------------------------------------------------------- #
# STRUCTS

# Wyckoff positions
struct WyckPos{D} <: AbstractVec{D}
    mult   :: Int
    letter :: Char
    qv     :: RVec{D} # associated with a single representative
end
qvec(wp::WyckPos)     = wp.qv
free(wp::WyckPos)     = free(qvec(wp))
constant(wp::WyckPos) = constant(qvec(wp))

multiplicity(wp::WyckPos) = wp.mult
label(wp::WyckPos) = string(multiplicity(wp), wp.letter)

function show(io::IO, ::MIME"text/plain", wp::WyckPos)
    print(io, wp.mult, wp.letter, ": ")
    show(io, MIME"text/plain"(), qvec(wp))
end

# Site symmetry groups
struct SiteGroup{D} <: AbstractGroup{D}
    num::Int
    wp::WyckPos{D}
    operations::Vector{SymOperation{D}}
    cosets::Vector{SymOperation{D}}
end
"""
    $(SIGNATURES)

Return the cosets of a `SiteGroup` `g`.

The cosets generate the orbit of the Wyckoff position `wyck(g)` (see 
[`orbit(::SiteGroup, ::WyckPos)`](@ref)) and furnish a left-coset decomposition of the
underlying space group, jointly with the operations in `g` itself.
"""
cosets(g::SiteGroup) = g.cosets

"""
    $(SIGNATURES)

Return the Wyckoff position associated with a `SiteGroup`.
"""
wyck(g::SiteGroup)   = g.wp

function summary(io::IO, g::SiteGroup)
    print(io, typeof(g), " #", num(g), " at ", label(wyck(g)), " = ")
    show(io, MIME"text/plain"(), qvec(wyck(g)))
    print(io, " with ", length(g), " operations")
end

# ---------------------------------------------------------------------------------------- #
# CONSTRUCTORS/GETTERS FOR WYCKPOS

"""
$(SIGNATURES)

Return the Wyckoff positions of space group `sgnum` in dimension `D` as a 
`Vector{WyckPos{D}`.

The positions are given in the conventional basis setting, following the conventions of the
Bilbao Crystallographic Server (from which the underlying data is sourced [1]).

## Example
```jldoctest
julia> wps = get_wycks(16, 2)
4-element Array{WyckPos{2},1}:
 6d: [α, β]
 3c: [0.5, 0.0]
 2b: [0.3333333333333333, 0.6666666666666666]
 1a: [0.0, 0.0]
```

## References
[1] Aroyo, et. al. Zeitschrift fuer Kristallographie (2006), 221, 1, 15-27.
"""
function get_wycks(sgnum::Integer, ::Val{D}) where D
    strarr = open((@__DIR__)*"/../data/wyckpos/$(D)d/"*string(sgnum)*".csv") do io
        DelimitedFiles.readdlm(io, '|', String, '\n')
    end
    mults   = parse.(Int, @view strarr[:,1])
    letters = only.(@view strarr[:,2])
    qvs     = RVec{D}.(@view strarr[:,3])

    return WyckPos{D}.(mults, letters, qvs)
end
get_wycks(sgnum::Integer, D) = get_wycks(sgnum, Val(D))

# ---------------------------------------------------------------------------------------- #
# METHODS

"""
    ∘(op::SymOperation, qv::RVec) --> RVec

Return the composition of `op` ``= \\{W|w\\}`` with a real-space vector `qv::RVec`.

The operation is taken to act directly, i.e. returns ``\\{W|w\\}```qv` ``= W```qv```+w``
rather than ``\\{W|w\\}^{-1}```qv` ``= W^{-1}```qv` ``- W^{-1}w``, which can instead be
obtained from `inv(op)∘qv`.
"""
function (∘)(op::SymOperation{D}, qv::RVec{D}) where D
    cnst, free = parts(qv)
    W, w       = unpack(op)

    cnst′ = W*cnst + w
    free′ = W*free

    return RVec{D}(cnst′, free′)
end
(∘)(op::SymOperation{D}, wp::WyckPos{D}) where D = WyckPos{D}(multiplicity(wp), wp.letter,
                                                              op∘qvec(wp))


"""
$(SIGNATURES)

Return the site symmetry group `g::SiteGroup` for a Wyckoff position `wp` in space group
`sg` (or with space group number `sgnum`; in this case, the dimensionality is inferred from
`wp`).

`g` contains as operations that are isomorphic to the those contained in `sg` (in the sense
that they might differ by lattice vectors) that leave the Wyckoff position `wp` invariant,
such that `all(op -> wp == op∘wp, g)` is true.

The returned `SiteGroup` also contains the coset representatives of the Wyckoff position
(that are again isomorphic to those featured in `sg`), accessible via [`cosets`](@ref),
which \\eg generate the orbit of the Wyckoff position (see
[`orbit(::SiteGroup, ::WyckPos)`](@ref)) and define a left-coset decomposition of `sg`
jointly with the elements in `g`.

## Example
```jldoctest sitegroup
julia> sgnum = 16; D = 2;

julia> wp = get_wycks(sgnum, D)[3] # pick a Wyckoff position
2b: [0.3333333333333333, 0.6666666666666666]

julia> sg = spacegroup(sgnum, D);

julia> g  = SiteGroup(sg, wp)
SiteGroup{2} #16 at 2b = [0.333333, 0.666667] with 3 operations:
 1 ────────────────────────────────── (x,y)
 {3⁺|1,1} ──────────────────── (-y+1,x-y+1)
 {3⁻|0,1} ───────────────────── (-x+y,-x+1)
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
julia> ops = [opʰ∘opᵍ for opʰ in cosets(g) for opᵍ in g];

julia> Set(sg) == Set(ops)
true
```
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
    qv       = qvec(wp)
    
    # both cosets and site symmetry group contains the identity operation, and the orbit 
    # automatically contains qv; add them outside loop
    siteops[1]  = cosets[1] = one(SymOperation{D})
    orbitqvs[1] = qv

    isite, icoset = 1,1
    for op in sg
        icoset == Ncoset && isite == Nsite && break # stop if all necessary representatives found
        isone(op) && continue # already added identity outside loop; avoid adding twice

        W, w = unpack(op) # rotation and translation
        qv′  = op∘qv
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

# `MulTable`s of `SiteGroup`s should be calculated with modτ = false always
function MultTable(g::SiteGroup; verbose::Bool=false)
    MultTable(operations(g); modτ=false, verbose=verbose)
end

function orbit(g::SiteGroup{D}, wp::WyckPos{D}) where D
    qv′s = cosets(g) .∘ Ref(wp)
end