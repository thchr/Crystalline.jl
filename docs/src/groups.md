# Groups
All groups in Crystalline are concrete instances of the abstract supertype [`Crystalline.AbstractGroup{D}`](@ref), referring to a group in `D` dimensions. `Crystalline.AbstractGroup{D}` is itself a subtype of `AbstractVector{SymOperation{D}}`.
Crystalline currently supports five group types: [`SpaceGroup`](@ref), [`PointGroup`](@ref), [`LittleGroup`](@ref), [`SubperiodicGroup`](@ref), [`SiteGroup`](@ref), and [`MSpaceGroup`](@ref).

## Example: space groups

The one, two, and three-dimensional space groups are accessible via [`spacegroup`](@ref), which takes the space group number `sgnum` and dimension `D` or `Val(D)` as input and returns a `SpaceGroup{D}`:
```@example spacegroup
using Crystalline

D     = 3  # dimension
sgnum = 16 # space group number (≤2 in 1D, ≤17 in 2D, ≤230 in 3D)
sg    = spacegroup(sgnum, D) # or `spacegroup(sgnum, Val(D))`
```
Where practical, `spacegroup` should be called with a `Val(D)` dimension to ensure type stability; here and elsewhere in the documentation, we will often use `D::Int` instead for simplicity.

By default, the returned operations are given in the conventional setting of the International Tables of Crystallography, Volume A (ITA). Conversion to a primitive basis (in the CDML setting) can be accomplished via [`primitivize`](@ref).

In addition to space groups, Crystalline.jl provides access to the operations of point groups ([`pointgroup`](@ref)), little groups ([`littlegroups`](@ref)), subperiodic groups ([`subperiodicgroup`](@ref); including rod, layer, and frieze groups), site symmetry groups ([`sitegroup`](@ref) and [`sitegroups`](@ref)), and magnetic space groups ([`mspacegroup`](@ref)).

### Multiplication tables
We can compute the multiplication table of a space group (under the previously defined notion of operator composition) using [`MultTable`](@ref):
```@example spacegroup
MultTable(sg)
```

Alternatively, exploiting overloading of the `*`-operator, "raw" multiplication tables can be constructed via a simple outer product:
```@example spacegroup
sg .* permutedims(sg) # equivalent to `reshape(kron(sg, sg), (length(sg), length(sg)))`
```

### Symmorphic vs. nonsymorphic space groups
To determine whether a space group is symmorphic or not, use [`issymmorph`](@ref) taking either a `SpaceGroup`, `LittleGroup`, or `SubperiodicGroup` (or a `SpaceGroup` identified by its number and dimensionality; in this case, using tabulated look-up).
To test whether a given `SymOperation` is symmorphic in a given centering setting, use [`issymmorph(::SymOperation, ::Char)`](@ref)

## Group generators
Generators of `SpaceGroup`s, `PointGroup`s, and `SubperiodicGroup`s are accessible via [`generators`](@ref), e.g.:
```@example spacegroup
ops = generators(sgnum, SpaceGroup{D})
```

To generate a group from a list of generators, we can use the [`generate`](@ref) method. As an example, we can verify that `ops` in fact returns symmetry operations identical to those in `sg`:
```@example spacegroup
generate(ops)
```

## Magnetic space groups
Magnetic space groups are accessible via [`mspacegroup`](@ref).
