# Groups
All groups in Crystalline are concrete instances of the abstract supertype [`AbstractGroup{D}`](@ref), referring to a group in `D` dimensions. `AbstractGroup{D}` is itself a subtype of `AbstractVector{SymOperation{D}}`.
Crystalline currently supports four group types: [`SpaceGroup`](@ref), [`LittleGroup`](@ref), [`PointGroup`](@ref), and [`SiteGroup`](@ref).

## Space groups

The one, two, and three-dimensional space groups are accessible via [`spacegroup`](@ref), which takes the space group number `sgnum` and dimensino `D` as input (ideally, the dimension is provided as a `Val{D}` for the sake of type stability) and returns a `SpaceGroup{D}` structure:
```@example spacegroup
using Crystalline

D     = 3  # dimension
sgnum = 16 # space group number (≤2 in 1D, ≤17 in 2D, ≤230 in 3D)
sg    = spacegroup(sgnum, D) # where practical, `spacegroup` should be called with a `Val{D}` dimension to ensure type stability; here we have D::Int instead for simplicity
```
By default, the returned operations are given in the conventional setting of the International Tables of Crystallography, Volume A (ITA). Conversion to a primitive basis (in the CDML setting) can be accomplished via [`primitivize`](@ref).

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
To determine whether a space group is symmorphic or not, use [`issymmorph`](@ref) taking either a `SpaceGroup`, a `LittleGroup`, or a space group identified by its number and dimensionality (in the latter case, using tabulated look-up).
To test whether a given `SymOperation` is symmorphic in a given centering setting, use [`issymmorph(::SymOperation, ::Char)`](@ref)

## Group generators
Generators of `SpaceGroup`s and `PointGroup`s are accessible via [`generators`](@ref), e.g.:
```@example spacegroup
ops = generators(sgnum, SpaceGroup{D})
```

To generate a group from a list of generators, we can use the [`generate`](@ref) method. As an example, we can verify that `ops` in fact returns symmetry operations identical to those in `sg`:
```@example spacegroup
generate(ops)
```