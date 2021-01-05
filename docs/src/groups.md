# Groups
All groups in Crystalline are concrete instances of the abstract supertype [`AbstractGroup{D}`](@ref), referring to a group in `D` dimensions. `AbstractGroup{D}` is itself a subtype of `AbstractVector{SymOperation{D}}`.
Crystalline currently supports four group types: [`SpaceGroup`](@ref), [`LittleGroup`](@ref), [`PointGroup`](@ref), and [`SiteGroup`](@ref).

## `SpaceGroup`

Access to space groups in a conventional setting is facilitated via:
```@example spacegroup
using Crystalline

D     = 3  # dimension
sgnum = 16 # space group number (≤2 in 1D, ≤17 in 2D, ≤230 in 3D)
sg    = spacegroup(sgnum, D) # where practical, `spacegroup` should be called with a `Val{D}` dimension to ensure type stability; here we have D::Int instead for simplicity
```
We can test that `sg` indeed forms a group under our previously defined composition operator, using [`MultTable`](@ref):
```@example spacegroup
MultTable(sg)
```