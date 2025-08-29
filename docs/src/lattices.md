# Isosurfaces with space group symmetry

Crystalline.jl implements a function [`levelsetlattice`](@ref) to generate symmetry-constrained periodic isosurfaces, following the approach described in Supplemental Section S3.D of [Phys. Rev. X **12**, 021066 (2022)](https://doi.org/10.1103/PhysRevX.12.021066) which relates and constrains the orbits of an expansion in reciprocal-lattice plane waves.

The resulting isosurfaces can be visualized using a 3D-capable backend of Makie.jl such as GLMakie.jl. 

## Example
To illustrate the functionality, we construct and visualize an isosurface for the double gyroid in space group 230. First, we build a "base", unparameterized surface for space group 230:

```@repl levelsetlattice
using Crystalline
flat = levelsetlattice(230, Val(3))
```
By default, [`levelsetlattice`](@ref) returns a `UnityFourierLattice`, with the "joint" coefficient of each orbit set to unity. These coefficients can be freely chosen, however, with each choice of coefficients resulting in a different symmetry-respecting surface. Imposing a set of coefficients can be accomplished with [`modulate(flat, modulation)`](@ref), where `modulation` is a vector of orbit-coefficients (random, if unspecified), which multiplies onto the coefficients of each orbit.

```@repl levelsetlattice
mflat = modulate(flat, [0, 1, 0.5])
```

The particular choice of modulation `[0, 1, 0.5]` above creates a double gyroid-like structure, but many other outcomes are possible, depending on the modulation.
We can visualize the resulting structure using GLMakie.jl (for which we must also supply a set of basis vectors):

```@repl levelsetlattice
using GLMakie
Rs = directbasis(230, Val(3))
```

```@example levelsetlattice
plot(mflat, Rs; filling = 0.3)
```

In the above, `filling` sets the filling fraction of the displayed lattice (see also [`filling2isoval`](@ref) and [`isoval2filling`](@ref) to map between filling fractions and isovalues).

By default, `levelsetlattice` (and `directbasis`) operates in the conventional basis. Conversion to a primitive basis (here, body-centered with centering type `I`) can be achieved via `primitivize`:

```@example levelsetlattice
plot(primitivize(mflat, 'I'), primitivize(Rs, 'I'); filling = 0.3)
```

## API

```@meta
CurrentModule = Crystalline
```

```@docs; canonical=false
UnityFourierLattice
ModulatedFourierLattice
levelsetlattice
modulate
filling2isoval
isoval2filling
primitivize(::AbstractFourierLattice, ::Char)
conventionalize(::AbstractFourierLattice, ::Char)
AbstractFourierLattice(::Any)
```