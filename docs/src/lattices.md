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

## Symmetry-safe truncations

The Fourier resolution of the expansion is set by the `idxmax` argument of [`levelsetlattice`](@ref). This choice is not merely a matter of resolution: if `idxmax` is too low, the space group's symmetry constraints may extinguish so many orbits that too few symmetry-inequivalent orbits remain — and then *every* modulation, however generic, produces a structure whose symmetry is strictly higher than the requested space group, belonging instead to a supergroup or to a lattice with a smaller effective unit cell.

Accordingly, the default `idxmax = :auto` selects the smallest isotropic truncation that is free of such accidental symmetry (never going below the previous default of 2), via [`default_idxmax`](@ref):

```@repl levelsetlattice
default_idxmax(230, Val(3))    # the double gyroid above
default_idxmax(228, Val(3))    # Fd-3c needs a much higher truncation
```

Only 14 (cubic) space groups require a higher truncation than 2: 196, 197, 201–204, 207, 209–211, 219, 222, 226, and 228.

A truncation supplied explicitly is checked against the same criterion, and warns if it is unsafe. The check is available directly as [`is_symmetry_safe`](@ref), and — since safety is not simply a matter of the *largest* index — it accepts anisotropic truncations that an isotropic threshold would reject:

```@repl levelsetlattice
is_symmetry_safe(27, (1, 1, 1))   # generically gains an inversion center
is_symmetry_safe(27, (2, 1, 1))   # safe, despite a lower total resolution
```

To query the minimal symmetry-safe truncations of a group directly, use [`minimal_idxmax`](@ref), whose `isotropic` field gives the smallest truncation with equal components and whose `general` field gives every truncation that is symmetry-safe but ceases to be so if any single component is lowered:

```@repl levelsetlattice
m = minimal_idxmax(228, Val(3))
m.isotropic
m.general
```

The elements of `general` are mutually incomparable — none is componentwise smaller than another — so there is usually no single smallest symmetry-safe truncation. Anisotropic choices can be considerably cheaper: for space group 228, `(1,3,5)` retains 4 free coefficients where the isotropic `(5,5,5)` retains 9, at the cost of a directionally biased resolution.

## API

```@meta
CurrentModule = Crystalline
```

```@docs; canonical=false
UnityFourierLattice
ModulatedFourierLattice
levelsetlattice
default_idxmax
minimal_idxmax
is_symmetry_safe
modulate
filling2isoval
isoval2filling
primitivize(::AbstractFourierLattice, ::Char)
conventionalize(::AbstractFourierLattice, ::Char)
AbstractFourierLattice(::Any)
```