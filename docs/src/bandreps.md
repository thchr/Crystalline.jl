# Elementary band representations

Crystalline.jl provides an interface to access the elementary band representations (EBRs) hosted by the Bilbao Crystallographic Server's [BANDREP](https://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl)  program via [`bandreps`](@ref).
Please cite the original research (listed in the associated docstrings).

As an example, we can obtain the all inequivalent EBRs in space group 219 (F-43c) with:
```@example ebrs
using Crystalline

brs = bandreps(219, 3) # space group 219 (dimension 3)
```
which returns a `BandRepSet`, which itself is an `AbstractVector` of `BandRep`s. This allows us to index into `brs` easily:
```@example ebrs
brs[1] # obtain the EBR induced by Wyckoff position 8a with irrep A
```

By default, `bandreps` returns the spinless EBRs with time-reversal symmetry.
This behavior can be controlled with the keyword arguments `spinful` (default, `false`) and `timereversal` (default, `true`).
By default, only minimal paths are included in the sampling of **k**-vectors; additional paths can be obtained by setting the keyword argument `allpaths = true` (default, `false`).

The distinct topological classes identifiable from symmetry can can be calculated via [`classification`](@ref), which uses the Smith normal form's principle factors:
```@example ebrs
classification(brs)
```
Which demonstrates that the symmetry indicator group of spinless particles with time-reversal symmetry in space group 219 is trivial.

## Topology and associated bases
The [`SymmetryBases.jl`](https://github.com/thchr/SymmetryBases.jl) package provides tools to analyze topology of symmetry vectors and compute associated Hilbert bases.

## API

```@meta
CurrentModule = Crystalline
```

```@docs
bandreps
classification
nontrivial_factors
basisdim
```