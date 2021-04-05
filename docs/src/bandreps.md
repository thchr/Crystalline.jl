# Elementary band representations

Crystalline.jl provides an interface to access the elementary band representations (EBRs) hosted by the [Bilbao Crystallographic Server's BANDREP program](https://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl) via [`bandreps`](@ref).
Please cite the original research (listed in the associated docstrings).

As an example, we can obtain the all inequivalent EBRs in space group 230 (Ia-3d) with:
```@example ebrs
using Crystalline

brs = bandreps(230, 3) # elementary band representations in space group 230 (dimension 3)
```
which returns a `BandRepSet`, which itself is an `AbstractVector` of `BandRep`s. This allows us to index into `brs` easily:
```@example ebrs
brs[2] # obtain the EBR induced by Wyckoff position 16a with irrep Aᵤ
```

By default, `bandreps` returns the spinless EBRs with time-reversal symmetry.
This behavior can be controlled with the keyword arguments `spinful` (default, `false`) and `timereversal` (default, `true`).
By default, only minimal paths are included in the sampling of **k**-vectors; additional paths can be obtained by setting the keyword argument `allpaths = true` (default, `false`).

The distinct topological classes identifiable from symmetry can can be calculated via [`classification`](@ref), which uses the Smith normal form's principle factors:
```@example ebrs
classification(brs)
```

Occasionally, it is helpful to extract a matrix representation of a `BandRepSet`: this can be achieved by [`matrix(::BandRepSet)`](@ref)

The [`SymmetryBases.jl`](https://github.com/thchr/SymmetryBases.jl) package provides tools to analyze topology of symmetry vectors and compute associated Hilbert bases.

## API
```@autodocs
Modules = [Crystalline]
Pages   = ["bandreps.jl"]
Order   = [:function]
Public  = true
Private = false
```