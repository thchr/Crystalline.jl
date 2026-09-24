# Elementary band representations

```@meta
CurrentModule = Crystalline
```

Crystalline.jl computes the elementary band representations (EBRs) across all maximal[^1] Wyckoff positions directly, via [`bandreps`](@ref). The results are checked against the EBRs tabulated by the Bilbao Crystallographic Server's [BANDREP](https://cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl) program in Crystalline's test suite, but those tables are not themselves part of the package.

[^1]: Note that the band representations returned by [`bandreps`](@ref) need not be _elementary_; i.e., a band representation returned by `bandreps` may be "composite" in the "exceptional" sense defined in the original topological quantum chemistry papers (e.g., https://arxiv.org/pdf/1709.01935). The inclusion of a non-elementary band representation into a set of elementary band representations makes no difference for the purposes of analyzing band topology or band connectivity using this set, however. I.e., the set of band representations returned by `bandreps` is usually equivalent to the set tabulated by the Bilbao Crystallographic Server, and is otherwise a strict superset. For the sake of simplicity, we will colloquially refer to the band representations returned by `bandreps` as EBRs, even though the set may technically contain non-elementary band representations.

As an example, we can obtain the all inequivalent EBRs in space group 219 (F-43c) with:
```@example ebrs
using Crystalline

brs = bandreps(219, Val(3)) # space group 219 (dimension 3)
```
which returns a `Collection{BandRep{3}}`, whose iterants are `BandRep{3}`s. We can inspect any individual vector in `brs`, e.g.:
```@example ebrs
brs[10] # obtain the EBR induced by Wyckoff position 8a with irrep A
```

By default, `bandreps` treats spinless systems; band representations of spinful systems are obtained with the keyword argument `spinful = Val(true)` (see [Double groups & spinful irreps](doublegroups.md)).
The presence or absence of time-reversal symmetry can be controlled with the keyword arguments `timereversal` (default, `true`).
By default, only maximal **k**-points are included in the projection onto little group irreps; additional **k**-points (e.g., high-symmetry lines and planes) can be obtained by setting the keyword argument `allpaths = true` (default, `false`).

A set of EBRs can be used as the basis for several analyses. For instance, we can use the EBRs to compute the symmetry indicator group, summarizing the distinct topological classes identifiable from symmetry. Crystalline.jl implements [`indicator_group`](@ref) and [`indicator_group_as_string`](@ref), which uses the Smith normal form's elementary factors to this end:
```@example ebrs
indicator_group_as_string(brs)
```
Which demonstrates that the symmetry indicator group of spinless particles with time-reversal symmetry in space group 219 is trivial.

```@meta; canonical=false
indicator_group
indicator_group_as_string
```

## Topological analysis

An EBR basis can also be used to analyze [`SymmetryVector`](@ref)s, including their topology and whether they fulfil compatibility relations.
```@meta; canonical=false
iscompatible
calc_topology
symmetry_indicators
```

## Associated bases
The [`SymmetryBases.jl`](https://github.com/thchr/SymmetryBases.jl) package provides additional tools to analyze fragile topology and to compute associated Hilbert bases.

## API

```@docs; canonical=false
bandreps
indicator_group
indicator_group_as_string
basisdim
```