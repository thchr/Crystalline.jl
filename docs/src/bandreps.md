# Elementary band representations

```@meta
CurrentModule = Crystalline
```

Crystalline.jl provides an interface to access the elementary band representations (EBRs) hosted by the Bilbao Crystallographic Server's [BANDREP](https://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl) program via [`bandreps`](@ref). Crystalline also provides [`calc_bandreps`](@ref), which computes the band representations across all maximal[^1] Wyckoff positions directly. Usually, `calc_bandreps` is preferable, as it returns more contextual information and aligns more simply with the existing conventions of Crystalline, from which it derives.

[^1]: Note that the band representations returned by [`calc_bandreps`](@ref) need not be _elementary_; i.e., a band representation returned by `calc_bandreps` may be "composite" in the "exceptional" sense defined in the original topological quantum chemistry papers (e.g., https://arxiv.org/pdf/1709.01935). The inclusion of a non-elementary band representation into a set of elementary band representations makes no difference for the purposes of analyzing band topology or band connectivity using this set, however. I.e., the set of band representations returned by `calc_bandreps` is usually equivalent to the set returned by `bandreps` (referencing the Bilbao Crystallographic Server tables), and is otherwise a strict superset. For the sake of simplicity, we will colloquially refer to the band representations returned by `calc_bandreps` as EBRs, even though the set may technically contain non-elementary band representations.

As an example, we can obtain the all inequivalent EBRs in space group 219 (F-43c) with:
```@example ebrs
using Crystalline

brs = calc_bandreps(219, Val(3)) # space group 219 (dimension 3)
```
which returns a `Collection{NewBandRep{3}}`, whose iterants are `NewBandRep{3}`s. We can inspect any individual vector in `brs`, e.g.:
```@example ebrs
brs[10] # obtain the EBR induced by Wyckoff position 8a with irrep A
```

Currently, `calc_bandreps` can only treat spinless systems; if spinful systems are required, use `bandreps`.
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
calc_bandreps
bandreps
indicator_group
indicator_group_as_string
basisdim
```