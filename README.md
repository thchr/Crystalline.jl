# Crystalline.jl

[![Documentation][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

Work-in-progress research package for crystalline symmetry, irreps, and bandreps.

## Installation

The package is not currently registred in the General registry, and additionally relies on a package not listed in General ([SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl)). 
As a result, installation is slightly more complicated than ordinarily (albeit, not by much).

To install, enter the `pkg>` prompt by typing `]` at the REPL. Then execute the following commands:
```julia
(@v1.5) pkg> registry add https://github.com/wildart/BoffinStuff.git # needed for the unregistred package SmithNormalForm (which Crystalline depends on)
(@v1.5) pkg> add https://github.com/thchr/Crystalline.jl
```
which will make Crystalline.jl available via 
```julia
julia> using Crystalline
```

## Quick-start

Crystalline.jl currently provides several functionalities for line groups, plane groups, and space groups, as well as crystallographic point groups.

Some example usage below:
```julia
julia> S"x,-y,-z" # construct a 3D `SymOperation` from its triplet form
2₁₀₀ ─────────────────────────── (x,-y,-z)
 ┌ 1  0  0 ╷ 0 ┐
 │ 0 -1  0 ┆ 0 │
 └ 0  0 -1 ╵ 0 ┘

julia> sg = spacegroup(16, Val(3)) # load the `SymOperation`s of the 3D space group #16 in a conventional setting
SpaceGroup{3} #16 (P222) with 4 operations:
 1 ──────────────────────────────── (x,y,z)
 2₀₀₁ ─────────────────────────── (-x,-y,z)
 2₀₁₀ ─────────────────────────── (-x,y,-z)
 2₁₀₀ ─────────────────────────── (x,-y,-z)

julia> lgirs = get_lgirreps(16, Val(3)); # load a dictionary of small irreps and their little groups for space group #16

julia> lgirs["Γ"] # inspect the small irreps at the Γ point
LGIrrep{3}: #16 (P222) at A = [α, 0.0, 0.5]
A₁ ─┬─────────────────────────────────────────────
    ├─ 1: ──────────────────────────────── (x,y,z)
    │     1.0
    │     
    ├─ 2₁₀₀: ─────────────────────────── (x,-y,-z)
    │     1.0
    └─────────────────────────────────────────────
A₂ ─┬─────────────────────────────────────────────
    ├─ 1: ──────────────────────────────── (x,y,z)
    │     1.0
    │     
    ├─ 2₁₀₀: ─────────────────────────── (x,-y,-z)
    │     -1.0
    └─────────────────────────────────────────────

julia> CharacterTable(lgirs["Γ"]) # build the character table for the Γ point small irreps
CharacterTable{3}: #16 (P222 at Γ = [0.0, 0.0, 0.0])
──────┬────────────────
      │ Γ₁  Γ₂  Γ₃  Γ₄ 
──────┼────────────────
    1 │  1   1   1   1 
 2₁₀₀ │  1  -1   1  -1
 2₀₁₀ │  1  -1  -1   1
 2₀₀₁ │  1   1  -1  -1
──────┴────────────────
```

Additional and related acessor functionality is included; e.g.  `pointgroup`, `get_pgirreps`, `bandreps`, `get_wycks`; see the [documentation][docs-dev-url] for a full description of the public API.

[ci-status-img]: https://github.com/thchr/Crystalline.jl/workflows/CI/badge.svg
[ci-status-url]: https://github.com/thchr/Crystalline.jl/actions
[docs-dev-img]:  https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:  https://thchr.github.io/Crystalline.jl/dev
[coverage-img]:  https://codecov.io/gh/thchr/Crystalline.jl/branch/master/graph/badge.svg
[coverage-url]:  https://codecov.io/gh/thchr/Crystalline.jl