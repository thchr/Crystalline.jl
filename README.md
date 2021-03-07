# Crystalline.jl

[![Documentation][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

Tools for crystalline symmetry implemented in the Julia language, with a focus on photonic crystals.

This package provides access e.g. to the symmetry operations of crystalline point groups, space groups, Wyckoff positions, their irreducible representations and band representations, as well as tools for their associated manipulation.

## Installation

The package is not currently registered in the General registry, but can be installed directly from this repository's URL.
To do so, enter the `pkg>` prompt by typing `]` at the Julia REPL and type:
```julia
pkg> add https://github.com/thchr/Crystalline.jl
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

julia> lgirs = get_lgirreps(16, Val(3)); # load a dictionary of small irreps and their little groups for space group #16, indexed by their k-point labels

julia> lgirs["A"] # inspect the small irreps at the A point
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

julia> CharacterTable(lgirs["Γ"]) # compute the character table for the small irreps at the Γ point
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