# Crystalline.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url] [![Documentation (dev)][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

Tools for crystalline symmetry implemented in the Julia language.

This package provides access e.g. to the symmetry operations of crystalline point groups, space groups, Wyckoff positions, their irreducible representations and band representations, as well as tools for their associated manipulation.

## Installation

The package is registered in the General registry and can be installed from the `pkg>` prompt (accessed by typing `]` at the Julia REPL) by executing:
```julia
pkg> add https://github.com/thchr/Crystalline.jl
```
whereafter Crystalline.jl can be loaded via
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

Additional and related acessor functionality is included; e.g. point group operations (`pointgroup`) and irreps (`get_pgirreps`), elementary band representations (`bandreps`), Wyckoff positions (`get_wycks`), physically real irreps (`realify`), transformation between conventional and primitive settings (`primitivize` and `conventionalize`), and Bravais lattice utilities and conventions.
For a full description of the public API, see the [documentation][docs-dev-url].

[ci-status-img]:   https://github.com/thchr/Crystalline.jl/workflows/CI/badge.svg
[ci-status-url]:   https://github.com/thchr/Crystalline.jl/actions
[docs-dev-img]:    https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:    https://thchr.github.io/Crystalline.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://thchr.github.io/Crystalline.jl/stable
[coverage-img]:    https://codecov.io/gh/thchr/Crystalline.jl/branch/master/graph/badge.svg
[coverage-url]:    https://codecov.io/gh/thchr/Crystalline.jl

#### Limitations
At present, the package's emphasis is on spinless systems (i.e., double groups and spinful irreps are not implemented).

#### Note
Crystalline.jl is a research package in active development: breaking changes are likely (but we will strive to follow semantic versioning).