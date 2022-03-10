# Crystalline.jl

[![Documentation (stable)][docs-stable-img]][docs-stable-url] [![Documentation (dev)][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]

Tools for crystalline symmetry implemented in the Julia language.

This package provides access e.g. to the symmetry operations of crystalline point groups, space groups, Wyckoff positions, their irreducible representations and band representations, as well as tools for their associated manipulation.

## Installation

The package is registered in the General registry and can be installed via Julia's package manager from the `pkg>` prompt (accessed by typing `]` at the Julia REPL):
```julia
pkg> add Crystalline
```
whereafter Crystalline.jl can be loaded via
```julia
julia> using Crystalline
```

## Functionality

Crystalline.jl currently provides several functionalities for line groups, plane groups, and space groups, as well as crystallographic point groups.

Example usage includes:
```julia
# construct a 3D `SymOperation` from its triplet form
julia> S"x,-y,-z"
2₁₀₀ ─────────────────────────── (x,-y,-z)
 ┌ 1  0  0 ╷ 0 ┐
 │ 0 -1  0 ┆ 0 │
 └ 0  0 -1 ╵ 0 ┘

# load the `SymOperation`s of the 3D space group ⋕16 in a conventional setting
julia> sg = spacegroup(16, Val(3))
SpaceGroup{3} ⋕16 (P222) with 4 operations:
 1
 2₀₀₁
 2₀₁₀
 2₁₀₀

# load a dictionary of small irreps and their little groups for space group ⋕16,
# indexed by their k-point labels; then inspect the small irreps at the A point
julia> lgirs = lgirreps(16, Val(3))
julia> lgirs["A"]
LGIrrep{3}: ⋕16 (P222) at A = [α, 0.0, 0.5]
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

# construct the character table for the small irreps at the Γ point
julia> characters(lgirs["Γ"])
CharacterTable{3}: ⋕16 (P222 at Γ = [0.0, 0.0, 0.0])
──────┬────────────────
      │ Γ₁  Γ₂  Γ₃  Γ₄ 
──────┼────────────────
    1 │  1   1   1   1 
 2₁₀₀ │  1  -1   1  -1
 2₀₁₀ │  1  -1  -1   1
 2₀₀₁ │  1   1  -1  -1
──────┴────────────────
```

Additional functionality includes e.g. point group operations (`pointgroup`) and irreps (`pgirreps`), elementary band representations (`bandreps`), Wyckoff positions (`wyckoffs`), conjugacy classes (`classes`), class-specific characters (`classcharacters`), group generators (`generators`), physically real irreps (`realify`), and Bravais lattice utilities and conventions (accessible via the lightweight stand-alone sub-package [Bravais.jl](https://github.com/thchr/Crystalline.jl/tree/master/Bravais)).
For a full description of the public API, see the [documentation][docs-dev-url].

### Current limitations
At present, the package's emphasis is on spinless systems (i.e., double groups and spinful irreps are not implemented).

## API stability
Crystalline.jl is a research package in active development: breaking changes are likely (but will respect semantic versioning).

## Citation

If you find this package useful in your reseach, please cite our arXiv paper:

- T. Christensen, H.C. Po, J.D. Joannopoulos, & M. Soljačić, *Location and topology of the fundamental gap in photonic crystals*, [arXiv:2106.10267 (2021)](https://arxiv.org/abs/2106.10267).

In addition, please cite any works explicitly referenced in documentation for individual methods that you use.


[ci-status-img]:   https://github.com/thchr/Crystalline.jl/workflows/CI/badge.svg
[ci-status-url]:   https://github.com/thchr/Crystalline.jl/actions
[docs-dev-img]:    https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:    https://thchr.github.io/Crystalline.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://thchr.github.io/Crystalline.jl/stable
[coverage-img]:    https://codecov.io/gh/thchr/Crystalline.jl/branch/master/graph/badge.svg
[coverage-url]:    https://codecov.io/gh/thchr/Crystalline.jl