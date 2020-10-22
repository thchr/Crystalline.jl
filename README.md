# Crystalline.jl

[![Documentation][docs-dev-img]][docs-dev-url] [![Build status][ci-status-img]][ci-status-url] [![Coverage][coverage-img]][coverage-url]
Work-in-progress research package for crystalline symmetry, irreps, and bandreps.

## Installation

The package is currently not registred in the General registry, and additionally relies on a package not listed in General ([SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl)). 
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

[ci-status-img]: https://github.com/thchr/Crystalline.jl/workflows/CI/badge.svg
[ci-status-url]: https://github.com/thchr/Crystalline.jl/actions
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://thchr.github.io/Crystalline.jl/dev
[coverage-img]: https://codecov.io/gh/thchr/Crystalline.jl/branch/master/graph/badge.svg
[coverage-url]: https://codecov.io/gh/thchr/Crystalline.jl