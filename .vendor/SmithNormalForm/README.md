# Smith Normal Form

[![Build Status](https://travis-ci.org/wildart/SmithNormalForm.jl.svg?branch=master)](https://travis-ci.org/wildart/SmithNormalForm.jl)
[![Coverage Status](https://coveralls.io/repos/wildart/SmithNormalForm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wildart/SmithNormalForm.jl?branch=master)


The [Smith normal form](https://en.wikipedia.org/wiki/Smith_normal_form) decomposition over integer domain implementation in Julia.

## Installation

For Julia 1.1+, add [BoffinStuff](https://github.com/wildart/BoffinStuff.git) registry in the package manager, and proceed with the installation:

```
pkg> registry add https://github.com/wildart/BoffinStuff.git
pkg> add SmithNormalForm
```

## Example

```julia
julia> using SmithNormalForm, LinearAlgebra

julia> M = [2 4 4; -6 6 12; 10 -4 -16]
3×3 Array{Int64,2}:
  2   4    4
 -6   6   12
 10  -4  -16

julia> F = smith(M)
Smith normal form:
[2 0 0; 0 6 0; 0 0 12]

julia> F.S
3×3 Array{Int64,2}:
  1   0  0
 -3   1  0
  5  -2  1

julia> F.T
3×3 Array{Int64,2}:
 1  2  2
 0  3  4
 0  1  1

julia> diagm(F)
3×3 Array{Int64,2}:
 2  0   0
 0  6   0
 0  0  12

julia> F.S*diagm(F)*F.T
3×3 Array{Int64,2}:
  2   4    4
 -6   6   12
 10  -4  -16
```
