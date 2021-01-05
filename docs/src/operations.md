```@meta
Author = "Thomas Christensen"
```

# `SymOperation{D}`

A `[SymOperation{D}](@ref)` is a representation of a spatial symmetry operation $g=\{W|w\}$, composed of a rotational $W$ and a translation part $w$.
The rotational and translation parts are assumed to share the same basis system; by default, operations returned by tools in Crystalline.jl will return operations in the conventional setting of International Tables of Crystallography, Volume A (ITA).

`SymOperation`s can be constructed in two ways, either by explicitly specifying the $W$ and $w$:

```@example operations
using Crystalline, StaticArrays
W, w = (@SMatrix [1 0 0; 0 0 1; 0 1 0]), (@SVector [0, 0.5, 0])
op = SymOperation(W, w)
```
or by its equivalent triplet form
```julia
op = SymOperation{3}("x,z+1/2,y")
```
There is also a string macro accessor `@S_str` that allows triplet input via `S"x,z+1/2,y"`.

In the above output, three equivalent notations for the symmetry operation are given: first, the Seitz notation {m₀₋₁₁|0,½,0}, then the triplet notation (x,z+1/2,y), and finally the explicit matrix notation.

## Components
The components of a `SymOperation` can be accessed - although it should rarely be needed:
```@example operations
rotation(op)
```
```@example operations
translation(op)
```
```@example operations
matrix(op)
```

## Operator composition
Composition of two operators $g_1$ and $g_2$ is defined by 
$$
g_1 ∘ g_2 = \{W_1|w_1\}∘\{W_2|w_2\} = \{W_1W_2|w_1 + W_1w_2\}
$$
We can compose two `SymOperation`s in Crystalline via:
```@example operations
op1 = S"z,x,y" # 3₁₁₁⁺
op2 = S"z,y,x" # m₋₁₀₁
op3 = op1 ∘ op2
```
Note that composition is taken modulo integer lattice translations by default, such that
```@example operations
op2′ = S"z,y,x+1" # {m₋₁₀₁|001}
op1 ∘ op2
```
rather than `S"x+1,z,y"`, which is the result of direct application of the above composition rule.
To compute "unreduced" composition, the alternative `[compose](@ref)` variant of `[∘](@ref)` can be used with an optional third argument `false`:
```@example operations
compose(op1, op2′, false)
```

## Operator inverses
The operator inverse is defined as $\{W|w\} = \{W^{-1}|-W^{-1}w\}$ and can be computed via
```@example operations
inv(op1) # inv(3₁₁₁⁺)
```

## Action of symmetry operators
TODO.