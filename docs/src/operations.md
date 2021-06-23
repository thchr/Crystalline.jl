```@meta
Author = "Thomas Christensen"
```

# Symmetry operations

A [`SymOperation{D}`](@ref) is a representation of a spatial symmetry operation $g=\{\mathbf{W}|\mathbf{w}\}$, composed of a rotational $\mathbf{W}$ and a translation part $\mathbf{w}$.
The rotational and translation parts are assumed to share the same basis setting; by default, operations returned by Crystalline.jl are in the conventional setting of the International Tables of Crystallography, Volume A (ITA).

`SymOperation`s can be constructed in two ways, either by explicitly specifying the $\mathbf{W}$ and $\mathbf{w}$:

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
The rotation and translation parts $\mathbf{W}$ and $\mathbf{w}$ of a `SymOperation{D}` $\{\mathbf{W}|\mathbf{w}\}$ can be accessed via [`rotation`](@ref) and [`translation`](@ref),  returning an `SMatrix{D, D, Float64}` and an `SVector{D, Float64}`, respectively.
The "augmented" matrix $[\mathbf{W}|\mathbf{w}]$ can similarly be obtained via [`matrix`](@ref).

## Operator composition
Composition of two operators $g_1$ and $g_2$ is defined by 
```math
g_1 \circ g_2 = \{\mathbf{W}_1|\mathbf{w}_1\} \circ \{\mathbf{W}_2|\mathbf{w}_2\} = \{\mathbf{W}_1\mathbf{W}_2|\mathbf{w}_1 + \mathbf{W}_1\mathbf{w}_2\}
```
We can compose two `SymOperation`s in Crystalline via:
```@example operations
op1 = S"z,x,y" # 3₁₁₁⁺
op2 = S"z,y,x" # m₋₁₀₁
op1 * op2
```
which is accessed by an overloaded call to `Base.*`, i.e. the multiplication operator (this enables us to also call derived methods of `*`, such as integer powers (e.g., `S"-y,x-y,z"^3 == S"x,y,z"`).
Note that composition is taken modulo integer lattice translations by default, such that
```@example operations
op2′ = S"z,y,x+1" # {m₋₁₀₁|001}
op1 * op2′ # equivalent to compose(op1, op2′, true)
```
rather than `S"x+1,z,y"`, which is the result of direct application of the above composition rule.
To compute "unreduced" composition, the more precise [`compose`](@ref) variant of [`*`](@ref) can be used with an optional third argument `false`:
```@example operations
compose(op1, op2′, false)
```

## Operator inverses
The operator inverse is defined as $\{\mathbf{W}|\mathbf{w}\} = \{\mathbf{W}^{-1}|-\mathbf{W}^{-1}\mathbf{w}\}$ and can be computed via
```@example operations
inv(op1) # inv(3₁₁₁⁺)
```

## Action of symmetry operators
A `SymOperation` can act on vectors in direct ([`RVec`](@ref)) or reciprocal ([`KVec`](@ref)) space.
When acting in reciprocal space, translation parts of a `SymOperation` have no effect.