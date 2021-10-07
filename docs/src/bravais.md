# Bravais types and basis systems

## Bravais.jl

```@meta
CurrentModule = Bravais
```

```@docs
crystal
crystalsystem
directbasis
bravaistype
reciprocalbasis
primitivebasismatrix
transform(::DirectBasis, ::AbstractMatrix)
transform(::ReciprocalBasis, ::AbstractMatrix)
primitivize(::Bravais.AbstractBasis, ::Integer)
primitivize(::DirectBasis, ::Char)
primitivize(::ReciprocalBasis, ::Char)
conventionalize(::DirectBasis, ::Char)
conventionalize(::ReciprocalBasis, ::Char)
```

## Crystalline.jl extensions of Bravais.jl functions

```@meta
CurrentModule = Crystalline
```

```@docs
transform
primitivize(::Crystalline.AbstractVec, ::Char)
conventionalize(::Crystalline.AbstractVec, ::Char)
```