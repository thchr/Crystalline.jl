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
primitivize(::Bravais.AbstractBasis, ::Integer)
primitivize(::DirectBasis, ::Char)
primitivize(::ReciprocalBasis, ::Char)
conventionalize(::DirectBasis, ::Char)
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