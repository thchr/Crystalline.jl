# Bravais types and basis systems

```@meta
CurrentModule = Crystalline
```

## Bravais.jl

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

## Extensions of Bravais.jl methods
```@docs
transform
primitivize(::Crystalline.AbstractVec, ::Char)
conventionalize(::Crystalline.AbstractVec, ::Char)
```