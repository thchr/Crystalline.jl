# Bravais.jl

Bravais types, basis systems, and transformations between conventional and primitive settings.

## API

```@meta
CurrentModule = Bravais
```

```@docs
crystal
crystalsystem
bravaistype
centering
DirectBasis
ReciprocalBasis
DirectPoint
ReciprocalPoint
directbasis
reciprocalbasis
primitivebasismatrix
transform
primitivize
conventionalize
```

## Crystalline.jl extensions of Bravais.jl functions

```@meta
CurrentModule = Crystalline
```

```@docs
transform(::KVec, ::AbstractMatrix{<:Real})
primitivize(::Crystalline.AbstractVec, ::Char)
conventionalize(::Crystalline.AbstractVec, ::Char)
```