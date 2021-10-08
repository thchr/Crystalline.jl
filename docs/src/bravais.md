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
directbasis
reciprocalbasis
primitivebasismatrix
transform(::DirectBasis, ::AbstractMatrix{<:Real})
transform(::ReciprocalBasis, ::AbstractMatrix{<:Real})
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
transform(::KVec, ::AbstractMatrix{<:Real})
primitivize(::Crystalline.AbstractVec, ::Char)
conventionalize(::Crystalline.AbstractVec, ::Char)
```