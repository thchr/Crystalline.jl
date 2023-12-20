# Bravais.jl

Bravais types, basis systems, and transformations between conventional and primitive settings.

## API

```@meta
CurrentModule = Bravais
```

### Types
```@docs
AbstractBasis
DirectBasis
ReciprocalBasis
AbstractPoint
DirectPoint
ReciprocalPoint
```

### Crystal systems & Bravais types
```@docs
crystalsystem
bravaistype
centering
```

### Basis construction
```@docs
crystal
directbasis
reciprocalbasis
```

### Transformations
```@docs
primitivebasismatrix
transform
primitivize
conventionalize
cartesianize
cartesianize!
latticize
latticize!
```

### Miscellaneous
```@docs
volume
metricmatrix
stack
```

## Crystalline.jl extensions of Bravais.jl functions

```@meta
CurrentModule = Crystalline
```

### `SymOperation`
```@docs
transform(::SymOperation, ::AbstractMatrix{<:Real}, ::Union{AbstractVector{<:Real}, Nothing}, ::Bool=true)
primitivize(::SymOperation, ::Char, ::Bool)
conventionalize(::SymOperation, ::Char, ::Bool)
```

### `AbstractVec`
```@docs
transform(::Crystalline.AbstractVec, ::AbstractMatrix{<:Real})
primitivize(::Crystalline.AbstractVec, ::Char)
conventionalize(::Crystalline.AbstractVec, ::Char)
```

### `AbstractFourierLattice`
```@docs
primitivize(::AbstractFourierLattice, ::Char)
conventionalize(::AbstractFourierLattice, ::Char)
```