module Bravais

using LinearAlgebra
using StaticArrays

# ---------------------------------------------------------------------------------------- #

using Base:
    @propagate_inbounds,
    ImmutableDict

import Base:
    parent,
    convert,
    getindex,
    size,
    IndexStyle

# ---------------------------------------------------------------------------------------- #

export crystal,
    crystalsystem,
    bravaistype,
    centering,
    primitivebasismatrix,
    directbasis,
    reciprocalbasis,
    AbstractBasis,
    volume,
    DirectBasis,
    ReciprocalBasis,
    AbstractPoint,
    DirectPoint,
    ReciprocalPoint,
    transform,
    primitivize,
    conventionalize,
    cartesianize

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("types.jl")
include("systems.jl")
include("transform.jl")

end # module