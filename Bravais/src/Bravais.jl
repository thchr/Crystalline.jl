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
    DirectBasis,
    ReciprocalBasis,
    transform,
    primitivize,
    conventionalize

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("types.jl")
include("systems.jl")
include("transform.jl")

end # module