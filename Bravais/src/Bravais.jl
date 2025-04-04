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
    centering_volume_fraction,
    directbasis,
    dualbasis,
    AbstractBasis,
    volume,
    metricmatrix,
    DirectBasis,
    ReciprocalBasis,
    AbstractPoint,
    DirectPoint,
    ReciprocalPoint,
    transform,
    primitivize,
    conventionalize,
    cartesianize,
    latticize,
    nigglibasis

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("types.jl")
include("show.jl")
include("systems.jl")
include("transform.jl")
include("niggli.jl")
include("deprecations.jl")

end # module