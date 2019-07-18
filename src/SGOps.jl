module SGOps
# packages
import Base: getindex, lastindex
using HTTP
using Gumbo
import PyPlot: plot, plot3D, plt
import LinearAlgebra: norm, dot
using Distributions
using JSON2

# included files and exports
include("types.jl") # defines useful types for space group symmetry analysis
export SpaceGroup, SymOperation, Crystal, # types
       # operations on ...
       matrix, shorthand, operations,     # ::SymOperation
       getindex, pg, translation, issymmorph,
       num, order,                        # ::SpaceGroup
       basis, dim, norms, angles          # ::Crystal

include("notation.jl")
export schoenflies, hermannmauguin, centering

include("symops.jl") # crawls symmetry operations from Bilbao
export get_symops, xyzt_op, issymmorph

include("bravais.jl")
export crystal, plot, crystalsystem, 
       bravaistype, primitivebasis, 
       gen_crystal

include("lattices.jl")
export gen_lattice,  symmetrize!, plotlattice

#include("crawl_kvecs.jl")

end # module
