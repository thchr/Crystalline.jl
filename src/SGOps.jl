module SGOps
# packages
using HTTP, Gumbo, LinearAlgebra, Distributions, JSON2, StaticArrays, Makie, TimerOutputs, Printf
import Base: getindex, lastindex, ∘, ==, string
import PyPlot: plot, plot3D, plt
import Statistics: quantile


# included files and exports
include("types.jl") # defines useful types for space group symmetry analysis
export SpaceGroup, SymOperation, Crystal, # types
       Irrep, MultTable, LGIrrep, KVec,
       # operations on ...
       matrix, xyzt, operations,          # ::SymOperation
       getindex, pg, translation, 
       issymmorph, ==,
       num, order,                        # ::SpaceGroup
       basis, dim, norms, angles,         # ::Crystal
       irreps, characters,                # ::Irrep
       label, isspecial, kstar,
       translations, findirrep,
       israyrep,
       string, parts                      # ::KVec

include("notation.jl")
export schoenflies, hermannmauguin, iuc, centering

include("symops.jl") # crawls symmetry operations from Bilbao
export get_symops, xyzt2matrix, matrix2xyzt, 
       issymmorph, littlegroup, starofk,
       multtable, isgroup, checkmulttable

include("bravais.jl")
export crystal, plot, crystalsystem, 
       bravaistype, primitivebasis, 
       gen_crystal, reciprocalbasis

include("../build/parse_isotropy_ir.jl")
export parseisoir, parselittlegroupirreps, 
       littlegroupirrep, klabel

include("lattices.jl")
export levelsetlattice, plotfourier, plotiso

#include("crawl_kvecs.jl")

end # module
