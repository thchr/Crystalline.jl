module SGOps
# packages
using HTTP, Gumbo, LinearAlgebra, Distributions, 
      JSON2, StaticArrays, Makie, TimerOutputs, 
      Printf, DelimitedFiles, SmithNormalForm,
      BSON
import Base: getindex, lastindex, string, isapprox,
             length, readuntil, vec, show, 
             +, -, ∘, ==
import LinearAlgebra: inv
import PyPlot: plot, plot3D, plt
import Statistics: quantile

# constant scalars
const DEFAULT_ATOL = 1e-12 # absolute tolerance for approximate equality
const NULL_ATOL = 1e-11    # absolute tolerance for nullspace 

# included files and exports
include("utils.jl") # useful utility methods (seldom needs exporting)

include("types.jl") # defines useful types for space group symmetry analysis
export SpaceGroup, SymOperation, Crystal, # types
       SGIrrep, MultTable, LGIrrep, KVec,
       BandRep, BandRepSet, 
       # operations on ...
       matrix, xyzt, operations,          # ::SymOperation
       getindex, rotation, translation, 
       issymmorph, ==,
       num, order,                        # ::SpaceGroup
       basis, dim, norms, angles,         # ::Crystal
       irreps, characters,                # ::SGIrrep
       label, isspecial, kstar,
       translations, findirrep,
       type, klabel,
       israyrep, kvec,                    # ::LGIrrep
       string, parts,                     # ::KVec
       vec, irreplabels, reps,            # ::BandRep & ::BandRepSet 
       isspinful

include("notation.jl")
export schoenflies, hermannmauguin, 
       iuc, centering, seitz

include("symops.jl") # crawls symmetry operations from Bilbao
export get_symops, xyzt2matrix, matrix2xyzt, 
       issymmorph, littlegroup, starofk,
       multtable, isgroup, checkmulttable,
       pointgroup, primitivize, conventionalize

include("bravais.jl")
export crystal, plot, crystalsystem, 
       bravaistype, primitivebasis, 
       gen_crystal, reciprocalbasis

include("../build/parse_isotropy_ir.jl")
export parseisoir, parselittlegroupirreps, 
       littlegroupirrep, klabel, herring,
       realify

include("lattices.jl")
export UnityFourierLattice, ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!, calcfourier,
       plotfourier, plotiso

include("bandrep.jl")
export crawlbandreps, dlm2struct, 
       bandreps, 
       matrix, classification, basisdim

end # module
