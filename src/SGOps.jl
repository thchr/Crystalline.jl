module SGOps
# packages
using HTTP, Gumbo, LinearAlgebra, Distributions, 
      JSON2, StaticArrays, Makie, TimerOutputs, 
      DelimitedFiles, SmithNormalForm,
      Meshing, JLD2
import Base: getindex, lastindex, string, isapprox,
             length, readuntil, vec, show, 
             +, -, ∘, ==
using Compat
import LinearAlgebra: inv
import PyPlot: plot, plot3D, plt
import Statistics: quantile

# constant scalars
const DEFAULT_ATOL = 1e-12 # absolute tolerance for approximate equality
const NULL_ATOL = 1e-11    # absolute tolerance for nullspace 

# included files and exports
include("utils.jl") # useful utility methods (seldom needs exporting)

include("types.jl") # defines useful types for space group symmetry analysis
export SymOperation, Crystal,               # types
       SGIrrep, MultTable, LGIrrep, KVec,
       BandRep, BandRepSet,
       SpaceGroup, PointGroup, LittleGroup,
       # operations on ...
       matrix, xyzt, operations,            # ::SymOperation
       getindex, rotation, translation, 
       issymmorph, ==,
       num, order,                          # ::AbstractGroup
       basis, dim, norms, angles,           # ::Crystal
       kstar, klabel, characters,           # ::AbstractIrrep
       label, type,
       isspecial, translations, findirrep,  # ::SGIrrep
       israyrep, kvec, irreps,              # ::LGIrrep
       string, parts,                       # ::KVec
       vec, irreplabels, reps,              # ::BandRep & ::BandRepSet 
       isspinful

include("notation.jl")
export schoenflies, hermannmauguin, 
       iuc, centering, seitz

include("symops.jl") # symmetry operations for space and plane groups
export get_sgops, xyzt2matrix, matrix2xyzt, 
       issymmorph, littlegroup, kstar,
       multtable, isgroup, checkmulttable,
       pointgroup, 
       primitivize, conventionalize, 
       reduce_ops

include("pointgroup.jl") # symmetry operations for crystallographic point groups
export get_pgops

include("bravais.jl")
export crystal, plot, crystalsystem, 
       bravaistype, primitivebasis, 
       gen_crystal, reciprocalbasis

include("irreps_reality.jl")
export herring, realify

include("../build/parse_isotropy_ir.jl")
export parseisoir, parselittlegroupirreps, 
       littlegroupirrep, klabel

include("littlegroup_irreps.jl")
export get_lgirreps, get_littlegroups,
       get_all_lgirreps

include("lattices.jl")
export UnityFourierLattice, ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!, calcfourier,
       plot

include("bandrep.jl")
export crawlbandreps, dlm2struct, 
       bandreps, 
       matrix, classification, basisdim

include("export2mpb.jl")
export prepare_mpbcalc, prepare_mpbcalc!

end # module
