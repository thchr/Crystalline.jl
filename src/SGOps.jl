module SGOps
# packages
using HTTP, Gumbo, LinearAlgebra, Distributions, 
      JSON2, StaticArrays, Makie, TimerOutputs,
      DelimitedFiles, SmithNormalForm,
import Base: getindex, lastindex, firstindex,  setindex!, IndexStyle, size, 
             eltype, length,                                            # indexing interface
             string, isapprox,
             readuntil, vec, show, 
             +, -, ∘, ==, ImmutableDict
using Compat
import LinearAlgebra: inv
import PyPlot: plot, plot3D, plt
import Statistics: quantile



# included files and exports
include("constants.jl")
export MAX_SGNUM

include("utils.jl") # useful utility methods (seldom needs exporting)
export get_kvpath

include("types.jl") # defines useful types for space group symmetry analysis
export SymOperation,                        # types
       DirectBasis, ReciprocalBasis,
       SGIrrep, MultTable, LGIrrep, 
       PGIrrep,
       KVec,
       BandRep, BandRepSet,
       SpaceGroup, PointGroup, LittleGroup,
       CharacterTable,
       # operations on ...
       matrix, xyzt, operations,            # ::SymOperation
       getindex, rotation, translation, 
       issymmorph, ==,
       num, order,                          # ::AbstractGroup
       norms, angles,                       # ::Basis
       kstar, klabel, characters,           # ::AbstractIrrep
       label, type, group,
       isspecial, translations,             # ::SGIrrep
       israyrep, kvec, irreps,              # ::LGIrrep
       find_lgirreps,
       dim, string, parts,                  # ::KVec
       vec, irreplabels, reps,              # ::BandRep & ::BandRepSet 
       isspinful

include("notation.jl")
export schoenflies, hermannmauguin, iuc,
       centering, seitz

include("symops.jl") # symmetry operations for space, plane, and line groups
export spacegroup, xyzt2matrix, matrix2xyzt,
       ∘, compose,
       issymmorph, littlegroup, kstar,
       multtable, isgroup, checkmulttable,
       pointgroup,
       primitivize, conventionalize,
       reduce_ops, transform,
       issubgroup, isnormal

include("pointgroup.jl") # symmetry operations for crystallographic point groups
export pointgroup, get_pgirreps

include("bravais.jl")
export crystal, plot, crystalsystem,
       bravaistype,
       directbasis, reciprocalbasis

include("irreps_reality.jl")
export herring, realify

include("../build/parse_isotropy_ir.jl")
export parseisoir, parselittlegroupirreps, 
       littlegroupirrep, klabel

include("special_representation_domain_kpoints.jl")
export ΦnotΩ_kvecs

include("littlegroup_irreps.jl")
export get_lgirreps, get_littlegroups,
       get_all_lgirreps

include("lattices.jl")
export UnityFourierLattice, ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!, calcfourier,
       plot

include("compatibility.jl")
export subduction_count, compatibility

include("bandrep.jl")
export crawlbandreps, dlm2struct,
       bandreps,
       matrix, classification, basisdim

include("export2mpb.jl")
export prepare_mpbcalc, prepare_mpbcalc!

end # module
