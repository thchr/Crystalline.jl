#__precompile__() # TODO: enable if we ever want to officially ship this - not worthwhile
                  #       while we're actively developing (increases time of "fresh" compiles)

module SGOps

# dependencies
using LinearAlgebra
using StaticArrays
using JSON2
using DelimitedFiles
using JLD2
using SmithNormalForm
using PrettyTables
using LightGraphs, MetaGraphs # → compatibility.jl
using Requires
using Compat
using Statistics: quantile

import Base: getindex, lastindex, firstindex, setindex!, # → indexing interface
             IndexStyle, size, eltype, length,           # ⤶
             string, isapprox, zero,
             readuntil, vec, show,
             +, -, ∘, ==, ImmutableDict
import LinearAlgebra: inv
import Random                 # → _Uniform in src/utils.jl
import Random: rand           # ⤶

# included files and exports
include("constants.jl")
export MAX_SGNUM

include("utils.jl") # useful utility methods (seldom needs exporting)
export get_kvpath

include("types.jl") # defines useful types for space group symmetry analysis
export SymOperation,                        # types
       DirectBasis, ReciprocalBasis,
       MultTable, LGIrrep, PGIrrep,
       KVec,
       BandRep, BandRepSet,
       SpaceGroup, PointGroup, LittleGroup,
       CharacterTable,
       # operations on ...
       matrix, xyzt,                        # ::SymOperation
       getindex, rotation, translation, 
       issymmorph, ==,
       num, order, operations,              # ::AbstractGroup
       norms, angles,                       # ::Basis
       kstar, klabel, characters,           # ::AbstractIrrep
       label, type, group,
       israyrep, kvec, irreps,              # ::LGIrrep
       isspecial, translations,
       find_lgirreps,
       dim, string, parts,                  # ::KVec
       vec, irreplabels, reps,              # ::BandRep & ::BandRepSet 
       isspinful

include("show.jl") # custom printing for structs defined in src/types.jl

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

include("symeigs2irrep.jl") # find irrep multiplicities from symmetry eigenvalue data
export find_representation

include("pointgroup.jl") # symmetry operations for crystallographic point groups
export pointgroup, get_pgirreps,
       PGS_IUCs

include("bravais.jl")
export crystal, crystalsystem,
       bravaistype,
       directbasis, reciprocalbasis

include("irreps_reality.jl")
export herring, realify

include("../build/parse_isotropy_ir.jl")
export parseisoir, parselittlegroupirreps, 
       littlegroupirrep

# Large parts of the functionality in special_representation_domain_kpoints.jl should not be
# in the core module, but belongs in a build file or similar. For now, the main goal of the
# file hasn't been achieved and the other methods are non-essential. So, we skip it.
#= 
# TODO: The `const ΦNOTΩ_KVECS_AND_MAPS = _ΦnotΩ_kvecs_and_maps_imdict()` call takes 15 s
#       precompile. It is a fundamentally awful idea to do it this way.
using CSV                     # → special_representation_domain_kpoints.jl
include("special_representation_domain_kpoints.jl")
export ΦnotΩ_kvecs
=#

include("littlegroup_irreps.jl")
export get_lgirreps, get_littlegroups,
       get_all_lgirreps

include("lattices.jl")
export UnityFourierLattice, ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!, calcfourier

include("compatibility.jl")
export subduction_count, compatibility

include("bandrep.jl")
export bandreps, matrix, classification, basisdim

include("export2mpb.jl")
export prepare_mpbcalc, prepare_mpbcalc!

# Optional code-loading, using Requires.
function __init__()
    
    # Plotting utitilities when PyPlot is loaded (also loads Meshing.jl)
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin  
        include("compat/pyplot.jl") # loads PyPlot and Meshing
        export plot, 
               plot_lattice_from_mpbparams, 
               mesh_3d_levelsetlattice
    end

end

end # module
