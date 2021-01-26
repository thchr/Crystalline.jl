#__precompile__() # TODO: enable if we ever want to officially ship this - not worthwhile
                  #       while we're actively developing (increases time of "fresh" compiles)

module Crystalline

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
using DocStringExtensions

using Base: OneTo, @propagate_inbounds 

import Base: getindex, setindex!,      # → iteration/AbstractArray interface
             IndexStyle, size, length, # ⤶
             iterate,
             string, isapprox, zero,
             readuntil, vec, show, summary,
             +, -, ∘, ==, ImmutableDict,
             isone, one
import LinearAlgebra: inv
import Random                 # → _Uniform in src/utils.jl
import Random: rand           # ⤶

# include submodules
include("SquareStaticMatrices.jl")
using .SquareStaticMatrices # exports SSqMatrix{D,T}

# included files and exports
include("constants.jl")
include("newfile.jl")
#include("Band_Representations.jl")
export MAX_SGNUM

include("utils.jl") # useful utility methods (seldom needs exporting)
export splice_kvpath, interpolate_kvpath

include("types.jl") # defines useful types for space group symmetry analysis
export SymOperation,                        # types
       DirectBasis, ReciprocalBasis,
       Reality, REAL, PSEUDOREAL, COMPLEX,
       MultTable, LGIrrep, PGIrrep,
       KVec, RVec,
       BandRep, BandRepSet,
       SpaceGroup, PointGroup, LittleGroup,
       CharacterTable,
       # operations on ...
       matrix, xyzt,                        # ::SymOperation
       getindex, rotation, translation, 
       issymmorph,
       num, order, operations,              # ::AbstractGroup
       norms, angles,                       # ::Basis
       kstar, klabel, characters,           # ::AbstractIrrep
       label, reality, group,
       israyrep, kvec, irreps,              # ::LGIrrep
       isspecial, translations,
       dim, parts,                          # ::KVec & RVec
       vec, irreplabels, klabels, kvecs,    # ::BandRep & ::BandRepSet 
       isspinful

include("show.jl") # custom printing for structs defined in src/types.jl

include("notation.jl")
export schoenflies, hermannmauguin, iuc,
       centering, seitz, mulliken

include("symops.jl") # symmetry operations for space, plane, and line groups
export @S_str, spacegroup, compose,
       issymmorph, littlegroup, kstar,
       pointgroup,
       primitivize, conventionalize, cartesianize,
       reduce_ops, transform,
       issubgroup, isnormal,
       generate

include("wyckoff.jl") # wyckoff positions and site symmetry groups
export get_wycks, WyckPos,
       multiplicity, vec,
       SiteGroup, orbit, cosets, wyck,
       findmaximal

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
export calc_reality, realify

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
export get_lgirreps, get_littlegroups

include("lattices.jl")
export ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!, calcfourier

include("compatibility.jl")
export subduction_count

include("bandrep.jl")
export bandreps, matrix, classification, basisdim

include("export2mpb.jl")
export prepare_mpbcalc, prepare_mpbcalc!

## __init__
# - open .jld2 data files, so we don't need to keep opening/closing them
# - optional code-loading, using Requires.
const DATA_DIR = dirname(@__DIR__)*"/data"
function __init__()
    # Open LGIrrep and LittleGroup data files for read access on package load (this saves
    # us a lot of time, compared to doing `jldopen` each time we need to e.g. call 
    # `get_lgirreps`, where the time for opening/closing would otherwise dominate)
    global LGIRREPS_JLDFILES = ImmutableDict((D=>JLD2.jldopen(DATA_DIR*"/lgirreps/$(D)d/irreps_data.jld2", "r")       for D in (1,2,3))...)
    global LGS_JLDFILES      = ImmutableDict((D=>JLD2.jldopen(DATA_DIR*"/lgirreps/$(D)d/littlegroups_data.jld2", "r") for D in (1,2,3))...)
    global PGIRREPS_JLDFILE  = JLD2.jldopen(DATA_DIR*"/pgirreps/3d/irreps_data.jld2", "r") # only has 3D data; no need for Dict
    # ensure we close files on exit
    atexit(() -> foreach(jldfile -> close(jldfile), values(LGIRREPS_JLDFILES)))
    atexit(() -> foreach(jldfile -> close(jldfile), values(LGS_JLDFILES)))
    atexit(() -> close(PGIRREPS_JLDFILE))

    # Plotting utitilities when PyPlot is loaded (also loads Meshing.jl)
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin  
        include("compat/pyplot.jl") # loads PyPlot and Meshing
        export plot, 
               plot_lattice_from_mpbparams, 
               mesh_3d_levelsetlattice
    end
end

end # module
