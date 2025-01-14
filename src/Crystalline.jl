module Crystalline

# dependencies
using LinearAlgebra
using StaticArrays
using DelimitedFiles
using JLD2
using PrettyTables
using Combinatorics           # → `find_isomorphic_parent_pointgroup` in pointgroup.jl
using Requires
using Reexport
using DocStringExtensions
import Graphs

using Base: OneTo, @propagate_inbounds

import Base: getindex, setindex!,      # → iteration/AbstractArray interface
             IndexStyle, size, copy,   # ⤶
             iterate,
             string, zero,
             readuntil, show, summary,
             *, +, -, ==, ImmutableDict,
             isone, one,
             convert, parent,
             sort!
import LinearAlgebra: inv


# include submodules
include("SquareStaticMatrices.jl")
using .SquareStaticMatrices # exports `SSqMatrix{D,T}`

# include vendored SmithNormalForm.jl package from ../.vendor/
include("../.vendor/SmithNormalForm/src/SmithNormalForm.jl")
using .SmithNormalForm: smith, Smith
export smith, Smith # export, so that loading Crystalline also namespaces these

@reexport using Bravais
import Bravais: primitivize, conventionalize, cartesianize, transform, centering
using Bravais: stack, all_centeringtranslations, centeringtranslation,
               centering_volume_fraction

include("surface_unitcell.jl")
export surface_basis # TODO: move to Bravais (but tricky cf. SmithNormalForm dependency)

# included files and exports
include("constants.jl")
export MAX_SGNUM, MAX_SUBGNUM, MAX_MSGNUM, MAX_MSUBGNUM, ENANTIOMORPHIC_PAIRS

include("utils.jl") # useful utility methods (seldom needs exporting)
export splice_kvpath, interpolate_kvpath

include("types.jl") # defines useful types for space group symmetry analysis
export SymOperation,                        # types
       DirectBasis, ReciprocalBasis,
       Reality, REAL, PSEUDOREAL, COMPLEX,
       Collection,
       MultTable, LGIrrep, PGIrrep, SiteIrrep,
       KVec, RVec,
       BandRep, BandRepSet,
       SpaceGroup, PointGroup, LittleGroup,
       CharacterTable,
       # operations on ...
       matrix, xyzt,                        # ::AbstractOperation
       rotation, translation, 
       issymmorph,                          # ::SymOperation
       num, order, operations,              # ::AbstractGroup
       klabel, characters,                  # ::AbstractIrrep
       classcharacters,
       label, reality, group,
       ⊕,
       israyrep,                            # ::LGIrrep
       isspecial,
       irdim,
       dim, parts,                          # ::KVec & RVec
       irreplabels, klabels,                # ::BandRep & ::BandRepSet 
       isspinful

include("types_symmetryvectors.jl")
export SymmetryVector, SymmetryVectors, NewBandRep, CompositeBandRep
export irreps, multiplicities, occupation

include("notation.jl")
export schoenflies, iuc, centering, seitz, mulliken

include("subperiodic.jl")
export SubperiodicGroup

include("magnetic/notation-data.jl")
include("magnetic/types.jl")
export MSymOperation, MSpaceGroup

include("tables/rotation_translation.jl")
include("tables/groups/pointgroup.jl")
include("tables/groups/spacegroup.jl")
include("tables/groups/subperiodicgroup.jl")
include("tables/groups/mspacegroup.jl")
include("tables/generators/pointgroup.jl")
include("tables/generators/spacegroup.jl")
include("tables/generators/subperiodicgroup.jl")

include("assembly/groups/pointgroup.jl")
include("assembly/groups/spacegroup.jl")
include("assembly/groups/subperiodicgroup.jl")
include("assembly/groups/mspacegroup.jl")
export pointgroup, spacegroup, subperiodicgroup, mspacegroup

include("assembly/generators/pointgroup.jl")
include("assembly/generators/spacegroup.jl")
include("assembly/generators/subperiodicgroup.jl")
export generate, generators

include("show.jl") # printing of structs from src/[types.jl, types_symmetry_vectors.jl]

include("orders.jl")

include("symops.jl") # symmetry operations for space, plane, and line groups
export @S_str, compose,
       issymmorph, littlegroup, orbit,
       reduce_ops,
       issubgroup, isnormal,
       cosets

include("conjugacy.jl") # construction of conjugacy classes
export classes, is_abelian

include("wyckoff.jl") # wyckoff positions and site symmetry groups
export wyckoffs, WyckoffPosition,
       multiplicity,
       SiteGroup, sitegroup, sitegroups,
       cosets,
       findmaximal,
       siteirreps

include("symeigs2irrep.jl") # find irrep multiplicities from symmetry eigenvalue data
export find_representation

include("pointgroup.jl") # symmetry operations for crystallographic point groups
export pgirreps, PG_IUCs, find_isomorphic_parent_pointgroup

include("irreps_reality.jl")
export realify, realify!, calc_reality

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
export lgirreps, littlegroups

include("lattices.jl")
export ModulatedFourierLattice,
       getcoefs, getorbits, levelsetlattice,
       modulate, normscale, normscale!

include("compatibility.jl")
export subduction_count, remap_to_kstar

include("bandrep.jl")
export bandreps, classification, nontrivial_factors, basisdim

include("calc_bandreps.jl")
export calc_bandreps

include("deprecations.jl")
export get_littlegroups, get_lgirreps, get_pgirreps, WyckPos, kvec, wyck, kstar

include("grouprelations/grouprelations.jl")
export maximal_subgroups, minimal_supergroups, conjugacy_relations

# some functions are extensions of base-owned names; we need to (re)export them in order to 
# get the associated docstrings listed by Documeter.jl
export position, inv, isapprox

# ---------------------------------------------------------------------------------------- #
# EXTENSIONS AND JLD-FILE INITIALIZATION
if !isdefined(Base, :get_extension)
    using Requires # load extensions via Requires.jl on Julia versions <v1.9
end
# define functions we want to extend and have accessible via `Crystalline.(...)` if an
# extension is loaded
function _create_isosurf_plot_data end # implemented on CrystallinePyPlotExt load

## __init__
# - open .jld2 data files, so we don't need to keep opening/closing them
# - optional code-loading, using Requires.

# store the opened jldfiles in `Ref{..}`s for type-stability's sake (need `Ref` since we
# need to mutate them in `__init__` but cannot use `global const` in a function, cf.
# https://github.com/JuliaLang/julia/issues/13817)
const LGIRREPS_JLDFILES = ntuple(_ -> Ref{JLD2.JLDFile{JLD2.MmapIO}}(), Val(3))
const LGS_JLDFILES      = ntuple(_ -> Ref{JLD2.JLDFile{JLD2.MmapIO}}(), Val(3))
const PGIRREPS_JLDFILE  = Ref{JLD2.JLDFile{JLD2.MmapIO}}()

const DATA_DIR = joinpath(dirname(@__DIR__), "data")

function __init__()
    # open `LGIrrep` and `LittleGroup` data files for read access on package load (this
    # saves a lot of time compared to `jldopen`ing each time we call e.g. `lgirreps`,
    # where the time for opening/closing otherwise dominates)
    for D in (1,2,3)
        global LGIRREPS_JLDFILES[D][] =
            JLD2.jldopen(DATA_DIR*"/irreps/lgs/$(D)d/irreps_data.jld2", "r")
        global LGS_JLDFILES[D][] =
            JLD2.jldopen(DATA_DIR*"/irreps/lgs/$(D)d/littlegroups_data.jld2", "r")
    end
    global PGIRREPS_JLDFILE[] = # only has 3D data; no need for tuple over dimensions
            JLD2.jldopen(DATA_DIR*"/irreps/pgs/3d/irreps_data.jld2", "r")

    # ensure we close files on exit
    atexit(() -> foreach(jldfile -> close(jldfile[]), LGIRREPS_JLDFILES))
    atexit(() -> foreach(jldfile -> close(jldfile[]), LGS_JLDFILES))
    atexit(() -> close(PGIRREPS_JLDFILE[]))

    # load extensions via Requires.jl on Julia versions <v1.9
    @static if !isdefined(Base, :get_extension)
        @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" begin  
            include("../ext/CrystallinePyPlotExt.jl") # loads PyPlot and Meshing
            export mesh_3d_levelsetlattice
        end
        @require GraphMakie = "1ecd5474-83a3-4783-bb4f-06765db800d2" begin
            include("../ext/CrystallineGraphMakieExt.jl")
        end
    end
end

# precompile statements
if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end # module
