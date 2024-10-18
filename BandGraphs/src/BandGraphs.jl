module BandGraphs

using Crystalline
using Crystalline: irdim, AbstractSymmetryVector
import JLD2
using PrettyTables
using LinearAlgebra
using Graphs: Graph, Edge, add_vertex!, add_edge!
using MetaGraphsNext
using Combinatorics: permutations # for `set_pinned_subgraphs!`
using BlockArrays: BlockArray, Block # for `assemble_adjacency`

# ---------------------------------------------------------------------------------------- #
# EXPORTS

export
    Connection,
    SubductionTable,
    Partition,
    SubGraph,
    BandGraph,
    chinese_postman,
    build_subgraphs,
    assemble_adjacency,
    assemble_degree,
    assemble_laplacian,
    assemble_graph,
    unfold_bandgraph,
    partition_graph,
    split_nonmaximal_nodes,
    permute_subgraphs,
    subduction_tables,
    plot_flattened_bandgraph,
    findall_separable_vertices,
    complete_split

# ---------------------------------------------------------------------------------------- #

# TODO: improve sharing of types between `BandGraphs` & `build/crawl_and_write_bandpaths.jl`
include("subduction-types.jl")
include("subduction-table.jl")

# ---------------------------------------------------------------------------------------- #
# DATA

const SUBDUCTIONSD_TR_3D = Dict{Int, Vector{SubductionTable{3}}}() # `timereversal = true`
const SUBDUCTIONSD_3D = Dict{Int, Vector{SubductionTable{3}}}()    # `timereversal = false`
const SUBDUCTIONSD_TR_2D = Dict{Int, Vector{SubductionTable{2}}}() # `timereversal = true`
const SUBDUCTIONSD_2D = Dict{Int, Vector{SubductionTable{2}}}()    # `timereversal = false`
function __init__()
    # TODO: Fix the awfulness here: the loading below only works if the dataset is created
    #       with the types from BandGraphs first - not a good circular thing to have...
    datapath = joinpath(dirname(@__DIR__), "data/connections/3d/subductions-tr.jld2")
    JLD2.jldopen(datapath, "r") do jldfile
        for (k,v) in jldfile["subductionsd"]
            SUBDUCTIONSD_TR_3D[k] = v
        end
    end
    datapath = joinpath(dirname(@__DIR__), "data/connections/3d/subductions.jld2")
    JLD2.jldopen(datapath, "r") do jldfile
        for (k,v) in jldfile["subductionsd"]
            SUBDUCTIONSD_3D[k] = v
        end
    end
    datapath = joinpath(dirname(@__DIR__), "data/connections/2d/subductions-tr.jld2")
    JLD2.jldopen(datapath, "r") do jldfile
        for (k,v) in jldfile["subductionsd"]
            SUBDUCTIONSD_TR_2D[k] = v
        end
    end
    datapath = joinpath(dirname(@__DIR__), "data/connections/2d/subductions.jld2")
    JLD2.jldopen(datapath, "r") do jldfile
        for (k,v) in jldfile["subductionsd"]
            SUBDUCTIONSD_2D[k] = v
        end
    end
end

"""
    subduction_tables(sgnum, ::Val{D}; timereversal=true)  --> Vector{SubductionTable{D}}

Return a vector of `SubductionTable`s from stored tabulations, with (`timereversal = true`)
or without (`timereversal = false`) time-reversal symmetry in space group `sgnum`.
"""
function subduction_tables(sgnum::Integer, ::Val{D}=Val(3); timereversal::Bool=true) where D
    if D == 3
        if timereversal
            return SUBDUCTIONSD_TR_3D[sgnum]::Vector{SubductionTable{3}}
        else
            return SUBDUCTIONSD_3D[sgnum]::Vector{SubductionTable{3}}
        end
    elseif D == 2
        if timereversal
            return SUBDUCTIONSD_TR_2D[sgnum]::Vector{SubductionTable{2}}
        else
            return SUBDUCTIONSD_2D[sgnum]::Vector{SubductionTable{2}}
        end
    else
        throw(DomainError(D, "dimension must be 1, 2, or 3"))
    end
end
function subduction_tables(sgnum::Integer, D::Integer; timereversal::Bool=true)
    return subduction_tables(sgnum, Val(D); timereversal)
end

# ---------------------------------------------------------------------------------------- #

include("types.jl")
include("subgraphs.jl")
include("matrices.jl")
include("graphs.jl")
include("subgraph-permutations.jl")
include("graph_routing.jl")
include("unfold.jl")
include("complete-split.jl")
include("subsetsum.jl")
include("separable.jl")

# ---------------------------------------------------------------------------------------- #
# EXTENSIONS: loaded via Requires.jl on Julia versions <v1.9; otherwise via the Pkg
#             extension system

if !isdefined(Base, :get_extension)
    using Requires
end

# we cannot export directly from an extension, so we must define & export methods here and
# then actually fill in extension module
function plot_flattened_bandgraph end
function make_vertices_dragable! end
export plot_flattened_bandgraph, make_vertices_dragable!
@static if !isdefined(Base, :get_extension)
    function __init__()
        @require GraphMakie="1ecd5474-83a3-4783-bb4f-06765db800d2" begin
            include("../ext/BandGraphsGraphMakieExt.jl")
        end
    end
end

# ---------------------------------------------------------------------------------------- #

end # module