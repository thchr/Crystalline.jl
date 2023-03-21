module BandGraphs

using Crystalline
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
    SymVector,
    Partition,
    SubGraph,
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
    SUBDUCTIONSD,
    plot_flattened_bandgraph

# ---------------------------------------------------------------------------------------- #

# TODO: improve sharing of types between `BandGraphs` & `build/crawl_and_write_bandpaths.jl`
include("subduction-types.jl")

# ---------------------------------------------------------------------------------------- #
# DATA

# TODO: FIX
const SUBDUCTIONSD = Dict{Int, Vector{SubductionTable{3}}}()
function __init__()
    # TODO: Fix the awfulness here: the loading below only works if the dataset is created
    #       with the types from BandGraphs first - not a good circular thing to have...
    datapath = joinpath(dirname(@__DIR__), "data/connections/3d/subductions-tr.jld2")
    JLD2.jldopen(datapath, "r") do jldfile
        for (k,v) in jldfile["subductionsd"]
            SUBDUCTIONSD[k] = v
        end
    end
end

# ---------------------------------------------------------------------------------------- #

include("types.jl")
include("subgraphs.jl")
include("matrices.jl")
include("graphs.jl")
include("subgraph-permutations.jl")
include("graph_routing.jl")
include("unfold.jl")

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