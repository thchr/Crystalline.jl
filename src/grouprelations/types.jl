# ---------------------- #
@enum GleicheKind::Int8 TRANSLATIONENGLEICHE KLASSENGLEICHE
@enum RelationKind::Int8 SUBGROUP SUPERGROUP

# ---------------------- #
struct ConjugacyTransform{D} # transform to a specific conjugacy class
    P :: Union{Nothing, SqSMatrix{D,Float64}}
    p :: Union{Nothing, SVector{D,Float64}}
end

# ---------------------- #
struct GroupRelation{D,AG} # all conjugacy classes of a subgroup
    parent  :: Int # the supergroup if num indicates a subgroup, and vice versa
    num     :: Int # group label; aka node label
    index   :: Int # index of H in G or vice versa, depending on context
    kind    :: GleicheKind
    classes :: Vector{ConjugacyTransform{D}} # isomorphisms over distinct conjugacy classes
end

function GroupRelation{AG}(parent::Int, num::Integer, index::Integer, kind::GleicheKind,
    transforms::Union{Nothing, Vector{ConjugacyTransform{D}}}
    ) where {D,AG}
return GroupRelation{D,AG}(parent, num, index, kind, transforms)
end

# ---------------------- #
struct GroupRelations{D,AG} <: AbstractVector{GroupRelation{D,AG}}
    num      :: Int # group label; aka node label
    children :: Vector{GroupRelation{D,AG}} # list of sub- or supergroups
end

function GroupRelations{AG}(num::Integer, gs::AbstractVector{GroupRelation{D,AG}}
    ) where {D,AG}
return GroupRelations{D,AG}(num, gs)
end
Base.getindex(rels::GroupRelations, i::Integer) = getindex(rels.children, i)
Base.size(rels::GroupRelations) = size(rels.children)

# ---------------------- #
struct GroupRelationGraph{D,AG} <: Graphs.AbstractGraph{Int}
    nums      :: Vector{Int} # group numbers associated w/ each graph vertex index
    fadjlist  :: Vector{Vector{Int}} # indexes into `nums` conceptually
    infos     :: Dict{Int, GroupRelations{D,AG}} # sub/supergroup relation details
    direction :: RelationKind # whether this is a sub- or supergroup relation
end
dim(::GroupRelationGraph{D}) where D = D

# ---------------------------------------------------------------------------------------- #
# Graphs.jl `AbstractGraph` interface for `GroupRelationGraph`

Graphs.vertices(gr::GroupRelationGraph) = eachindex(gr.nums)
Base.eltype(::GroupRelationGraph) = Int
Graphs.edgetype(::GroupRelationGraph) = Graphs.SimpleEdge{Int}
function Graphs.has_edge(gr::GroupRelationGraph, src::Integer, dst::Integer)
    return Graphs.has_vertex(gr, src) && dst ∈ gr.fadjlist[src]
end
Graphs.has_vertex(gr::GroupRelationGraph, src::Integer) = src ∈ Graphs.vertices(gr)
function Graphs.outneighbors(gr, src::Integer) # every index "outgoing" upon `src`
    Graphs.has_vertex(gr, src) || throw(DomainError(src, "`src` is not in graph"))
    return gr.fadjlist[src]
end
function Graphs.inneighbors(gr, src::Integer) # every index "incident" upon `src`
    Graphs.has_vertex(gr, src) || throw(DomainError(src, "`src` is not in graph"))
    idxs = findall(∋(src), gr.fadjlist)
    return idxs
end
Graphs.ne(gr::GroupRelationGraph) = sum(length, gr.fadjlist)
Graphs.nv(gr::GroupRelationGraph) = length(Graphs.vertices(gr))
Graphs.is_directed(::Union{Type{<:GroupRelationGraph}, GroupRelationGraph}) = true 

Graphs.edges(gr::GroupRelationGraph) = GroupRelationGraphEdgeIter(gr)
struct GroupRelationGraphEdgeIter{D,AG}
    gr :: GroupRelationGraph{D,AG}
end
function Base.iterate(it::GroupRelationGraphEdgeIter, 
                      state::Union{Nothing, Tuple{Int,Int}}=(1,1))
    # breadth-first iteration
    i,j = state
    if length(it.gr.fadjlist[i]) < j
        if length(it.gr.fadjlist) > i
            return iterate(it, (i+1, 1)) # go to next vertex
        else
            return nothing # no more edges or vertices
        end
    else
        return Graphs.Edge(i, it.gr.fadjlist[i][j]), (i, j+1)
    end
end
Base.length(it::GroupRelationGraphEdgeIter) = Graphs.ne(it.gr)