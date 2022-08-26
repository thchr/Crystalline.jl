
using Graphs

include(joinpath((@__DIR__), "types.jl"))
include(joinpath((@__DIR__), "show.jl"))
include(joinpath((@__DIR__), "init_data.jl"))

# ---------------------------------------------------------------------------------------- #

function group(g::GroupRelation{D,AG}, classidx::Int=1) where {D,AG}
    # TODO: instantiate a group of type AG with number `n.num`
    H = spacegroup(g.num, Val{D}())
    # transform `H` to the setting defined by `transforms[conjugacyclass]`
    t = g.classes[classidx]
    if isnothing(t.P) && isnothing(t.p)
        return H
    else
        P, p = something(t.P), something(t.p)
        P⁻¹ = inv(P)
        return typeof(H)(num(H), transform.(H, Ref(P⁻¹), Ref(-P⁻¹*p)))
    end
end

# technically, this concept of taking a "all the recursive children" of a graph vertex is
# called the "descendants" of the graph in 'networkx' (https://networkx.org/documentation/stable/reference/algorithms/dag.html)
# similarly, going the other way (supergroups) is called the "ancestors").
function maximal_subgroups(num::Integer, ::Type{SpaceGroup{D}}=SpaceGroup{3};
                                kind::GleicheKind = TRANSLATIONENGLEICHE) where D
    D ∈ (1,2,3) || _throw_invalid_dim(D)
    subgroupsd = (D == 3 ? SUBGROUPSD_3D :
                  D == 2 ? SUBGROUPSD_2D :
                  D == 1 ? SUBGROUPSD_1D : 
                  error("unreachable")    ) :: Dict{Int, GroupRelations{D, SpaceGroup{D}}}

    infos = Dict{Int, GroupRelations{D,SpaceGroup{D}}}()
    vertices = Int[num]
    remaining_vertices = Int[num]
    fadjlist = Vector{Int}[]
    while !isempty(remaining_vertices)
        num′ = pop!(remaining_vertices)
        
        info′ = subgroupsd[num′]
        children = filter(rel->rel.kind==kind, info′)
        infos[num′] = valtype(infos)(num′, children)

        push!(fadjlist, Int[])
        for child in children
            k = findfirst(==(child.num), vertices)
            if isnothing(k)
                push!(vertices, child.num)
                child.num ∉ remaining_vertices && pushfirst!(remaining_vertices, child.num)
                k = length(vertices)
            end
            push!(fadjlist[end], k)
        end
    end
    return GroupRelationGraph(vertices, fadjlist, infos, SUBGROUP)
end