
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
# called the "descendants" of the graph in 'networkx'; similarly, going the other way 
# (supergroups) is called the "ancestors":
# see https://networkx.org/documentation/stable/reference/algorithms/dag.html)
for (f, rel_shorthand, rel_kind, doctxt) in
                            ((:maximal_subgroups,   :MAXSUB, :SUBGROUP,  "maximal sub"),
                             (:minimal_supergroups, :MINSUP, :SUPERGROUP, "minimal super"))
    relations_3D = Symbol(rel_shorthand, :_RELATIONSD_3D)
    relations_2D = Symbol(rel_shorthand, :_RELATIONSD_2D)
    relations_1D = Symbol(rel_shorthand, :_RELATIONSD_1D)
    @eval begin
        @doc """
            $($f)(num::Integer, AG::Type{<:AbstractGroup}=SpaceGrop{3}; kind)
        
        Returns the graph structure of the $($doctxt)groups as a `GroupRelationGraph` for a
        group of type `AG` and number `num`.

        ## Visualization

        The resulting group structure can be plotted using Makie.jl (e.g., GLMakie.jl) using
        `plot(::GroupRelationGraph)`:
        ```jl
        julia> using Crystalline
        julia> gr = $($f)(112, SpaceGroup{3})
        julia> using GraphMakie, GLMakie
        julia> plot(gr)
        ```
        ## Keyword arguments

        - `kind` (default, `Crystalline.TRANSLATIONENGLEICHE`): to return klassengleiche
          relations, set `kind = Crystalline.KLASSENGLEICHE`). For klassengleiche
          relationships, only a selection of reasonably low-index relationships are
          returned.

        ## Data sources

        The group relationships returned by this function were retrieved from the Bilbao
        Crystallographic Server's [MAXSUB](https://www.cryst.ehu.es/cryst/maxsub.html)
        program. Please cite the original reference work associated with MAXSUB:

        - Aroyo et al., [Z. Kristallogr. Cryst. Mater. **221**, 15 (2006)](https://doi.org/10.1524/zkri.2006.221.1.15).
        """
        function $f(num::Integer, ::Type{SpaceGroup{D}}=SpaceGroup{3};
                    kind::GleicheKind = TRANSLATIONENGLEICHE) where D
            D ∈ (1,2,3) || _throw_invalid_dim(D)
            relationsd = (D == 3 ? $relations_3D :
                        D == 2 ? $relations_2D :
                        D == 1 ? $relations_1D : error("unreachable")
                        ) :: Dict{Int, GroupRelations{D, SpaceGroup{D}}}

            infos = Dict{Int, GroupRelations{D,SpaceGroup{D}}}()
            vertices = Int[num]
            remaining_vertices = Int[num]
            fadjlist = Vector{Int}[]
            while !isempty(remaining_vertices)
                num′ = pop!(remaining_vertices)
                
                info′ = relationsd[num′]
                children = filter(rel->rel.kind==kind, info′)
                infos[num′] = valtype(infos)(num′, children)

                push!(fadjlist, Int[])
                for child in children
                    k = findfirst(==(child.num), vertices)
                    if isnothing(k)
                        push!(vertices, child.num)
                        if child.num ∉ remaining_vertices
                            pushfirst!(remaining_vertices, child.num)
                        end
                        k = length(vertices)
                    end
                    push!(fadjlist[end], k)
                end
            end
            return GroupRelationGraph(vertices, fadjlist, infos, $rel_kind)
        end # function $f
    end # begin @eval
end # for (f, ...)


"""
    conjugacy_relations(gr::GroupRelationGraph, sgnumᴳ, sgnumᴴ) --> Vector{<:ConjugacyTransform}

Given a graph `gr` representing a sub- or supergroup relation, return the possible
transformations that make up the conjugacy classes between the groups `sgnumᴳ` and `sgnumᴴ`.

The returned transforms `ts` bring `G` into the setting of `H` via `transform.(G, Ref(t.P),
Ref(t.p))`, where `t` is an element of `ts`, i.e., one possible conjugacy transform, and
`t.P` and `t.p` denote the associated transform's rotation and translation, respectively. 

Note that the transforms need not preserve volume: accordingly, some operations may be
redundant after transformation (use [`reduce_ops`](@ref) or `unique!` to remove these).

## Example
Consider the following example, which looks at the subgroup relationship between space
groups ⋕168 and ⋕3:
```jldoctest conjugacy_relation
julia> gr = maximal_subgroups(168)
GroupRelationGraph (subgroups) of SpaceGroup{3} ⋕168 with 4 vertices:
 ⋕168
 ├─►⋕143ᵀ (index 2)  ╌╌►⋕1ᵀ (index 3)
 └─►⋕3ᵀ (index 3)  ╌╌►⋕1ᵀ (index 2)

julia> sg168 = spacegroup(168)
SpaceGroup{3} ⋕168 (P6) with 6 operations:
 1
 3₀₀₁⁺
 3₀₀₁⁻
 2₀₀₁
 6₀₀₁⁻
 6₀₀₁⁺

julia> sg3 = spacegroup(3)
SpaceGroup{3} ⋕3 (P2) with 2 operations:
 1
 2₀₁₀
```
Note that the symmetry operation 2₀₁₀ (from ⋕3) and 2₀₀₁ (from ⋕168) differ by a
transformation; even though they are isomorphic, this is not clearly reflected because they
are in different settings. We can use `conjugacy_relations` to find the transformations that
brings ⋕168 into the setting of ⋕3:
```jldoctest conjugacy_relation
julia> ts = conjugacy_relations(gr, 168, 3) # possible transforms from ⋕168 to ⋕3
1-element Vector{Crystalline.ConjugacyTransform{3}}:
 P = [0 0 1; 1 0 0; 0 1 0]

julia> sg168′ = transform.(sg168, Ref(ts[1].P), Ref(ts[1].p))
6-element Vector{SymOperation{3}}:
 1
 3₀₁₀⁺
 3₀₁₀⁻
 2₀₁₀
 6₀₁₀⁻
 6₀₁₀⁺

julia> issubgroup(sg168, sg3) # not identified as subgroup, since settings differ
false

julia> issubgroup(sg168′, sg3) # settings now agree, and subgroup relationship is evident
true
```
Here, there is only one possible transformation: in general, however, there may be many.
"""
function conjugacy_relations(
            gr::GroupRelationGraph{D},
            sgnumᴳ::Integer,
            sgnumᴴ::Integer) where D

    idxᴳ = findfirst(==(sgnumᴳ), gr.nums)
    idxᴴ = findfirst(==(sgnumᴴ), gr.nums)
    idxᴳ !== nothing || error(lazy"`gr` does not contain group $sgnumᴳ")
    idxᴴ !== nothing || error(lazy"`gr` does not contain group $sgnumᴴ")
    has_path(gr, idxᴳ, idxᴴ) || error(lazy"`gr` does not connect group $sgnumᴳ to $sgnumᴴ")

    P₀ = one(SMatrix{D, D, Float64})
    p₀ = zero(SVector{D, Float64})
    sgnumᴳ == sgnumᴴ && return [ConjugacyTransform{D}(P₀, p₀)] # sgnumᴳ == sgnumᴴ case

    idxsᴳ²ᴴs = all_simple_paths(gr, idxᴳ, idxᴴ)
    classes = Vector{ConjugacyTransform{D}}()
    for idxsᴳ²ᴴ in idxsᴳ²ᴴs
        accumulate_classes_recur!(classes, gr, P₀, p₀, idxsᴳ²ᴴ)
    end

    return unique!(classes)
end

# recursively combine the composition of the distinct conjugacy classes as we traverse the
# path of transformations; modifies `classes`
function accumulate_classes_recur!(
            classes::Vector{ConjugacyTransform{D}}, gr::GroupRelationGraph{D}, 
            P::SMatrix{D, D}, p::SVector{D}, idxsᴳ²ᴴ) where D
    
    if length(idxsᴳ²ᴴ) == 1 # base case
        t = ConjugacyTransform{D}(P, p)
        t ∉ classes && push!(classes, t)
        return classes
    end

    idxᴳ′ = idxsᴳ²ᴴ[1]
    idxᴴ′ = idxsᴳ²ᴴ[2]
    sgnumᴳ′ = gr.nums[idxᴳ′]
    sgnumᴴ′ = gr.nums[idxᴴ′]
    rels = gr[sgnumᴳ′] # ::GroupRelations
    iᴴ′ = something(findfirst(r->r.num==sgnumᴴ′, rels.children))
    rel = rels[iᴴ′] # ::GroupRelation
    rel.kind == TRANSLATIONENGLEICHE || error("Klassengleiche relations are not supported")
    classes′ = rel.classes
    for ctransform in classes′
        # composition of transformations work differently for sub- and supergroup relations,
        # because they are applied in a different order
        if gr.direction == SUBGROUP
            p′ = P * ctransform.p + p
            P′ = P * ctransform.P
        else # SUPERGROUP
            p′ = ctransform.P * p + ctransform.p
            P′ = ctransform.P * P
        end
        accumulate_classes_recur!(classes, gr, P′, p′, @view idxsᴳ²ᴴ[2:end])
    end
end