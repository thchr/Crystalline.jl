# `show` & `summary`
is_tgleiche(kind::GleicheKind) = kind == TRANSLATIONENGLEICHE

function Base.summary(io::IO, rels::GroupRelations{D,AG}) where {D,AG}
    print(io, "GroupRelations{", AG, "} ⋕", rels.num, " (", iuc(rels.num, D), 
              ") with ", length(rels), " elements")
end

function Base.summary(io::IO, rel::GroupRelation{D,AG}) where {D,AG}
    print(io, "GroupRelation{", AG, "} ⋕", rel.num)
    printstyled(io, is_tgleiche(rel.kind) ? "ᵀ" : "ᴷ"; 
                    color = is_tgleiche(rel.kind) ? :green : :red)
    print(io, " (", iuc(rel.num, D), ") with ", length(rel.classes), " conjugacy class", 
              length(rel.classes) == 1 ? "" : "es")
end

function Base.show(io::IO, ::MIME"text/plain", rels::GroupRelations{D,AG}) where {D,AG}
    AG <: SpaceGroup || error("correct show not yet implemented for this type")
    summary(io, rels)
    print(io, ":")
    for rel in values(rels)
        print(io, "\n ")
        _print_child(io, rel.num, rel.kind, rel.index)
    end
end

function _print_child(
            io::IO, 
            num::Int,
            kind::Union{Nothing, GleicheKind},
            index::Union{Nothing, Int})
    print(io, "⋕", num)
    if !isnothing(kind)
        printstyled(io, is_tgleiche(kind) ? "ᵀ" : "ᴷ";
                        color = is_tgleiche(kind) ? :green : :red)
    end
    if !isnothing(index)
        printstyled(io, " (index ", index, ")"; color=:light_black)
    end
end

function Base.show(io::IO, rels::GroupRelations{D,AG}) where {AG,D}
    print(io, "[")
    Nrels = length(rels)
    for (i, rel) in enumerate(values(rels))
        print(io, rel.num)
        printstyled(io, is_tgleiche(rel.kind) ? "ᵀ" : "ᴷ";
                        color = is_tgleiche(rel.kind) ? :green : :red)
        i == Nrels || print(io, ", ")
    end
    print(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", g::GroupRelation)
    summary(io, g)
    print(io, ":")
    for (i, ct) in enumerate(g.classes)
        printstyled(io, "\n ", Crystalline.supscriptify(string(i)), "⁾ ", color=:light_black)
        print(io, ct)
    end
end

function Base.show(io::IO, t::ConjugacyTransform{D}) where D
    print(io, "P = ")
    if isnothing(t.P) && isnothing(t.p)
        print(io, 1)
    else
        P, p = something(t.P), something(t.p)
        if !isone(P)
            print(io, replace(string(rationalize.(float(P))), r"//1([ |\]|;])"=>s"\1", "//"=>"/", "Rational{Int64}"=>""))
        else
            print(io, 1)
        end
        if !iszero(p)
            print(io, ", p = ",
                      replace(string(rationalize.(float(p))), r"//1([,|\]])"=>s"\1", "//"=>"/", "Rational{Int64}"=>""))
        end
    end
end

function Base.summary(io::IO, gr::GroupRelationGraph{D, AG}) where {D,AG}
    print(io, "GroupRelationGraph",
              " (", gr.direction == SUPERGROUP ? "super" : "sub", "groups)",
              " of ", AG, " ⋕", first(gr.nums), 
              " with ", length(gr.nums), " vertices")
end
function Base.show(io::IO, gr::GroupRelationGraph)
    summary(io, gr)
    print(io, ":")
    num = first(gr.nums) # find "base" group number
    _print_graph_leaves(io, gr, num, Bool[], nothing, nothing)
end
function _print_graph_leaves(
            io::IO, 
            gr::GroupRelationGraph, 
            num::Int, 
            is_last::Vector{Bool}, 
            kind::Union{Nothing, GleicheKind}, 
            index::Union{Nothing, Int}
            )

    # indicate "structure"-relationship of graph "leaves"/children via box-characters
    indent = length(is_last)
    print(io, "\n ")
    if indent ≥ 2
        for l in 1:indent-1
            if !is_last[l]
                printstyled(io, "│"; color=:light_black)
                print(io, "  ")
            else
                print(io, "   ")
            end
        end
    end
    if indent > 0
        printstyled(io, last(is_last) ? "└" : "├", "─►"; color=:light_black)
    end

    # print info about child
    _print_child(io, num, kind, index)

    # now move on to children of current child; print info about them, recursively
    if !isnothing(kind) && kind == KLASSENGLEICHE
        # avoid printing recursively for klassengleiche as this implies infinite looping
        # NB: in principle, this means we don't show the entire structure of 
        #     `gr.infos`; but it is also really not that clear how to do this
        #     properly for klassengleiche relations which are not a DAG but often 
        #     cyclic

    # only one child and child is trivial group
    elseif length(gr.infos[num].children) == 1 && gr.infos[num].children[1].num == 1
        child = gr.infos[num].children[1]
        printstyled(io, "  ╌╌►"; color=:light_black)
        _print_child(io, child.num, child.kind, child.index)
    
    # multiple children/nontrivial children
    else
        N_children = length(gr.infos[num].children)
        for (i,child) in enumerate(gr.infos[num].children)
            child_is_last = vcat(is_last, i==N_children)
            _print_graph_leaves(io, gr, child.num, child_is_last, child.kind, child.index)
        end
    end
end