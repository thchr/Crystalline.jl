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
        print(io, '\n')
        _print_child(io, rel.num, rel.kind, rel.index)
    end
end

function _print_child(io::IO, num::Int, kind::GleicheKind, index::Int)
    print(io, " ⋕", num)
    printstyled(io, is_tgleiche(kind) ? "ᵀ" : "ᴷ", color = is_tgleiche(kind) ? :green : :red)
    print(io, " (index ", index, ")")
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