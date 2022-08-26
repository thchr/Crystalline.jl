
using Crystalline, HTTP, Gumbo, LinearAlgebra

# getting the maximal subgroups of t type

# getting the matching transformation matrix
# https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-tranmax?way=up&super=2&sub=1&index=2&client=maxsub&what=&path=&series=&conj=all&type=t
# simpler: 
# https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-tranmax?way=up&super=2&sub=1&type=t

""" 
    crawl_maximal_subgroups(numᴳ::Integer; subgroup_kind=:tk))

Obtains the maximal subgroups for a given space group G, with number `numᴳ`, by crawling
the Bilbao server, i.e. find the space groups {Hⱼ} such that Hⱼ < G with _no_ Hⱼ < Hₖ.
Only works for 3D space groups.

The keyword argument `subgroup_kind` specifies whether to return the translationengleiche
("t"), klassengleiche ("k"), or either ("tk") subgroup kind.
"""
function crawl_maximal_subgroups(numᴳ::Integer, D::Integer=3;
                                 subgroup_kind::Symbol=:tk, verbose::Bool=true)
    if subgroup_kind ∉ (:t, :k, :tk)
        throw(DomainError(subgroup_kind, "the subgroup type must be :t, :k, or :tk"))
    end
    htmlraw = http_request_maximal_subgroups(numᴳ, D)

    # grab relevant rows of table
    body_html = last(children(htmlraw.root))
    table_html = if D == 3
        children(body_html[6][1][1])[2:end]
    elseif D == 2
        children(children(body_html[4][3][1])[1])[2:end]
    end
    # grab relevant columns of table & extract xyzt-forms of generators
    numsᴴ = parse.(Int, _grab_column(table_html, 2)) :: Vector{Int}
    indices = parse.(Int, _grab_column(table_html, 4)) :: Vector{Int}
    kinds = Symbol.(_grab_column(table_html, 5)) :: Vector{Symbol}

    if subgroup_kind != :tk
        subgroup_kind_indices = findall(==(subgroup_kind), kinds)
        keepat!(numsᴴ, subgroup_kind_indices)
        keepat!(indices, subgroup_kind_indices)
        keepat!(kinds, subgroup_kind_indices)
    end
    Ppss = find_transformation_matrices.(numsᴴ, numᴳ, D, kinds, indices) # vector of {P|p} transformations
    verbose && println(stdout, "crawled ", D == 3 ? "space" : "plane", " group ", numᴳ,
                               " (", length(numsᴴ), " subgroups)")
    return [(;num, index, kind, Pps) for (num, index, kind, Pps) in zip(numsᴴ, indices, kinds, Ppss)]
end

function _grab_column(table_html::Vector{HTMLNode}, i::Integer)
    strip.(Gumbo.text.(only.(children.(getindex.(table_html, i))))) :: Vector{SubString{String}}
end

function http_request_maximal_subgroups(numᴳ::Integer, D::Integer=3)
    (numᴳ < 1 || numᴳ > 230) && error(DomainError(numᴳ))
    bilbaourl = "https://www.cryst.ehu.es/cgi-bin/"
    program_spec = if D == 3
        "cryst/programs/nph-lxi?client=maxsub&way=up&type=t&gnum="
    elseif D == 2
        "plane/programs/nph-plane_maxsub?gnum="
    else
        throw(DomainError(D))
    end
    contents = HTTP.request("GET", bilbaourl * program_spec * string(numᴳ))
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #

function find_transformation_matrices(numᴴ::Integer, numᴳ::Integer, D::Integer,
                                      kind::Symbol, index::Integer)
    htmlraw = http_request_subgroup_transformation(numᴴ, numᴳ, D, kind, index)

    # grab relevant rows of table(s)
    body_html = children(last(children(htmlraw.root)))
    multiple_conjugacy_classes = (Gumbo.text(body_html[5]) == 
                                  "  The subgroups of the same type can be divided in\n  ")
    if !multiple_conjugacy_classes
        # all transformations in the single (only) conjugacy class
        tables_html = [children(body_html[5][5][2])[2:end]]
        # Note that a single conjugacy class may still have multiple distinct transformations
        # (e.g., mapping H onto different conjugate parts of G); we are generally only interested
        # in the first mapping, so we just grab the first row (the first transformation; 
        # usually the simplest). 
    else
        # transformations to multiple distinct conjugacy classes: include a representative
        # of each conjugacy class. Different conjugacy class transformations map to
        # different conjugacy classes of the supergroup; so we should include both when
        # we think about implications for topology (because the symmetry eigenvalues in
        # distinct conjugacy classes can differ)
        table_idx = 5
        tables_html = Vector{HTMLNode}[]
        while (table_idx ≤ length(children(body_html[9])) && 
               (subbody_html = body_html[9][table_idx];
               begins_with_conjugacy_str = startswith(Gumbo.text(subbody_html[1]), "Conjugacy class");
               begins_with_conjugacy_str || table_idx == 5) )
            subbody_idx = 1+begins_with_conjugacy_str
            push!(tables_html, children(subbody_html[subbody_idx])[2:end]) # transformation
            table_idx += 2
        end
    end
    # The transformations {P|p} consist of a rotation part P and a translation part p;
    # notably, the rotation part P may also contain a scaling part (i.e. detP is not
    # necessarily 1)
    Pps = Vector{Tuple{Matrix{Rational{Int}}, Vector{Rational{Int}}}}(undef, length(tables_html))
    for (i,table_html) in enumerate(tables_html)
        row_html = table_html[1]
        Pp_str = replace(replace(Gumbo.text(row_html[2]), r"\n *"=>'\n'), r" +"=>',')
        Pp_str_rows = split.(split(Pp_str, '\n'), ',')
        P = [parse_maybe_fraction(Pp_str_rows[row][col]) for row in 1:D, col in 1:D]
        p = [parse_maybe_fraction(Pp_str_rows[row][D+1]) for row in 1:D]
        Pps[i] = (P, p)
    end

    # the transformation is intended as acting _inversely_ on the elements of H
    # (equivalently, directly on the elements of g∈G that have isomorphic elements in H);
    # that is, for every h∈H, there exists a mapped h′∈H′ with H′<G in the exact sense, with
    #       h′ = {P|p} h {P|p}⁻¹ = `transform(h, inv(P), -inv(W)*p)`
    # since det P is not necessarily 1, the transformation can collapse multiple elements
    # onto eachother, such that some h′ become equivalent in the new lattice basis; this is
    # e.g. the case in subgroup 167 of space group 230; to fix this, the result could be 
    # passed to `reduce_ops`. Hence, the relation we get it:
    #       H′ = [transform(h, inv(P), -inv(P)*p) for h in H]
    #       H′_reduced = reduce_ops(H′, centering(G))
    #       G_reduced = reduce_ops(G, centering(G))
    #       issubgroup(H′_reduced, G_reduced) == true
    return Pps
end

function parse_maybe_fraction(s)
    slash_idxs = findfirst("/", s)
    if isnothing(slash_idxs)
        return Rational{Int}(parse(Int, s)::Int)
    else
        slash_idx = slash_idxs[1] # only takes up one codepoint
        parse(Int, s[1:slash_idx-1])::Int // parse(Int, s[slash_idx+1:end])::Int
    end
end

function http_request_subgroup_transformation(numᴴ::Integer, numᴳ::Integer, # H < G
                                              D::Integer, kind::Symbol, index::Integer) 
    (numᴳ < 1 || numᴳ > 230) && error(DomainError(numᴳ))
    (numᴴ < 1 || numᴴ > 230) && error(DomainError(numᴴ))

    baseurl = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-tranmax?way=up&"
    url = baseurl * "type=$(kind)&super=$(numᴳ)&sub=$(numᴴ)&index=$(index)"
    if D == 2
        url *= "&client=planesubmax&subtype=plane"
    end

    contents = HTTP.request("GET", url)
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #
# Crawl and save all the data (~10-15 min)
@time begin
subgroup_data = ([# === 1D === (done manually)
                  # line group 1
                  [(num = 1, index = 2, kind = :k, Pps = [([2//1;;], [0//1])]),
                   (num = 1, index = 3, kind = :k, Pps = [([3//1;;], [0//1])]) ],
                  # line group 2
                  [(num = 1, index = 2, kind = :t, Pps = [([1//1;;], [0//1])]),
                   (num = 2, index = 2, kind = :k, Pps = [([2//1;;], [0//1]), 
                                                          ([2//1;;], [1//2])]),
                   (num = 2, index = 3, kind = :k, Pps = [([3//1;;], [0//1])]) ]  ],
                  # === 2D ===
                  crawl_maximal_subgroups.(1:17, 2),
                  # === 3D ===
                  crawl_maximal_subgroups.(1:230, 3))
end

using JLD2
jldsave("data/spacegroup_subgroups_data.jld2"; subgroups_data=subgroups_data)

# ---------------------------------------------------------------------------------------- #
# Test that H and G are indeed subgroups under the provided transformations

using Test
@testset for D in 1:3
    println("=== D = ", D, " ===")
    subgroups_data_D = subgroups_data[D]
    for numᴳ in 1:MAX_SGNUM[D]
        println("   G = ", numᴳ)
        G = spacegroup(numᴳ, D)
        G_reduced = reduce_ops(G, centering(G))
        subgroup_data = subgroups_data_D[numᴳ]
        numsᴴ = getindex.(subgroup_data, Ref(:num))
        Ppss  = getindex.(subgroup_data, Ref(:Pps))
        kinds = getindex.(subgroup_data, Ref(:kind))
        for (numᴴ, kind, Pps) in zip(numsᴴ, kinds, Ppss)
            kind == :k && continue # skip klassengleiche
            println("     H = ", numᴴ)
            H = spacegroup(numᴴ, D)
            for (c, (P, p)) in enumerate(Pps)
                length(Pps) > 1 && println("    Conjugacy class ", Char(Int('a')-1+c))
                H′ = transform.(H, Ref(inv(P)), Ref(-inv(P)*p))
                H′_reduced = reduce_ops(H′, centering(G))
                @test issubgroup(G_reduced, H′_reduced)
            end
        end
    end
end
