using HTTP, Gumbo

# Small convenience script to crawl and subsequently write all the xyzt forms of the
# general positions/generators of subperiodic groups.

# ---------------------------------------------------------------------------------------- #
## Functions for crawling the generators of 3D space groups from Bilbao

""" 
    crawl_subperiodic_xyzt(num::Integer, D::Integer=3)

Obtains the general positions or generators in xyzt format for a given subperiodic group
number `num` by crawling the Bilbao server. 
"""
function crawl_subperiodic_xyzt(num::Integer, D::Integer=3, P::Integer=2,
                                kind::String="operations")
    htmlraw = crawl_subperiodic_html(num, D, P, kind)

    # layout of html page depends on whether kind is `"operations"` or `"generators"`, so we
    # pick a different html element in the structure depending on this
    idx = kind == "generators" ? 5 : kind == "operations" ? 4 : error(DomainError(kind))
    # grab relevant rows of table
    table_html = children(children(children(children(last(children(htmlraw.root)))[idx])[1])[1])[3:end]
    # grab relevant columns of table & extract xyzt-forms of generators
    gens_str = strip.(text.(only.(children.(getindex.(table_html, 2)))))
    
    return gens_str
end

function subperiodic_url(num::Integer, D::Integer, P::Integer, kind::String)
    sub = if D==3 && P==2
        (num < 1 || num > 80) && error(DomainError(num))
        "layer"
    elseif D==3 && P==1
        (num < 1 || num > 75) && error(DomainError(num))
        "rod"
    elseif D==2 && P==1
        (num < 1 || num > 7) && error(DomainError(num))
        "frieze"
    else
        error("cannot crawl subperiodic groups of dimensionality $D and periodicity $P")
    end

    return "https://cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen?" *
           "what=$(kind == "generators" ? "gen" : kind == "operations" ? "gp" : error())" *
           "&gnum=$(num)&subtype=$(sub)&"
end

function crawl_subperiodic_html(num::Integer, D::Integer, P::Integer, kind::String)
    contents = HTTP.request("GET", subperiodic_url(num, D, P, kind))
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #
## Use crawling functions & write information to:
#      `/test/data/xyzt-operations/generators/subperiodic/{layer|rod|frieze}/`
#      `/test/data/xyzt-operations/subperiodic/{layer|rod|frieze}/`
for (D, P, sub, maxnum) in [(3,2,"layer",80), (3,1,"rod",75), (2,1,"frieze",7)]
    println(sub)
    for kind in ["operations", "generators"]
        println("  ", kind)
        for num in 1:maxnum
            println("    ", num)
            gens_str = crawl_subperiodic_xyzt(num, D, P, kind)
            filename = (@__DIR__)*"/../test/data/xyzt-$kind/subperiodic/$sub/"*string(num)*".csv"
            open(filename; write=true, create=true, truncate=true) do io
                first = true
                for str in gens_str
                    # skip identity operation unless it's the only element; otherwise redundant
                    if kind == "operations" ||
                       (((D == 3 && str != "x,y,z") || (D == 2 && str != "x,y")) ||
                        length(gens_str) == 1)
                        if first
                            first = false
                        else
                            write(io, '\n')
                        end
                        write(io, str)
                    end
                end
            end
        end
    end
end