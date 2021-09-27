using Crystalline, HTTP, Gumbo

# Small convenience script to crawl and subsequently write all the xyzt forms of the
# generators of the 3D space-groups.

# ---------------------------------------------------------------------------------------- #
## Functions for crawling the generators of 3D space groups from Bilbao

""" 
    crawl_generators_xyzt(sgnum::Integer, D::Integer=3)

Obtains the generators in xyzt format for a given space group number `sgnum` by crawling
the Bilbao server. Only works for `D = 3`.
"""
function crawl_generators_xyzt(sgnum::Integer, D::Integer=3)
    htmlraw = crawl_generators_html(sgnum, D)

    # grab relevant rows of table
    table_html = children(children(children(children(last(children(htmlraw.root)))[6])[1])[1])[3:end]
    # grab relevant columns of table & extract xyzt-forms of generators
    gens_str = strip.(text.(only.(children.(getindex.(table_html, 2)))))
    
    return gens_str
end

function crawl_generators_html(sgnum::Integer, D::Integer=3)
    D != 3 && error("only 3D space group operations are crawlable")
    (sgnum < 1 || sgnum > 230) && error(DomainError(sgnum))

    baseurl = "http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?what=&gnum="
    contents = HTTP.request("GET", baseurl * string(sgnum))
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #
## Use crawling functions & write information to `/data/`
for sgnum in 1:230
    println(sgnum)
    gens_str = crawl_generators_xyzt(sgnum)
    filename = (@__DIR__)*"/../data/generators/sgs/3d/"*string(sgnum)*".csv"
    open(filename; write=true, create=true, truncate=true) do io
        first = true
        for str in gens_str
            # skip identity operation unless it's the only element; otherwise redundant
            if str != "x,y,z" || length(gens_str) == 1
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