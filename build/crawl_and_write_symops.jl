using Crystalline, HTTP, Gumbo

#= Small convenience script to crawl and subsequently write all the xyzt
   forms of the symmetry operations of the 230 three-dimensional space-
   groups. Enables us to just read the symmetry data from the hard-disk
   rather than constantly querying the Bilbao server =#
for sgnum in 1:230
    sgops_str = Crystalline.crawl_sgops_xyzt(sgnum)
    filename = (@__DIR__)*"/../data/sgops/3d/"*string(sgnum)*".csv"
    open(filename; write=true, create=true, truncate=true) do io
        foreach(str -> write(io, str, '\n'), sgops_str)
    end 
end

# ----- NOW-REDUNANT FUNCTIONS FOR CRAWLING 3D SPACE GROUPS FROM BILBAO -----
""" 
    crawl_sgops_xyzt(sgnum::Integer, D::Integer=3)

Obtains the symmetry operations in xyzt format for a given space group
number `sgnum` by crawling the Bilbao server; see `spacegroup` for 
additional details. Only works for `D = 3`.
"""
function crawl_sgops_xyzt(sgnum::Integer, D::Integer=3)
    htmlraw = crawl_sgops_html(sgnum, D)

    ops_html = children.(children(last(children(htmlraw.root)))[4:2:end])
    Nops = length(ops_html)
    sgops_str = Vector{String}(undef,Nops)

    for (i,op_html) in enumerate(ops_html)
        sgops_str[i] = _stripnum(op_html[1].text) # strip away the space group number
    end
    return sgops_str
end

function crawl_sgops_html(sgnum::Integer, D::Integer=3)
    D != 3 && error("only 3D space group operations are crawlable")
    (sgnum < 1 || sgnum > 230) && error(DomainError(sgnum))

    baseurl = "http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?what=text&gnum="
    contents = HTTP.request("GET", baseurl * string(sgnum))
    return parsehtml(String(contents.body))
end

function _stripnum(s)
    if occursin(' ', s) # if the operation "number" is included as part of s
        _,s′ = split(s, isspace; limit=2)
    end
    return String(s′) # ensure we return a String, rather than possibly a SubString
end