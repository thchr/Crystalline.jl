using Crystalline, HTTP, Gumbo

# Small convenience script to crawl and subsequently write all the xyzt forms of the
# generators of the 3D space-groups.

# ---------------------------------------------------------------------------------------- #
## Functions for crawling the generators of 3D space groups from Bilbao

""" 
    crawl_generators_xyzt(pgnum::Integer, D::Integer=3)

Obtains the generators in xyzt format for a given point group number `pgnum` by crawling
the Bilbao server. Only works for `D = 3`.
"""
function crawl_pg_generators_xyzt(pgnum::Integer, D::Integer=3, 
                                  include_pgiuc::Union{Nothing, String}=nothing)
    htmlraw = crawl_pg_generators_html(pgnum, D, include_pgiuc)

    # grab relevant rows of table
    table_html = children(children(children(children(last(children(htmlraw.root)))[7])[1])[1])[3:end]
    # grab relevant columns of table & extract xyzt-forms of generators
    gens_str = strip.(text.(only.(children.(getindex.(table_html, 2)))))
    
    return gens_str
end

function crawl_pg_generators_html(pgnum::Integer, D::Integer=3,
                                  include_pgiuc::Union{Nothing, String}=nothing)
    baseurl = begin
        if D == 3
            "https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-point_genpos?w2do=gens&num="
        elseif D == 2
            "https://www.cryst.ehu.es/cgi-bin/plane/programs/nph-point_plane-genpos?w2do=gens&num="
            include_pgiuc === nothing || error(DomainError(setting, "setting requests cannot be made separately in 2D"))
        else
            error(DomainError(D, "only 2D and 3D point group generators are crawlable"))
        end
    end
    url = baseurl * string(pgnum)
    if include_pgiuc !== nothing
        url *= "&what=" * include_pgiuc
    end
    contents = HTTP.request("GET", url)
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #
## Use crawling functions & write information to `/data/generators/pgs/`
let D = 3
    for (pgnum, pgiucs) in enumerate(Crystalline.PG_NUM2IUC[D])
        for (setting, pgiuc) in enumerate(pgiucs)
            if pgiuc ∈ ("-62m", "-6m2")
                # need to flip the meaning of `setting` here to match sorting on Bilbao
                setting = setting == 1 ? 2 : (setting == 2 ? 1 : error())
            end
            if setting == 1
                gens_str = crawl_pg_generators_xyzt(pgnum, D)
            else
                gens_str = crawl_pg_generators_xyzt(pgnum, D, pgiuc)
            end

            unmangled_pgiuc = Crystalline.unmangle_pgiuclab(pgiuc)
            filename = (@__DIR__)*"/../data/generators/pgs/$(D)d/$(unmangled_pgiuc).csv"
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
    end
end

# Dimension 2 is following a different scheme than dimension 3 on Bilbao, so just correct
# for that manually
let D = 2
    for (pgnum, pgiuc) in enumerate(Crystalline.PG_IUCs[D])
        gens_str = crawl_pg_generators_xyzt(pgnum, D)
        
        unmangled_pgiuc = Crystalline.unmangle_pgiuclab(pgiuc)
        filename = (@__DIR__)*"/../data/generators/pgs/$(D)d/$(unmangled_pgiuc).csv"
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
end