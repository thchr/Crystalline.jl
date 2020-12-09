using Crystalline, HTTP, Gumbo
using StaticArrays, LinearAlgebra

# ---------------------------------------------------------------------------------------- #
# CRAWL 3D WYCKOFF POSITIONS FROM BILBAO (2D and 1D obtained manually...)

const BILBAO_URL = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/"
const WYCK_URL_BASE_3D = BILBAO_URL*"nph-normsets?from=wycksets&gnum="
wyck_url(sgnum) = WYCK_URL_BASE_3D*string(sgnum)

""" 
    crawl_wyckpos(sgnum::Integer, D::Integer=3)

Obtains the Wyckoff positions for a given space group number `sgnum` by crawling the Bilbao
server; see `WyckPos` for additional details. 

TODO: Only works for `D = 3` currently. 
TODO: Does not yet compute the associated site symmetry groups.
"""
function crawl_wyckpos(sgnum::Integer, ::Val{D}=Val{3}()) where D
    D== 3 || throw("not implemented")
    htmlraw = crawl_wyckpos_html(sgnum, D)

    wycks_html = children.(children(children(children(last(children(htmlraw.root)))[3])[5][1]))
    
    Nwyck = length(wycks_html)-1
    wycks = Vector{WyckPos{D}}(undef, Nwyck)

    for (i,el) in enumerate(@view wycks_html[2:end])
        letter, mult_str, sitesym_str, qv_str, _ = getfield.(first.(getfield.(el, Ref(:children))), Ref(:text))

        qv = RVec{D}(qv_str)
        mult = parse(Int, mult_str) 

        wycks[i] = WyckPos{D}(mult, only(letter), qv)
    end
    return wycks
end

function crawl_wyckpos_html(sgnum::Integer, D::Integer=3)
    if D != 3; error("We do not crawl plane group data; see json files instead; manually crawled.") end
    if sgnum < 1 || sgnum > 230; error(DomainError(sgnum)); end

    if D == 3
        contents = HTTP.request("GET", wyck_url(sgnum))
        return parsehtml(String(contents.body))
    else
        error("We did not yet implement 2D plane groups")
    end
end

function _write_wyckpos_3d(sgnum::Integer, D::Integer=3)
    wps = crawl_wyckpos(sgnum, D)
    open((@__DIR__)*"/../data/wyckpos/3d/"*string(sgnum)*".csv", "w+") do io
        for (idx, wp) in enumerate(wps)
            qstr = strip(string(qvec(wp)), ('[',']'))
            for repl in ('α'=>'x', 'β'=>'y', 'γ'=>'z', " "=>"")
                qstr = replace(qstr, repl)
            end
            print(io, wp.mult, '|', wp.letter, '|', qstr)
            idx ≠ length(wps) && println(io)
        end
    end
end

foreach(_write_wyckpos_3d, 1:230)
