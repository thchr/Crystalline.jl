using Crystalline, HTTP, Gumbo
using StaticArrays, LinearAlgebra

# ---------------------------------------------------------------------------------------- #
# CRAWL 3D WYCKOFF POSITIONS FROM BILBAO (2D and 1D obtained manually...)

const BILBAO_URL = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/"
const WYCK_URL_BASE_3D = BILBAO_URL*"nph-normsets?from=wycksets&gnum="
wyck_3d_url(sgnum) = WYCK_URL_BASE_3D*string(sgnum)

""" 
    crawl_wyckpos_3d(sgnum::Integer)

Obtains the 3D Wyckoff positions for a given space group number `sgnum` by crawling the
Bilbao Crystallographic Server; returns a vector of `WyckoffPosition{3}`.
"""
function crawl_wyckpos_3d(sgnum::Integer)
    htmlraw = crawl_wyckpos_3d_html(sgnum)

    wps_html = children.(children(children(children(last(children(htmlraw.root)))[3])[5][1]))
    
    N_wps = length(wps_html)-1
    wps = Vector{WyckoffPosition{3}}(undef, N_wps)

    for (i,el) in enumerate(@view wps_html[2:end])
        letter, mult_str, sitesym_str, rv_str, _ = 
                                getfield.(first.(getfield.(el, Ref(:children))), Ref(:text))

        rv = RVec{3}(rv_str)
        mult = parse(Int, mult_str) 

        wps[i] = WyckoffPosition{3}(mult, only(letter), rv)
    end
    return wps
end

function crawl_wyckpos_3d_html(sgnum::Integer)
    (sgnum < 1 || sgnum > 230) && error(DomainError(sgnum))

    contents = HTTP.request("GET", wyck_3d_url(sgnum))
    return parsehtml(String(contents.body))
end

# ---------------------------------------------------------------------------------------- #
# CRAWL & SAVE/WRITE 3D WYCKOFF POSITIONS TO `data/wyckpos/3d/...`

function _write_wyckpos_3d(sgnum::Integer)
    wps = crawl_wyckpos_3d(sgnum)
    open((@__DIR__)*"/../data/wyckpos/3d/"*string(sgnum)*".csv", "w+") do io
        for (idx, wp) in enumerate(wps)
            qstr = strip(string(parent(wp)), ('[',']'))
            for repl in ('α'=>'x', 'β'=>'y', 'γ'=>'z', " "=>"")
                qstr = replace(qstr, repl)
            end
            print(io, wp.mult, '|', wp.letter, '|', qstr)
            idx ≠ length(wps) && println(io)
        end
    end
end

# actually crawl and write 3D Wyckoff positions
foreach(1:MAX_SGNUM[3]) do sgnum
    println(sgnum)
    _write_wyckpos_3d(sgnum)
end