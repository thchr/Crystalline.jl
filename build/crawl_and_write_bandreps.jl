using Crystalline, HTTP
import ProgressMeter: @showprogress

# crawling functionality
const BANDREP_URL="http://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl"
"""
    crawlbandreps(sgnum::Integer, allpaths::Bool=false, timereversal::Bool=true)
                                                       --> ::String (valid HTML table)
    
Crawls a band representation from the Bilbao database. This is achieved by sending a "POST"
request to the `<form> ... </form>` part of their HTML page.

**Input**:
- `sgnum`    : Space group number 
- `allpaths` : Whether to include only maximal **k**-points (true) or 
                all, i.e. general, **k**-points (true)
- `timereversal` : Whether the band representation assumes time-reversal symmetry (`true`)
                   or not (`false`)Type of bandrep, as string, either

Note that we do not presently try to crawl Wyckoff-type inputs.

**Output**: 
- `table_html` : Valid HTML <table> ... </table> for the requested band representation.
"""
function crawlbandreps(sgnum::Integer, allpaths::Bool=false, timereversal::Bool=true)
    # with (true ⇒ "Elementary TR") or without (false ⇒ "Elementary") time-reversal symmetry
    brtype = timereversal ? "Elementary TR" : "Elementary"
    brtypename = lowercasefirst(filter(!isspace, brtype))
    
    inputs = "super=$(sgnum)&"*                     # inputs to <form> ... </form>, in POST mode; 
             "$(brtypename)=$(brtype)&"*            # inputs are supplied as name=value pairs, 
             "nomaximal=$(allpaths ? "yes" : "no")" # with distinct inputs separated by '&' 
                                                    # (see e.g. https://stackoverflow.com/a/8430062/9911781)
    response = HTTP.post(BANDREP_URL, [], inputs)
    body = String(response.body)
    startstring = r"\<tr\>\<td(?:.*?)\>Wyckoff pos."
    startidx=findfirst(startstring, body)[1]
    stopidx=findlast("</td></tr></table>",body)[end]
    table_html="<table>"*body[startidx:stopidx]

    return table_html
end

# parsing functionality
function html2dlm(body::String, oplus::Union{String,Char}='⊕')
    dlm='|'
    # list of replacements, mostly using regexes; see https://regexr.com for inspiration.
    replacepairlist= (
        r"\<tr\>(.*?)\<\/tr\>" => SubstitutionString("\\1\n"), # rows; newline-separated (requires ugly hack around https://github.com/JuliaLang/julia/issues/27125 to escape \n)
        "<tr>"=>"\n",                           # <tr> tags are not always paired with </tr>; second pass here to remove stragglers
        r"\<td(?:.*?)\>(.*?)\<\/td\>" => SubstitutionString("\\1$(dlm)"), # columns; comma-separated
        r"\<sup\>(.*?)\<\/sup\>" => x->Crystalline.supscriptify(x[6:end-6]), # superscripts (convert to unicode)
        r"\<sub\>(.*?)\<\/sub\>" => x->Crystalline.subscriptify(x[6:end-6]), # subscripts   (convert to unicode)
        r"\<i\>(.*?)\<\/i\>"=>s"\1",            # italic annotations
        r"\<center\>(.*?)\<\/center\>"=>s"\1",  # centering annotations
        "<br>"=>"",                             # linebreak tag in html; no </br> tag exists
        "<font size=\"5\">&uarr;</font>"=>"↑",  # induction arrow
        r"\<font style\=\"text-decoration:overline;\"\>(.*?)\<\/font\>"=>s"\1ˢ", # spinful irrep
        "&oplus;"=>oplus,                       # special symbols # ⊕
        "&Gamma;"=>'Γ',                                           # Γ
        "&Sigma;"=>'Σ',                                           # Σ
        "&Lambda;"=>'Λ',                                          # Λ
        "&Delta;"=>'Δ',                                           # Δ
        "GP"=>'Ω',                                                # GP to Ω
        "&nbsp;"=>"",                                             # non-breaking space
        '*'=>"",                                                  # * (for Σ k-points in e.g. sg 146: sending `u,-2*u,0`⇒`u,-2u,0`)
        r"\<form.*?\"Decomposable\"\>\<\/form\>"=>"Decomposable", # decomposable bandreps contain a link to decompositions; ignore it
        "$(dlm)Decomposable"=>"$(dlm)true",    # simplify statement of composability
        "$(dlm)Indecomposable"=>"$(dlm)false",
        "\\Indecomposable"=>"",
        "<table>"=>"","</table>"=>"",           # get rid of table tags
        "$(dlm)\n"=>"\n",                       # if the last bits of a line is ", ", get rid of it
        r"\n\s*\Z"=>"",                         # if the last char in the string is a newline (possibly with spurious spaces), get rid of it
        r"\n\s+"=>"\n",                         # remove spurious white/space at start of any lines
        #r"((_.){2,})"=>(x)->"_{"*replace(x,"_"=>"")*"}",  # tidy up multi-index subscripts
        #r"((\^.){2,})"=>(x)->"^{"*replace(x,"^"=>"")*"}"  # tidy up multi-index superscripts
    )

    for replacepair in replacepairlist
        body = replace(body, replacepair)
    end

    return body
end

#=
# utilities for conversion between different textual/array/struct representations
html2array(body) = Crystalline.dlm2array(html2dlm(body))

function html2struct(body::String, sgnum::Integer, allpaths::Bool=false, 
                     spinful::Bool=false, timereversal::Bool=true)
    Crystalline.array2struct(html2array(body), sgnum, allpaths, spinful, timereversal)
end
=#

# crawl & write bandreps to disk
function writebandreps(sgnum, allpaths, timereversal=true)
    paths_str = allpaths ? "allpaths" : "maxpaths"
    brtype_str = timereversal ? "ElementaryTR" : "Elementary"

    BR_dlm = html2dlm(crawlbandreps(sgnum, allpaths, timereversal), '⊕')

    filename = (@__DIR__)*"/../data/bandreps/3d/$(brtype_str)/$(paths_str)/$(string(sgnum)).csv"
    open(filename; write=true, create=true, truncate=true) do io
        write(io, BR_dlm)
    end 
end

# run to crawl everything... (takes ∼5 min)
for allpaths in [true, false]
    for timereversal in [true, false]
        @info "allpaths=$allpaths, timereversal=$timereversal"
        @showprogress 0.1 "Crawling ..." for sgnum in 1:MAX_SGNUM[3]
            writebandreps(sgnum, allpaths, timereversal)
        end
        # @sync for sgnum in 1:MAX_SGNUM[3] # spuriously fails on HTTP requests; not worth it
        #     @async writebandreps(sgnum, allpaths, timereversal)
        # end
    end
end