struct BandRepTrait end

# crawling functionality
const BANDREP_URL="http://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl"
"""
    crawlbandreps(sgnum::Integer, allpaths::Bool=false, brtype::String="Elementary TR")
                                                       --> ::String (valid HTML table)
    
    Crawls a band representation (BR) from the Bilbao database. This is achieved
    by sending a 'POST' request to the <form> ... </form> part of their HTML page.

    Input:
    - sgnum    : Space groupnumber 
    - allpaths : Whether to include only maximal k-points (true) or 
                 all, i.e. general, k-points (true)
    - brtype   : Type of bandrep, as string, either
                 ∘ "Elementary TR" : Elementary BR, with time-reversal assumed
                 ∘ "Elementary"    : Elementary BR, without time-reversal symmetry
                 (we do not presently try to crawl Wyckoff-type inputs)

    Output     : Valid HTML <table> ... </table> for the requested BR.
"""
function crawlbandreps(sgnum::Integer, allpaths::Bool=false, brtype::String="Elementary TR")
    if brtype != "Elementary TR" && brtype != "Elementary"
        error(ArgumentError("String 'brtype' must be 'Elementary TR' or 'Elementary', was '$(brtype)'"))
    else
        brtypename = lowercasefirst(filter(!isspace, brtype))
    end
    
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
function html2dlm(body::String, oplus::Union{String,Char}="⊕ ", ::BandRepTrait=BandRepTrait())
    dlm='|'
    # list of replacements, mostly using regexes; see https://regexr.com for inspiration.
    replacepairlist= (
        r"\<tr\>(.*?)\<\/tr\>" => SubstitutionString("\\1\n"), # rows; newline-separated (requires ugly hack around https://github.com/JuliaLang/julia/issues/27125 to escape \n)
        "<tr>"=>"\n",                           # <tr> tags are not always paired with </tr>; second pass here to remove stragglers
        r"\<td(?:.*?)\>(.*?)\<\/td\>" => SubstitutionString("\\1$(dlm)"), # columns; comma-separated
        r"\<sup\>(.*?)\<\/sup\>" => x->supscriptify(x[6:end-6]), # superscripts (convert to unicode)
        r"\<sub\>(.*?)\<\/sub\>" => x->subscriptify(x[6:end-6]), # subscripts   (convert to unicode)
        r"\<i\>(.*?)\<\/i\>"=>s"\1",            # italic annotations
        r"\<center\>(.*?)\<\/center\>"=>s"\1",  # centering annotations
        "<br>"=>"",                             # linebreak tag in html; no </br> tag exists
        "<font size=\"5\">&uarr;</font>"=>"↑",  # induction arrow
        r"\<font style\=\"text-decoration:overline;\"\>(.*?)\<\/font\>"=>s"\1ˢ", # spinful irrep
        "&oplus;"=>oplus,                       # special symbols # ⊕
        "&Gamma;"=>"Γ",                                           # Γ
        "&Sigma;"=>"Σ",                                           # Σ
        "&Lambda;"=>"Λ",                                          # Λ
        "&Delta;"=>"Δ",                                           # Δ
        "&nbsp;"=>"",                                             # non-breaking space
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

# utilities for conversion between different textual/array/struct representations
function html2array(body)
    str = html2dlm(body, BandRepTrait()) 
    M = dlm2array(str);
end

function html2struct(body::String, sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, 
                     ::BandRepTrait=BandRepTrait())
    M = html2array(body)
    array2struct(M, sgnum, allpaths, spinful, BandRepTrait())
end

dlm2array(str::String) = DelimitedFiles.readdlm(IOBuffer(str), '|', String, '\n')
dlm2array(io::IO) = DelimitedFiles.readdlm(io, '|', String, '\n')

# utilities for creation of BandRepStruct 
function dlm2struct(str::Union{String,IO}, sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, 
                    ::BandRepTrait=BandRepTrait())
    M = dlm2array(str);
    array2struct(M, sgnum, allpaths, spinful, BandRepTrait())
end

function array2struct(M::Matrix{String}, sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, 
                      ::BandRepTrait=BandRepTrait()) 

    klist =  permutedims(mapreduce(x->String.(split(x,":")), hcat, M[4:end,1])) # 1ˢᵗ col is labels, 2ⁿᵈ col is coordinates as strings
    klabs, kvs = (@view klist[:,1]), KVec.(@view klist[:,2])

    temp = split_paren.(@view M[1,2:end]) 
    wyckpos, sitesym = getindex.(temp, 1), getindex.(temp, 2) # we really shouldn't be doing this type of trickery ...

    temp .= split_paren.(@view M[2,2:end]) # same size, so reuse array
    label, dim = getindex.(temp, 1), parse.(Int64, getindex.(temp, 2))

    decomposable = parse.(Bool, vec(@view M[3,2:end]))

    irreptags = collect(eachcol(@view M[4:end, 2:end])) 
    for col in irreptags 
        col .= replace.(col, Ref(r"\([1-9]\)"=>""))  # get rid of irrep dimension info
    end

    BRs = BandRep.(wyckpos, sitesym, label, dim, decomposable, 
                   map(isspinful, irreptags), irreptags)
    if !spinful 
        filter!(x->!x.spinful, BRs)        
    end
    return BandRepSet(sgnum, BRs, kvs, klabs, allpaths, spinful)
end

# main "getter" function; reads data from csv files
function bandreps(sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, brtype::String="Elementary TR")
    paths_str = allpaths ? "allpaths" : "maxpaths"
    filename = (@__DIR__)*"/../data/bandreps/3d/$(filter(!isspace, brtype))/$(paths_str)/$(string(sgnum)).csv"
    open(filename) do io
        BRS = dlm2struct(io, sgnum, allpaths, spinful, BandRepTrait())
    end 
end

# Convert a BandRepSet into a matrix representation
#= 
function matrix(BRS::BandRepSet)
    unique(reps(BRS)))
    nbrs = length(BRS)
    A = Matrix{Int64}(nbrs)
end 
=#

# misc minor utility functions
function split_paren(str::AbstractString)
    openpar = findfirst(x->x=='(', str) # index of the opening parenthesis
    before_paren = SubString(str, firstindex(str), prevind(str, openpar))
    inside_paren = SubString(str, nextind(str, openpar), prevind(str, lastindex(str)))
    return before_paren, inside_paren
end

isspinful(col::AbstractVector) = any(x->occursin(r"\\bar|ˢ", x), col)

function searchpriornumerals(coord, pos₂)
    pos₁ = copy(pos₂)
    while (prev = prevind(coord, pos₁)) != 0 # not first character
        if isnumeric(coord[prev]) || coord[prev] == '.'
            pos₁ = prev
        elseif coord[prev] == '+' || coord[prev] == '-'
            pos₁ = prev
            break
        else 
            break
        end
    end
    prefix = coord[pos₁:prevind(coord, pos₂)] 
    if !any(isnumeric, prefix) # in case there's no numerical prefix, it must be unity
        prefix *= "1"
    end

    return prefix
end
