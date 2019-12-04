struct BandRepTrait end

# crawling functionality
const BANDREP_URL="http://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl"
"""
    crawlbandreps(sgnum::Integer, allpaths::Bool=false, brtype::String="Elementary TR")
                                                       --> ::String (valid HTML table)
    
Crawls a band representation (BR) from the Bilbao database. This is achieved
by sending a 'POST' request to the <form> ... </form> part of their HTML page.

**Input**:
- `sgnum`    : Space group number 
- `allpaths` : Whether to include only maximal **k**-points (true) or 
                all, i.e. general, **k**-points (true)
- `brtype`   : Type of bandrep, as string, either
    - "Elementary TR" : Elementary BR, with time-reversal assumed
    - "Elementary"    : Elementary BR, without time-reversal symmetry 
    (we do not presently try to crawl Wyckoff-type inputs)

**Output**: Valid HTML <table> ... </table> for the requested BR.
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
function html2dlm(body::String, oplus::Union{String,Char}='⊕', ::BandRepTrait=BandRepTrait())
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
    wyckpos, sitesym = getindex.(temp, 1), getindex.(temp, 2) # wyckoff position and site symmetry point group of bandrep

    temp .= split_paren.(@view M[2,2:end]) # same size, so reuse array
    label, dim = getindex.(temp, 1), parse.(Int64, getindex.(temp, 2)) # label of bandrep

    decomposable = parse.(Bool, vec(@view M[3,2:end])) # whether bandrep can be further decomposed

    brtags = collect(eachcol(@view M[4:end, 2:end])) # set of irreps that jointly make up the bandrep
    for br in brtags 
        br .= replace.(br, Ref(r"\([1-9]\)"=>""))  # get rid of irrep dimension info
    end
    if !spinful 
        delidxs = findall(map(isspinful, brtags))
        for vars in (brtags, wyckpos, sitesym, label, dim, decomposable)
            deleteat!(vars, delidxs)    
        end
    end
    irreplabs, irrepvecs = get_irrepvecs(brtags)              

    BRs = BandRep.(wyckpos, sitesym, label, dim, decomposable, 
                   map(isspinful, brtags), irrepvecs, brtags)

    return BandRepSet(sgnum, BRs, kvs, klabs, irreplabs, allpaths, spinful)
end


function get_irrepvecs(brtags)
    Nklabs = length(first(brtags)) # there's equally many (composite) irrep tags in each band representation
    irreplabs = Vector{String}()
    for kidx in Base.OneTo(Nklabs)
        irreplabs_at_kidx = Vector{String}()
        for tag in getindex.(brtags, kidx) # tag could be a combination like Γ1⊕2Γ₂ (or something simpler, like Γ₁)
            for irrep in split(tag, '⊕')
                irrep′ = filter(!isdigit, irrep) # filter off any multiplicities
                if irrep′ ∉ irreplabs_at_kidx
                    push!(irreplabs_at_kidx, irrep′)
                end
            end
        end
        sort!(irreplabs_at_kidx)
        append!(irreplabs, irreplabs_at_kidx)
    end

    irrepvecs = [zeros(Int64, length(irreplabs)) for _=Base.OneTo(length(brtags))]
    for (bridx, tags) in enumerate(brtags)
        for (kidx,tag) in enumerate(tags)
            for irrep in split(tag, '⊕') # note this irrep tag may contain numerical prefactors!
                buf = IOBuffer(irrep)
                prefac_str = readuntil(buf, !isdigit)
                seek(buf, ncodeunits(prefac_str)) # go back to first non-digit position in buffer
                if isempty(prefac_str)
                    prefac = Int64(1)
                else
                    prefac = parse(Int64, prefac_str)
                end
                irrep′ = read(buf, String) # the rest of the irrep buffer is the actual cdml label
                close(buf)
                irrepidx = findfirst(==(irrep′), irreplabs) # find position in irreplabs vector
                irrepvecs[bridx][irrepidx] = prefac
            end
        end
    end
    return irreplabs, irrepvecs
end


# main "getter" function; reads data from csv files
function bandreps(sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, brtype::String="Elementary TR")
    paths_str = allpaths ? "allpaths" : "maxpaths"
    filename = (@__DIR__)*"/../data/bandreps/3d/$(filter(!isspace, brtype))/$(paths_str)/$(string(sgnum)).csv"
    open(filename) do io
        BRS = dlm2struct(io, sgnum, allpaths, spinful, BandRepTrait())
    end 
end



# Converts a BandRepSet into a matrix representation, with distinct band reps
# along the rows, and irreducible representations along the columns
function matrix(BRS::BandRepSet)
    M = Matrix{Int64}(undef, length(BRS), length(irreplabels(BRS)))
    @inbounds for (i,BR) in enumerate(reps(BRS))
        for (j,v) in enumerate(vec(BR))
            M[i,j] = v
        end
    end
    return M
end 


"""
    classification(BRS::BandRepSet) --> String

Calculate the symmetry indicator classification of a band representation 
set, meaning the index-classification inferrable on the basis of symmetry
alone.

Technically, the calculation answers a question like "what direct product 
of Zₙ groups is the the quotient group Xᵇˢ = {BS}/{AI} isomorphic to?".
See Po, Watanabe, & Vishwanath, Nature Commun. 8, 50 (2017).
"""
function classification(BRS::BandRepSet)
    Λ = smith(matrix(BRS)).SNF # get the diagonal components of the Smith normal decomposition
    Λ .= abs.(Λ) # the sign has no significance; can be absorbed in T or S if M = T⁻¹ΛS
    nontriv_idx = findall(x-> !(isone(x) || iszero(x)), Λ)
    if isempty(nontriv_idx)
        return "Z₁"
    else
        return ("Z"*join(subscriptify.(string.(sort(@view Λ[nontriv_idx]))), "×Z"))
    end
end

"""
    basisdim(BRS::BandRepSet) --> Int64

Computes the dimension of the (linearly independent parts) of a 
band representation set. This is dᵇˢ = dᵃⁱ in the notation of 
Po, Watanabe, & Vishwanath, Nature Commun. 8, 50 (2017). In other words,
this is the number of linearly independent basis vectors that span the 
expansions of a band structure or atomic insulator viewed as symmetry-data.
""" 
function basisdim(BRS::BandRepSet)
    Λ = smith(matrix(BRS)).SNF
    nnz = sum(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix
    return nnz
end


"""
    wyckbasis(BRS::BandRepSet) --> Vector{Vector{Int64}}

Computes the (band representation) basis for bands generated by local
orbitals placed at the Wyckoff positions. Any band representation that
can be expanded on this basis with positive integer coefficients 
correspond to a trivial insulators (i.e. deformable to atomic limit).
Conversely, bands that cannot are topological, either fragily (some  
negative coefficients) or strongly (fractional coefficients).
"""
function wyckbasis(BRS::BandRepSet) 
    # Compute Smith normal form: for an n×m matrix A with integer elements,
    # find matrices S, diagm(Λ), and T (of size n×n, n×m, and m×m, respectively)
    # with integer elements such that A = S*diagm(Λ)*T. Λ is a vector 
    # [λ₁, λ₂, ..., λᵣ, 0, 0, ..., 0] with λⱼ divisible by λⱼ₊₁ and r ≤ min(n,m).
    # The matrices T and S have integer-valued pseudo-inverses.
    F = smith(matrix(BRS)) 
    S,S⁻¹,T,T⁻¹,Λ = F.S,F.Sinv, F.T, F.Tinv, F.SNF
    nnz = sum(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix
    return T[Base.OneTo(nnz), :], Λ, T⁻¹[Base.OneTo(nnz), :]
end

# misc minor utility functions
isspinful(br::AbstractVector{T} where T<:AbstractString) = any(x->occursin(r"\\bar|ˢ", x), br)

function split_paren(str::AbstractString)
    openpar = findfirst(==('('), str) # index of the opening parenthesis
    before_paren = SubString(str, firstindex(str), prevind(str, openpar))
    inside_paren = SubString(str, nextind(str, openpar), prevind(str, lastindex(str)))
    return before_paren, inside_paren
end

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
