using HTTP, Gumbo, Cascadia, Crystalline, LinearAlgebra

function crawl_bandpaths(sgnum::Integer, timereversal::Bool=true; maxfails=10)
    body = crawl_bandpaths_html(sgnum, timereversal; maxfails)
    return extract_bandpaths_tables(body)
end

function crawl_bandpaths_html(sgnum::Integer, timereversal::Bool=true; maxfails=10)
    brtype = timereversal ? "Elementary TR" : "Elementary"
    brtypename = lowercasefirst(filter(!isspace, brtype))
    
    inputs = "super=$(sgnum)&"*                    # inputs to <form> ... </form>, in POST mode; 
             "$(brtypename)=$(brtype)&"*           # inputs are supplied as name=value pairs, 
             "tr=" * (timereversal ? "yes" : "no") # with distinct inputs separated by '&' 
                                                   # (https://stackoverflow.com/a/8430062/9911781)

    response = retrying_http_request(inputs, maxfails)

    return parsehtml(String(response.body))
end

function retrying_http_request(inputs, maxfails=10)
    for _ in 1:maxfails
        try 
            return HTTP.post("https://www.cryst.ehu.es/cgi-bin/cryst/programs/bandpaths.pl",
                             [], inputs)
        catch
        end
    end
    error("failed to retrieve data over $maxfails attempts")
end

function extract_bandpaths_tables(body::HTMLDocument)
    # return a 2-tuple `ts` of "tables":
    #  `ts[1]`: the k-point connectivity table
    #  `ts[2]`: the associated irrep compatibility/subduction data
    root=body.root
    # we need the 2nd and 3rd table-environments of `body`
    select = sel"table"
    matches = eachmatch(select, root)
    bandpath_htmltable, subduction_htmltable = matches[2], matches[3]

    return convert_html_table(bandpath_htmltable), convert_html_table(subduction_htmltable)
end

function convert_html_table(t::HTMLElement{:table})
    table_rows = Vector{Vector{String}}[]
    for (i, row) in enumerate(eachmatch(sel"tr", t)) # over rows
        i == 1 && continue # skip header row
        row_elements = Vector{String}[]
        for element in eachmatch(sel"td", row)       # over columns in row
            push!(row_elements, element_to_strings(element))
        end
        push!(table_rows, row_elements)
    end

    return permutedims(reduce(hcat, table_rows))
end

function element_to_strings(e::HTMLElement)
    # read the text in `e`, and if that text contains any <br/> tags, split on those tags
    # and return a vector of strings; if there are no <br/> tags, just return a 1-element
    # vector with the text content of `e`
    strs = String[]
    io = IOBuffer()
    for c in e.children
        if c isa HTMLElement{:br}
            push!(strs, String(take!(io)))
        else
            print(io, c)
        end
    end
    push!(strs, String(take!(io)))

    return html2unicode_postprocess.(strs)
end

function html2unicode_postprocess(s::AbstractString)
    replacepairlist = (
        r"\<sup\>(.*?)\<\/sup\>" => x->Crystalline.supscriptify(x[6:end-6]),     # superscripts (convert to unicode)
        r"\<sub\>(.*?)\<\/sub\>" => x->Crystalline.subscriptify(x[6:end-6]),     # subscripts   (convert to unicode)
        r"\<font style\=\"text-decoration:overline;\"\>(.*?)\<\/font\>"=>s"\1ˢ", # spinful irrep
        r"\([0-9]\)" => "", # get rid of irrep multiplicity information
        "⊕"=>"+",          # convert annoying direct add symbols to plain add
        "GP"=>"Ω",          # we denote the general point Ω, not GP
        isspace=>"")        # remove extraneous space

    return replace(replace(s, replacepairlist...), " "=>"")
end


# ---------------------------------------------------------------------------------------- #
# Types

# defines LabeledKVec{D}, Connection{D}, SubductionTable{D}
#include("../BandGraphs/src/subduction-types.jl")
using BandGraphs: LabeledKVec, Connection, SubductionTable

# ---------------------------------------------------------------------------------------- #
# Connections between k-vectors

function parse_labeledkpoint(s::AbstractString)
    label, kstr = split(s, ":")[1:2]
    return LabeledKVec(Symbol(label), KVec{3}(kstr))
end
function parse_connections(bandpaths_data::Matrix{<:AbstractString})
    connections = Vector{Connection{3}}(undef, size(bandpaths_data, 1)*2)
    for (i, row) = enumerate(eachrow(bandpaths_data))
        # left to center
        kᴳ¹ = parse_labeledkpoint(row[1]) # left-hand maximal k-point in table
        kᴴ  = parse_labeledkpoint(row[2]) # non-maximal k-point in table
        kᴳ² = parse_labeledkpoint(row[3]) # right-hand maximal k-point in table
        connections[2i-1] = Connection(kᴳ¹, kᴴ)
        connections[2i]   = Connection(kᴳ², kᴴ)
    end
    
    return unique!(connections)
end

function parse_connections(bandpaths_data::Matrix{<:AbstractVector{<:AbstractString}})
    return parse_connections(only.(bandpaths_data))
end

# ---------------------------------------------------------------------------------------- #
# function to find element of little group with highest-order screw (or, failing any screws,
# any glide operation); need this to figure out the k-point of a monodromy-related shifted
# k-point

#  function find_highest_order_nonsymmorphic_operation(sgnum::Integer)
#      issymmorph(sgnum, 3) && return nothing
#      sg = spacegroup(sgnum, Val(3))
#      N = round(Int, inv(det(Bravais.primitivebasismatrix(centering(sgnum, 3)))))
#      sg = sg[1:div(length(sg), N)] # only include the first "non-centering copies" operations
#  
#      return find_highest_order_nonsymmorphic_operation(sg, centering(sgnum, 3))
#  end
function find_highest_order_nonsymmorphic_operation(g::Crystalline.AbstractGroup{3}, cntr)
    ns = Crystalline.rotation_order.(g)
    ops = iterated_composition.(g, abs.(ns))
    for (j,op) in enumerate(ops)
        if !(rotation(op) ≈ I) && !(rotation(op) ≈ -I)
            error("unexpectedly did not get identity or inversion operation on iteration")
        end
    end
    ts = translation.(ops)
    uns = sort!(unique(ns); rev=true)
    for n in uns
        idxs = something(findall(==(n), ns))
        for i in idxs
            t = ts[i]
            if norm(t) > Crystalline.DEFAULT_ATOL
                t′ = float.(rationalize.(t, tol=1e-2)) # hack to get rid of floating point errors
                if !all(isinteger, primitivize(RVec(t′), centering(g)).cnst)
                    error("obtained non-integer translation $t for operation $(g[i])")
                end
                return (g[i], t′)
            end
        end
    end

    error("did not find any nonsymmorphic operations")
end

function iterated_composition(op, n)
    n < 0 && error("negative n not handled")
    n == 0 && return one(op)
    op′ = op
    for _ in 1:n-1
        op′ = compose(op′, op, false)
    end
    return op′
end

# ---------------------------------------------------------------------------------------- #
# Irrep subduction

function parse_subductions(
            subduction_data::Matrix{<:AbstractVector{<:AbstractString}},
            num::Integer;
            spinful::Bool = false)

    spinful_check = spinful ? contains("ˢ") : !contains("ˢ")
    
    Ncols = size(subduction_data, 2)
    kᴴ_colidx = Ncols == 5 ? 3 :
                Ncols == 7 ? 4 :
                error("unhandled subduction_data table dimensions")
    kᴳ_colidxs = (1, Ncols)
    rules_colidxs = (2, Ncols-1)

    tables = Vector{SubductionTable{3}}()
    for row in eachrow(subduction_data)
        for (kᴳ_colidx, rules_colidx) in zip(kᴳ_colidxs, rules_colidxs)
            c, irlabsᴳ, irlabsᴴ, table = 
                something(parse_subductions_of_columns_in_row(
                    row, kᴳ_colidx, kᴴ_colidx, rules_colidx, spinful_check))
            push!(tables, SubductionTable(num, c, irlabsᴳ, irlabsᴴ, table, false))
        end
    end

    if Ncols == 7 # monodromy additions
        lgs = littlegroups(num, 3)
        for row in eachrow(subduction_data)
            for (kᴳ_colidx, rules_colidx) in zip(kᴳ_colidxs, (3, Ncols-2))
                tmp = parse_subductions_of_columns_in_row(
                    row, kᴳ_colidx, kᴴ_colidx, rules_colidx, spinful_check)
                if !isnothing(tmp)
                    # we have a set of monodromy-derived subduction rules to complement the
                    # plain rules
                    c, irlabsᴳ, irlabsᴴ, table = tmp

                    # now, we want to infer the k-point associated with the translated high-
                    # symmetry k-point `kᴳ′` and the associated connection `c′`
                    kᴳ, kᴴ = c.kᴳ, c.kᴴ
                    
                    kᴴlab = kᴴ.label
                    lg = lgs[string(kᴴlab)]
                    # FIXME/TODO: the below approach to get the "translated" k-coordinate
                    #             doesn't seem right: need to fix eventually, but not really
                    #             important per se
                    t = find_highest_order_nonsymmorphic_operation(lg, centering(num, 3))[2]
                    kᴳ′ = LabeledKVec(Symbol(kᴳ.label, '′'), kᴳ.kv + t)
                    c′ = Connection(kᴳ′, kᴴ)
                    irlabsᴳ′ = map(irlab->replace(irlab, string(kᴳ.label) => string(kᴳ′.label)), irlabsᴳ)
                    push!(tables, SubductionTable(num, c′, irlabsᴳ′, irlabsᴴ, table, true))
                end
            end
        end
    end


    return unique!(t->t.c, tables)
end

function parse_subductions_of_columns_in_row(
            row, kᴳ_colidx, kᴴ_colidx, rules_colidx, spinful_check)
    kᴳ = parse_labeledkpoint(only(row[kᴳ_colidx]))
    kᴴ = parse_labeledkpoint(only(row[kᴴ_colidx]))
    c = Connection(kᴳ, kᴴ)
    subduction_rules = row[rules_colidx]

    # check if we are in a "monodromy" column; there may be no rules then
    length(subduction_rules) == 1 && isempty(subduction_rules[1]) && return nothing

    # okay: there are some rules, now we process them
    # step 1: find the irreps and sort them
    irlabsᴳ, irlabsᴴ = Vector{String}(), Vector{String}()
    for subduction_rule in subduction_rules
        irlabᴳ, irlabsᴴ_expr = split(subduction_rule, "→")
        spinful_check(irlabᴳ) && irlabᴳ ∉ irlabsᴳ && push!(irlabsᴳ, irlabᴳ)
        irlabsᴴ′ = split(replace(irlabsᴴ_expr, r"[0-9]"=>""), "+")
        for irlabᴴ in irlabsᴴ′
            spinful_check(irlabᴴ) && irlabᴴ ∉ irlabsᴴ && push!(irlabsᴴ, irlabᴴ)
        end
    end
    sort!(irlabsᴳ); sort!(irlabsᴴ)

    # step 2: in the above sorting, find the matrix representation of subduction rules
    table = zeros(Int, length(irlabsᴳ), length(irlabsᴴ))
    for subduction_rule in subduction_rules
        irlabᴳ, irlabsᴴ_expr = split(subduction_rule, "→")
        spinful_check(irlabᴳ) || continue
        i = something(findfirst(==(irlabᴳ), irlabsᴳ))
        irlabsᴴ_parts = split(irlabsᴴ_expr, "+")
        for irlabᴴ_part in irlabsᴴ_parts
            label_idx = findfirst(isletter, irlabᴴ_part)
            irlabᴴ = irlabᴴ_part[label_idx:end]
            spinful_check(irlabᴴ) || continue
            j = findfirst(==(irlabᴴ), irlabsᴴ)

            mult_idx = prevind(irlabᴴ_part, label_idx)               
            mult = mult_idx == 0 ? 1 : parse(Int, irlabᴴ_part[1:mult_idx])
            
            table[i,j] = mult
        end
    end

    return c, irlabsᴳ, irlabsᴴ, table
end

# ---------------------------------------------------------------------------------------- #

connectionsd = Dict{Int, Vector{Connection{3}}}()
subductionsd = Dict{Int, Vector{SubductionTable{3}}}()
data = Dict{Int, Any}()
for sgnum in 3:230
    # Note: SGs 1 & 2 fail since there are no connections
    println("sgnum ", sgnum)
    tmp = crawl_bandpaths(sgnum, true)
    data[sgnum] = tmp
    bandpath_data, subduction_data = tmp
    connectionsd[sgnum] = parse_connections(bandpath_data)
    subductionsd[sgnum] = parse_subductions(subduction_data, sgnum; spinful=false)
end
for sgnum in 1:2
    connectionsd[sgnum] = Connection{3}[]
    subductionsd[sgnum] = SubductionTable{3}[]
end

using JLD2
jldsave("BandGraphs/data/connections/3d/subductions-tr.jld2"; subductionsd=subductionsd)
jldsave("BandGraphs/data/connections/3d/connections.jld"; connectionsd=connectionsd)