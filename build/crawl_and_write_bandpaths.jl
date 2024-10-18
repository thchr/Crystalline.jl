using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using HTTP, Gumbo, Cascadia, LinearAlgebra
using Crystalline, BandGraphs

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

    bandpath_data   = convert_html_table(bandpath_htmltable)
    subduction_data = convert_html_table(subduction_htmltable)
    canonicalize_subduction_data_structure!(subduction_data)

    return bandpath_data, subduction_data
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

function canonicalize_subduction_data_structure!(subduction_data::Matrix)
    # The structure of the `subduction_data` table, as returned from BCS, is a bit finnicky.
    # If there are no monodromy additions, the table has 5 columns, structured as:
    #  | kᴳ₁ | kᴳ₁-kᴴ rules | | kᴴ | kᴳ₂-kᴴ rules | kᴳ₂ |
    # If there _are_ monodromy additions, the table has 7 columns, structured as:
    #  | kᴳ₁ | kᴳ₁-kᴴ rules |  kᴳ₁′-kᴴ rules | kᴴ | (...)² | kᴳ₂ | 
    # where by kᴳ₁′ and kᴳ₂′ we denote a monodromy-derived k-vector and where (...)²
    # represents represents 2 columns of rules, with a variable structure that depends on
    # whether there is a monodromy rule or not (whose presence is only optional):
    #  A. Has monodromy rule:  (...)² = | kᴳ₂-kᴴ rules      | kᴳ₂′-kᴴ rules |
    #  B. No monodromy rule:   (...)² = | empty string ("") | kᴳ₂-kᴴ rules  |
    # This is annoying because the structure is variable and depends on the presence of
    # monodromy rules.
    # To fix this, we simply reorder the columns to ensure the fixed following structure:
    # | kᴳ₁ | kᴳ₁-kᴴ rules | kᴳ₁′-kᴴ rules | kᴴ | kᴳ₂′-kᴴ rules | kᴳ₂-kᴴ rules | kᴳ₂ |
    Ncols = size(subduction_data, 2)
    if Ncols == 7
        for i in 1:size(subduction_data, 1)
            tᵢ₅ = subduction_data[i, 5]
            if !isempty(first(tᵢ₅))
                # there is a monodromy rule in this row, so we must flip columns 5 & 6
                tᵢ₆ = subduction_data[i, 6]
                subduction_data[i, 5] = tᵢ₆
                subduction_data[i, 6] = tᵢ₅
            end
        end
    end

    # now with `subduction_data` is in a consistent, simple form and we can parse it in
    # `parse_subductions` without worrying about a variable location for the different data 
    return subduction_data
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
function find_highest_order_nonsymmorphic_operation(
            g::Crystalline.AbstractGroup{D}, cntr::Char=centering(g, D)) where D

    ns = Crystalline.rotation_order.(g)
    ops = iterated_composition.(g, abs.(ns))
    for (j,op) in enumerate(ops)
        if !(rotation(op) ≈ I) && !(rotation(op) ≈ -I)
            error("unexpectedly did not get identity or inversion operation on iteration")
        end
    end
    ts = Crystalline.translation.(ops)
    uns = sort!(unique(ns); rev=true)
    for n in uns
        idxs = something(findall(==(n), ns))
        for i in idxs
            t = ts[i]
            if norm(t) > Crystalline.DEFAULT_ATOL
                t′ = float.(rationalize.(t, tol=1e-2)) # hack to get rid of floating point errors
                if !all(isinteger, primitivize(RVec(t′), cntr).cnst)
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

    # First, we do the "standard" rules (not monodromy-derived)
    tables = Vector{SubductionTable{3}}()
    for row in eachrow(subduction_data)
        for (kᴳ_colidx, rules_colidx) in zip(kᴳ_colidxs, rules_colidxs)
            c, irlabsᴳ, irlabsᴴ, table = 
                something(parse_subductions_of_columns_in_row(
                    row, kᴳ_colidx, kᴴ_colidx, rules_colidx, spinful_check))
            push!(tables, SubductionTable(num, c, irlabsᴳ, irlabsᴴ, table, false))
        end
    end

    # Now we do the monodromy additions
    if Ncols == 7
        lgs = littlegroups(num, 3)
        for row in eachrow(subduction_data)
            for (kᴳ_colidx, rules_colidx) in zip(kᴳ_colidxs, (3, Ncols-2))
                tmp = parse_subductions_of_columns_in_row(
                    row, kᴳ_colidx, kᴴ_colidx, rules_colidx, spinful_check)
                if !isnothing(tmp)
                    # we have a set of monodromy-derived subduction rules to complement the
                    # plain rules
                    c, irlabsᴳ, irlabsᴴ, table = tmp

                    # determine the monodromy-related k-point and the associated connection
                    c′ = find_monodromy_related_connection(c, lgs)

                    irlabsᴳ′ = map(irlabsᴳ) do irlab
                        replace(irlab, string(c.kᴳ.label) => string(c′.kᴳ.label))
                    end
                    push!(tables, SubductionTable(num, c′, irlabsᴳ′, irlabsᴴ, table, true))
                end
            end
        end
    end

    return unique!(t->t.c, tables)
end

# infer the k-point associated with the shifted high-symmetry k-point `kᴳ′` and return the
# associated connection `c′` between `kᴳ′` and `kᴴ`
function find_monodromy_related_connection(c, lgs)
    kᴳ, kᴴ = c.kᴳ.kv, c.kᴴ.kv

    kᴴlab = c.kᴴ.label
    kᴳlab = c.kᴳ.label
    kᴳ′lab = Symbol(kᴳlab, '′')
    lgᴴ = lgs[string(kᴴlab)]

    cntr = centering(num(lgᴴ), 3)
    kᴳₚ = primitivize(kᴳ, cntr)
    kᴴₚ = primitivize(kᴴ, cntr)
    # FIXME/TODO: the below approach to get the "translated" k-coordinate
    #             is still not right: must fix
    #             (Yes, it doesn't make sense to just add a _real-space_
    #              translation to the k-vector. Need to get k⋅t = -2πn with 
    #              smallest possible k (and in primitive setting))

    t = find_highest_order_nonsymmorphic_operation(lgᴴ, cntr)[2]
    tₚ  = primitivize(RVec(t), cntr).cnst

    free_parts = vec(transpose(tₚ) * kᴴₚ.free)
    idxs = findall(v -> norm(v) > Crystalline.DEFAULT_ATOL, free_parts)
    idx = if length(idxs) == 1
        last(idxs)
    elseif length(idxs) == 2
        # we cannot uniquely determine which k-point to extend to in this situation;
        # for now, we just pick _a_ solution, hoping that it will end up being the one that
        # will give us the same subduction table (TO BE VERIFIED)
        last(idxs)
    else
        # hoping that this scenario never occurs
        _error_info("unexpected situation when trying to solve k⋅t = 1; too many free parameters", kᴳlab, kᴴlab, kᴳ, kᴳₚ, kᴴ, kᴴₚ, "⋅", "⋅", t, tₚ, cntr)
    end
    kᴳ′ₚ = search_n2π_solution(free_parts, idx, kᴳₚ, kᴴₚ)
    if kᴳ′ₚ === nothing
        _error_info("failed to find a kᴳ′ that differs from from kᴳ by a primitive reciprocal lattice vector", kᴳlab, kᴴlab, kᴳ, kᴳₚ, kᴴ, kᴴₚ, kᴳ′, kᴳ′ₚ, t, tₚ, cntr)
    end
    Δkᴳₚ = (kᴳ′ₚ - kᴳₚ).cnst
    if !(dot(Δkᴳₚ, tₚ) ≈ round(dot(Δkᴳₚ, tₚ)))
        _error_info("Δkᴳ⋅t ∉ ℤ", kᴳlab, kᴴlab, kᴳ, kᴳₚ, kᴴ, kᴴₚ, kᴳ′, kᴳ′ₚ, t, tₚ, cntr)
    end

    kᴳ′ = conventionalize(kᴳ′ₚ, cntr)
    c′ = Connection(LabeledKVec(kᴳ′lab, kᴳ′), c.kᴴ)

    return c′
end

function search_n2π_solution(free_parts, idx′, kᴳₚ, kᴴₚ, n::Integer=1, nmax::Integer=6)
    # solve (kᴳ-kᴳ′)⋅t = -n2π for least positive integer `n`, stopping at `nmax`
    n>nmax && return nothing
    αβγ_idx′ = -n/free_parts[idx′]
    kᴳ′ₚ = KVec(kᴳₚ.cnst + kᴴₚ.free[:,idx′] * αβγ_idx′)
    # test that kᴳ′ differs from kᴳ by a primitive reciprocal lattice vector
    Δkᴳₚ = (kᴳ′ₚ - kᴳₚ).cnst
    if all(v->norm(round(v)-v)<Crystalline.DEFAULT_ATOL, Δkᴳₚ)
        return kᴳ′ₚ
    else
        # try to see if the solution could be fixed by using the remaining degrees of
        # freedom in kᴴₚ
        # NB: Below is not a very careful approach; it's a hail Mary & we should probably do
        #     do better to be safe. For now, we settled for this.
        idx′′ = findfirst(i->i≠idx′ && Crystalline.freeparams(kᴴₚ)[i], 1:3)
        if !isnothing(idx′′)
            v = kᴴₚ.free[:,idx′′]
            for c in (1.0, 0.5, -0.5, -1.0)
                kᴳ′ₚ = KVec(kᴳ′ₚ.cnst + c*v)
                Δkᴳₚ = (kᴳ′ₚ - kᴳₚ).cnst
                if all(v->norm(round(v)-v)<Crystalline.DEFAULT_ATOL, Δkᴳₚ)
                    return kᴳ′ₚ
                end
            end
        end
    end
    # continue searching for a higher `n` solution
    search_n2π_solution(free_parts, idx′, kᴳₚ, kᴴₚ, n+1, nmax)
end

function pseudo_inverse_smith(A)
    # compute the pseudo-inverse of a matrix `A` using the Smith normal form
    # (https://en.wikipedia.org/wiki/Pseudoinverse#Using_the_Smith_normal_form)
    F = Crystalline.smith(A)
    A⁺ = F.Tinv * pinv(diagm(F)) * F.Sinv
    return A⁺
end

function _error_info(msg, kᴳlab, kᴴlab, kᴳ, kᴳₚ, kᴴ, kᴴₚ, kᴳ′, kᴳ′ₚ, t, tₚ, cntr)
    printstyled(
        """
        $kᴳlab → $kᴴlab (cntr: $cntr)
          kᴳ:  $kᴳ \t\tkᴳₚ:  $kᴳₚ
          kᴳ′: $kᴳ′\t\tkᴳ′ₚ: $kᴳ′ₚ
          kᴴ:  $kᴴ \t\tkᴴₚ:  $kᴴₚ
          t:   $t  \ttₚ:   $tₚ\n"""; color=:yellow)
    error(msg)
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

## --------------------------------------------------------------------------------------- #

timereversal = false
connectionsd = Dict{Int, Vector{Connection{3}}}()
subductionsd = Dict{Int, Vector{SubductionTable{3}}}()
data = Dict{Int, Any}()
for sgnum in 3:230
    # Note: SGs 1 & 2 fail since there are no connections
    println("sgnum ", sgnum)
    tmp = crawl_bandpaths(sgnum, timereversal)
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
dir = joinpath((@__DIR__), "..", "BandGraphs/data/connections/3d/")
jldsave(joinpath(dir, "subductions$(timereversal ? "-tr" : "").jld2"); 
        subductionsd=subductionsd)
jldsave(joinpath(dir, "connections$(timereversal ? "-tr" : "").jld2"); 
         connectionsd=connectionsd)