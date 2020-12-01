using JSON2, Crystalline, HTTP, Gumbo
using StaticArrays

struct RVec{D} <: Crystalline.AbstractVec
    cnst::SVector{D, Float64}
    free::Matrix{Float64}      # should do the tuple of tuples trick...
end
struct WyckPos{D}
    letter::Char
    mult::Int
    q::RVec{D} # just keep a single representative; assume assoc. w/ some spacegroup setting
    sitesym::Crystalline.GenericGroup{D}
    #cosets::Vector{SymOperation{D}}
end

# ---------------------------------------------------------------------------------------- %
# essentially a copy of the KVec parsing mechanism and the KVec printing for now...
function _strip_split(str::AbstractString)
    str = filter(!isspace, strip(str, ['(',')','[',']'])) # tidy up string (remove parens & spaces)
    return split(str,',')
end
RVec(str::AbstractString) = (xyz = _strip_split(str); _RVec(xyz, Val(length(xyz))))
RVec{D}(str::AbstractString) where D = (xyz = _strip_split(str); _RVec(xyz, Val{D}()))
function _RVec(xyz::Vector{<:SubString}, ::Val{D}) where D
    length(xyz) == D || throw(DimensionMismatch("Dimension D doesn't match input string"))
    cnst = zero(MVector{D, Float64})
    free = zeros(Float64, D, D)
    for (i, coord) in enumerate(xyz)
        # --- "free" coordinates, free[i,:] ---
        for (j, matchgroup) in enumerate((('α','u','x'),('β','v','y'),('γ','w','z')))
            pos₂ = findfirst(∈(matchgroup), coord)
            if !isnothing(pos₂)
                prefix = Crystalline.searchpriornumerals(coord, pos₂)
                free[i,j] = parse(Float64, prefix)
            end
        end
        
        # --- "fixed"/constant coordinate, cnst[i] ---
        m = match(r"(?:\+|\-)?(?:(?:[0-9]|/|\.)+)(?!(?:[0-9]|\.)*[αuxβvyγwz])", coord)
        # regex matches any digit sequence, possibly including slashes, that is _not_
        # followed by one of the free-part identifiers. If there's a '+' or '-' before
        # the first digit, it is stored in the first capture slot. The digit sequence
        # is stored in the second capture slot. The third capture slot is redundant.
        # We do not allow arithmetic aside from division here, obviously: any extra numbers 
        # terms are ignored.
        if m===nothing # no constant terms
            if last(coord) ∈ ('α','u','β','v','γ','w','x','y','z') # free-part only case
                continue # cnst[i] is zero already
            else
                throw(ErrorException("Unexpected parsing error in constant term"))
            end
        else # exploit that we require constant parts to come before free parts
            cnst[i] = Crystalline.parsefraction(m.captures[2])
            m.captures == '-' && (cnst[i] *= -1)
        end
    end
    cnst = convert(SVector{D}{Float64}, cnst)
    return RVec{D}(cnst, free)
end

# ---------------------------------------------------------------------------------------- %
const BILBAO_URL = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/"
const WYCK_URL_BASE_3D = BILBAO_URL*"nph-normsets?from=wycksets&gnum="
wyck_url(sgnum) = WYCK_URL_BASE_3D*string(sgnum)

# ----- NOW-REDUNANT FUNCTIONS FOR CRAWLING 3D SPACE GROUPS FROM BILBAO -----
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
        letter, mult_str, sitesym_str, q_str, _ = getfield.(first.(getfield.(el, Ref(:children))), Ref(:text))

        println(mult_str, letter, ": ", filter(!isspace, q_str), " (SS = ", sitesym_str, ")")

        q = RVec{D}(q_str)
        mult = parse(Int, mult_str) 

        # TODO: compute site symmetry group!
        sitesym = Crystalline.GenericGroup{D}(SymOperation{D}[])
        # TODO: compute coset representatives for other positions in orbit

        wycks[i] = WyckPos{D}(only(letter), mult, q, sitesym)
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