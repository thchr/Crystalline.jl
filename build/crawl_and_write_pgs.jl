using Crystalline, HTTP, Gumbo, JLD2

const BILBAO_PG_URL_BASE = "https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_out.pl?tipogrupo=spg&pointspace=point&"

function bilbao_pgs_specifier(parent_sgnum::Integer, pgnum::Integer, pgiuc::String)
    # in principle, it is sufficient to only supply the appropraite parent space group
    # `parent_sgnum`. Note that it is _not_ sufficient to only supply `pgnum` due to 
    # different setting variation possibilities. The inclusion of `pgiuc` appear mostly
    # immaterial
    return "num=$(parent_sgnum)&super=$(pgnum)&symbol=$(pgiuc)"
    # the url to a given point group is BILBAO_PG_URL_BASE*bilbao_pgs_specifier(..)
end

function findfirst_matching_parent_sgnum(pgiuc::String)
    # The 1D and 2D cases could be inferred from the 3D case, and Bilbao only lists the 3D
    # case, so this is restricted to the 3D case (note though, that the settings are not 
    # shared between 1D/2D and 3D, unfortunately (e.g. "m"))
    for sgnum in 1:MAX_SGNUM[3]
        sg = spacegroup(sgnum, Val(3))
        pg = Crystalline.find_parent_pointgroup(sg)
        label(pg) == pgiuc && return sgnum, num(pg) # return parent_sgnum, pgnum
    end
    throw(DomainError(pgiuc, "requested label cannot be found"))
end

function bilbao_pgs_url(pgiuc::String)

    parent_sgnum, pgnum = findfirst_matching_parent_sgnum(pgiuc)
    url = BILBAO_PG_URL_BASE*bilbao_pgs_specifier(parent_sgnum, pgnum, pgiuc)
    
    return url
end

if !isdefined(@__MODULE__,:PGS_HTML_REPLACEMENTS)
    const PGS_HTML_REPLACEMENTS = (
        r"<font style=\\\"text-decoration: overline;\\\">(.*?)</font>"=>s"-\1",
                                                                # ↑ overbar ⇒ minus sign
        r"<sub>(.*)</sub>"=>x->Crystalline.subscriptify(x[6:end-6]),  # subscripts
        r"<sup>(.*)</sup>"=>x->Crystalline.supscriptify(x[6:end-6]),  # superscripts
    )
end

function _parse_pgs_elements(s::AbstractString)
    if first(s) === 'e'
        # the form of ϕ_str will always be "(-)iNπ/M" with integers N and M: to just get the 
        # paesing job done, I resort to the following hardcoded approach:
        ϕ_str_nom, ϕ_str_den = getfield(match(r"e<sup>(.*?)/(.*)</sup>", s), :captures)
        ϕ_str_nom′ = replace(replace(ϕ_str_nom, 'i'=>""), 'π'=>"") # strip i and π
        ϕ_nom′ = if ϕ_str_nom′ == "-"  # special handling for N = 1
           -1.0
        elseif ϕ_str_nom′ == "" 
            1.0
        else
            parse(Float64, ϕ_str_nom′)
        end
        ϕ = π*ϕ_nom′/parse(Float64, ϕ_str_den)
        return cis(ϕ)
    else
        # first, handle corner case for parse(ComplexF64, s): cannot parse "±i", but "±1i" 
        # is OK (substitution below hacks around fact that s"\11i" refers to capture 11,
        # rather than capture 1 + string 1i)
        if occursin('i', s)
            s = replace(s, r"(\+|\-)i"=>m->m[1]*"1i") 
            if s == "i"; s = "1i"; end # not covered by regex above
        end
        return parse(ComplexF64, s)
    end
end

""" 
    crawl_pgirs(pgiuc::String, D::Integer=3, consistency_checks::Bool=true)

Crawl the point group symmetry operations as well as the associated point group irreps for
the the point group characterized by the IUC label `pgiuc` and dimension `D` from Bilbao's
website. Only works for `D = 3`, since Bilbao doesn't tabulate 1D/2D point group irreps 
explicitly. Returns a `Vector{PGIrrep{D}}`.
"""
function crawl_pgirs(pgiuc::String, D::Integer=3; consistency_checks::Bool=true)
    D ≠ 3 && throw(DomainError(D, "dimensions D≠3 not crawlable"))

    # --- crawl html data from Bilbao ---
    contents = HTTP.request("GET", bilbao_pgs_url(pgiuc))
    body = parsehtml(String(contents.body))
    # pick out the table that contains the irreps, including header
    table_html = children(last(children(body.root)))[end-7] # -7 part is hardcoded
    # strip out the header row and split at rows (i.e. across pgops)
    # header is [N | pgop matrix | pgop seitz | irrep 1 | irrep 2 | ...])
    header_html = table_html[1][1]
    rows_html = children(table_html[1])[2:end]
    # split each row into columns
    rowcols_html = children.(rows_html)
    Nirs = length(first(rowcols_html))-3

    # --- extract point group operators ---
    # extract pgops' matrix form from 2nd column
    pgops_html = getindex.(rowcols_html, 2)
    pgops_str = replace.(replace.(string.(pgops_html), 
                                        Ref(r".*?<td align=\"right\">(.*?)</td>"=>s"\1")), # extra matrix elements & remove "front"
                               Ref(r"</tr>.*"=>"")) # remove "tail"
    pgops_strs = split.(pgops_str, " ", keepempty=false) # matrix elements, row by row
    pgops_matrices = [collect(reshape(parse.(Int, strs), (3,3))') for strs in pgops_strs]
    pgops     = SymOperation{3}.(hcat.(pgops_matrices, Ref(zeros(Int, 3))))
    # build point group PointGroup{3} from extracted pgops
    pgnum = Crystalline.pointgroup_iuc2num(pgiuc, 3)
    pg = PointGroup{3}(pgnum, pgiuc, pgops) # note that sorting matches pointgroup(...)
    if consistency_checks
        # test that pgops match results from pointgroup(pgiuc, 3), incl. matching sorting
        pg′ = pointgroup(pgiuc, Val(3))
        @assert operations(pg) == operations(pg′)

        #= # ~~~ dead code ~~~
        # extract pgops' seitz notation from 3rd columns 
        seitz_html = children.(first.(children.(getindex.(rowcols_html, 3))))
        seitz_str = [join(string.(el)) for el in seitz_html]
        for replacepair in PGS_HTML_REPLACEMENTS
            seitz_str .= replace.(seitz_str, Ref(replacepair))
        end
        # | test whether we extract pgops and seitz notation mutually consistently; this
        # | isn't really possible to test for neatly in general, unfortunately, due to an
        # | indeterminacy of conventions in cases like 2₁₋₁₀ vs 2₋₁₁₀; at the moment.
        # | We keep the code in case it needs to be used for debugging at some point.
        @assert all(replace.(seitz.(pgops), Ref(r"{(.*)\|0}"=>s"\1")) .== seitz_str)
        =#
    end
    
    # --- extract irreps ---
    irreps_html = [string.(getindex.(rowcols_html, i)) for i in 4:3+Nirs] # indexed first over distinct irreps, then over operators
    irreps_mats = Vector{Vector{Matrix{ComplexF64}}}(undef, Nirs)
    for (i,irs_html) in enumerate(irreps_html)
        # strip non-essential "pre/post" html; retain matrix/scalar as nested row/columns
        matrix_tabular = 
            first.(getfield.(match.(r".*?<table><tbody>(.*?)</tbody></table>.*", irs_html),
                                    :captures))
        matrix_els_regexs = 
            eachmatch.(r"<td align=\\\"(?:left|right|center)\\\">\s*(.*?)</td>", 
                       matrix_tabular)
        matrix_els_strs = [first.(getfield.(els, :captures)) for els in matrix_els_regexs]
        # matrix_elements′ is an array of "flattened" matrices with string elements; some of 
        # these elements can be of the form e<sup>(.*?)</sup> for an exponential, the rest
        # are integers. Now we parse these strings, taking care to allow exponential forms
        matrix_els = [_parse_pgs_elements.(els) for els in matrix_els_strs] # ::Vector{C64}
        Dⁱʳ = isqrt(length(first(matrix_els)))
        # a vector of matrices, one for each operation in pgops, for the ith irrep
        irreps_mats[i] = [collect(transpose(reshape(els, Dⁱʳ, Dⁱʳ))) for els in matrix_els]
    end

    # --- extract irrep labels and irrep reality ---
    #                                                          ↓ ::Vector{Vector{HTMLNode}}
    labs_and_realities_html = children.(getindex.(children(header_html)[4:end],1))
    irrreps_labs_html = join.(getindex.(labs_and_realities_html,         # ::Vector{String}
                                Base.OneTo.(lastindex.(labs_and_realities_html) .- 1)))
    irreps_labs = replace.(replace.(irrreps_labs_html, Ref(PGS_HTML_REPLACEMENTS[2])), 
                    Ref(PGS_HTML_REPLACEMENTS[3]))
    irreps_labs .= replace.(irreps_labs, Ref("GM"=>"Γ"))
    irreps_realities = parse.(Int64, strip.(string.(getindex.(labs_and_realities_html, 
                                                    lastindex.(labs_and_realities_html))),
                                            Ref(['(', ')'])))
    
    # --- return a vector of extracted point group irreps ::Vector{PGIrrep{3}} ---
    return PGIrrep{3}.(irreps_labs, Ref(pg), irreps_mats, irreps_realities)
end




function __crawl_and_write_3d_pgirreps()
    savepath = (@__DIR__)*"/../data/pgirreps/3d/"


    # we only save the point group irreps.
    filename_irreps = savepath*"/pgirreps_data"

    JLD2.jldopen(filename_irreps*".jld2", "w") do irreps_file
        for pgiuc in PGS_IUCs[3]
            # ==== crawl and prepare irreps data ====
            pgirs = crawl_pgirs(pgiuc, 3; consistency_checks=true)
            matrices = [pgir.matrices for pgir in pgirs]
            types = [type(pgir) for pgir in pgirs]
            cdmls = [label(pgir) for pgir in pgirs]

            # ==== save irreps ====
            # we do not save the point group operations anew; they are already stored in 
            # "data/pgops/..."; note that we explicitly checked the sorting and equivalence 
            # of operations when pgirs was crawled above (cf. flag `consistency_checks=true`)
            unmangled_pgiuc = Crystalline.unmangle_pgiuclab(pgiuc) # replace '/'s by '_slash_'s
            irreps_file[unmangled_pgiuc*"/matrices"] = matrices
            irreps_file[unmangled_pgiuc*"/types"] = types
            irreps_file[unmangled_pgiuc*"/cdmls"] = cdmls
        end
    end
    return filename_irreps
end





# ======================================================
#= 
pgirs = Vector{Vector{PGIrrep{3}}}(undef, length(PGS_IUCs[3]))
for (idx, pgiuc) in enumerate(PGS_IUCs[3])
    println(pgiuc)
    pgirs[idx] = crawl_pgirs(pgiuc)
end
=#
         
# Save to local format
__crawl_and_write_3d_pgirreps()

