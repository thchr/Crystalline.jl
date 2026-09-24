# The elementary band representation (EBR) tables of the Bilbao Crystallographic Server's
# [BANDREP program](https://cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl), and the parsing
# needed to read them.
#
# Crystalline computes its own EBRs with `bandreps`; these tables exist only to validate
# that computation against an independent reference, and so live with the tests rather than
# in the package. The underlying CSV files are a lazy artifact (see `Artifacts.toml`).
module BilbaoBandReps

using Crystalline
using DelimitedFiles: readdlm
using Pkg.Artifacts: ensure_artifact_installed
using Base: OneTo, @propagate_inbounds

export BilbaoBandRep, BilbaoBandRepSet, bilbao_bandreps

# The CSV tables are a lazy artifact rather than part of the repository. `@artifact_str`
# cannot resolve it from here: `find_artifacts_toml` walks up from the calling file and stops
# at the first `Project.toml`, which is `test/Project.toml`, so name the package's own.
const BANDREPS_DATA_DIR = Ref{String}()
function bandreps_datadir()
    if !isassigned(BANDREPS_DATA_DIR)
        BANDREPS_DATA_DIR[] = ensure_artifact_installed(
                "bandreps", joinpath(pkgdir(Crystalline), "Artifacts.toml"))
    end
    return BANDREPS_DATA_DIR[]
end

# ---------------------------------------------------------------------------------------- #
# Types

"""
    BilbaoBandRep <: AbstractVector{Int}

A single elementary band representation, as tabulated by BANDREP.
"""
struct BilbaoBandRep <: AbstractVector{Int}
    wyckpos :: String       # Wyckoff position that induces the BR
    sitesym :: String       # site symmetry point group of Wyckoff position (IUC notation)
    label   :: String       # symbol ρ↑G, with ρ the irrep of the site symmetry group
    dim     :: Int          # dimension (i.e., number of bands) of the band representation
    spinful :: Bool         # whether the BR involves spinful irreps ("\\bar"'ed irreps)
    irvec   :: Vector{Int}  # references `irlabs`; nonzero entries are in the band rep
    irlabs  :: Vector{String} # labels, as in the parent `BilbaoBandRepSet`
end
Base.position(br::BilbaoBandRep) = br.wyckpos
Crystalline.label(br::BilbaoBandRep) = br.label
Crystalline.irreplabels(br::BilbaoBandRep) = br.irlabs
Crystalline.dim(br::BilbaoBandRep) = br.dim

Base.size(br::BilbaoBandRep) = (length(br.irvec) + 1,) # + 1 to include the filling
@propagate_inbounds function Base.getindex(br::BilbaoBandRep, i::Int)
    return i == length(br.irvec)+1 ? dim(br) : br.irvec[i]
end
Base.IndexStyle(::Type{<:BilbaoBandRep}) = IndexLinear()

"""
    BilbaoBandRepSet <: AbstractVector{BilbaoBandRep}

The elementary band representations of a single space group, as tabulated by BANDREP.
"""
struct BilbaoBandRepSet <: AbstractVector{BilbaoBandRep}
    sgnum        :: Int                  # space group number
    bandreps     :: Vector{BilbaoBandRep}
    kvs          :: Vector{<:KVec}       # 𝐤-points
    klabs        :: Vector{String}       # associated 𝐤-labels (CDML notation)
    irlabs       :: Vector{String}       # (sorted) CDML irrep labels at _all_ 𝐤-points
    spinful      :: Bool                 # whether spinful irreps are included
    timereversal :: Bool                 # whether time-reversal symmetry is assumed
end
Crystalline.num(brs::BilbaoBandRepSet) = brs.sgnum
Crystalline.klabels(brs::BilbaoBandRepSet) = brs.klabs
Crystalline.irreplabels(brs::BilbaoBandRepSet) = brs.irlabs
Crystalline.isspinful(brs::BilbaoBandRepSet) = brs.spinful

Base.size(brs::BilbaoBandRepSet) = (length(brs.bandreps),)
@propagate_inbounds Base.getindex(brs::BilbaoBandRepSet, i::Int) = brs.bandreps[i]
Base.IndexStyle(::Type{<:BilbaoBandRepSet}) = IndexLinear()

Base.stack(brs::BilbaoBandRepSet) = reduce(hcat, brs)

function Base.show(io::IO, ::MIME"text/plain", brs::BilbaoBandRepSet)
    print(io, "BilbaoBandRepSet (⋕", num(brs), "): ",
              length(brs), " BandReps, ",
              "sampling ", length(irreplabels(brs)), " LGIrreps ",
              "(", isspinful(brs) ? "spinful" : "spinless", " ",
              brs.timereversal ? "w/" : "w/o", " TR)")
end

# ---------------------------------------------------------------------------------------- #
# Conversion from Crystalline's own band representations, for comparison

function Base.convert(::Type{BilbaoBandRep}, br::Crystalline.BandRep)
    return BilbaoBandRep(label(position(br.siteir)),
                         br.siteir.pglabel,
                         label(br.siteir)*"↑G",
                         occupation(br),
                         isspinful(br),
                         collect(br)[1:end-1],
                         irreplabels(br))
end

function Base.convert(::Type{BilbaoBandRepSet},
                      brs::Collection{<:Crystalline.BandRep})
    return BilbaoBandRepSet(num(brs),
                            convert.(Ref(BilbaoBandRep), brs),
                            [position(lgirs) for lgirs in irreps(brs)],
                            klabels(brs),
                            irreplabels(brs),
                            isspinful(first(brs)),
                            first(brs).timereversal)
end

# ---------------------------------------------------------------------------------------- #
# Parsing of the BANDREP tables

# conversion between delimited text files and array representations
dlm2array(io::IO) = readdlm(io, '|', String, '\n')
dlm2array(str::String) = dlm2array(IOBuffer(str))

function dlm2struct(str::Union{String,IO}, sgnum::Integer, allpaths::Bool=false,
                    spinful::Bool=false, timereversal::Bool=true)
    M = dlm2array(str)
    array2struct(M, sgnum, allpaths, spinful, timereversal)
end

function array2struct(M::Matrix{String}, sgnum::Integer, allpaths::Bool=false,
                      spinful::Bool=false, timereversal::Bool=true)

    klist =  permutedims(mapreduce(x->String.(split(x,":")), hcat, M[4:end,1])) # 1ˢᵗ col is labels, 2ⁿᵈ col is coordinates as strings
    klabs, kvs = (@view klist[:,1]), KVec.(@view klist[:,2])

    temp = split_paren.(@view M[1,2:end])
    wyckpos, sitesym = getindex.(temp, 1), getindex.(temp, 2) # wyckoff position and site symmetry point group of bandrep

    temp .= split_paren.(@view M[2,2:end]) # same size, so reuse array
    label, dim = getindex.(temp, 1), parse.(Int, getindex.(temp, 2)) # label of bandrep

    # whether M contains info on decomposability; we don't use this anymore, but need to
    # know to parse the rest of the contents correctly (we used to always include this info
    # but might not in the future; so protect against this)
    has_decomposable_info = M[3,1] == "Decomposable"
    # decomposable = parse.(Bool, vec(@view M[3,2:end])) # whether BR can be BR-decomposed

    # set of irreps that jointly make up the bandrep
    brtags = collect(eachcol(@view M[3+has_decomposable_info:end, 2:end]))
    for br in brtags
        br .= replace.(br, Ref(r"\([1-9]\)"=>""))  # get rid of irrep dimension info
    end
    # A band rep set can either reference single-valued or double-valued irreps, not both;
    # thus, we "throw out" one of the two here, depending on `spinful`.
    if spinful  # double-valued irreps only (spinful systems)
        delidxs = findall(map(!has_spinful_tag, brtags))
    else        # single-valued irreps only (spinless systems)
        delidxs = findall(map(has_spinful_tag, brtags))
    end
    for vars in (brtags, wyckpos, sitesym, label, dim)
        deleteat!(vars, delidxs)
    end
    irlabs, irvecs = get_irrepvecs(brtags)

    brs = BilbaoBandRep.(wyckpos, sitesym, label, dim, map(has_spinful_tag, brtags), irvecs,
                         Ref(irlabs))

    return BilbaoBandRepSet(sgnum, brs, kvs, klabs, irlabs, spinful, timereversal)
end

function get_irrepvecs(brtags)
    Nklabs = length(first(brtags)) # there's equally many (composite) irrep tags in each band representation
    irlabs = Vector{String}()
    for kidx in OneTo(Nklabs)
        irlabs_at_kidx = Vector{String}()
        for tag in getindex.(brtags, kidx) # tag could be a combination like Γ1⊕2Γ₂ (or something simpler, like Γ₁)
            for irrep in split(tag, '⊕')
                irrep′ = filter(!isdigit, irrep) # filter off any multiplicities
                if irrep′ ∉ irlabs_at_kidx
                    push!(irlabs_at_kidx, irrep′)
                end
            end
        end
        sort!(irlabs_at_kidx)
        append!(irlabs, irlabs_at_kidx)
    end

    irvecs = [zeros(Int, length(irlabs)) for _ in OneTo(length(brtags))]
    for (bridx, tags) in enumerate(brtags)
        for (kidx,tag) in enumerate(tags)
            for irrep in split(tag, '⊕') # note this irrep tag may contain numerical prefactors!
                buf = IOBuffer(irrep)
                prefac_str = readuntil(buf, !isdigit)
                seek(buf, ncodeunits(prefac_str)) # go back to first non-digit position in buffer
                if isempty(prefac_str)
                    prefac = Int(1)
                else
                    prefac = parse(Int, prefac_str)
                end
                ir′ = read(buf, String) # the rest of the irrep buffer is the actual cdml label
                close(buf)
                iridx = findfirst(==(ir′), irlabs) # find position in irlabs vector
                irvecs[bridx][iridx] = prefac
            end
        end
    end
    return irlabs, irvecs
end

has_spinful_tag(tags::AbstractVector{<:AbstractString}) = any(x->occursin(r"\\bar|ˢ", x), tags)

function split_paren(str::AbstractString)
    openpar = something(findfirst(==('('), str)) # index of the opening parenthesis
    before_paren = SubString(str, firstindex(str), prevind(str, openpar))
    inside_paren = SubString(str, nextind(str, openpar), prevind(str, lastindex(str)))
    return before_paren, inside_paren
end

# ---------------------------------------------------------------------------------------- #
# Loading

"""
    bilbao_bandreps(sgnum::Integer, D::Integer=3;
                    allpaths::Bool=false, spinful::Bool=false, timereversal::Bool=true)
                                                                    --> BilbaoBandRepSet

Return the elementary band representations (EBRs) tabulated by the Bilbao Crystallographic
Server's [BANDREP program](https://cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl) for space
group `sgnum` and dimension `D`.

## Keyword arguments

- `allpaths`: include a minimal sufficient set (`false`, default) or all (`true`)
              **k**-vectors.
- `spinful`: single- (`false`, default) or double-valued (`true`) irreps, as appropriate for
             spinless and spinful particles, respectively. Only available for `D=3`.
- `timereversal`: assume presence (`true`, default) or absence (`false`) of time-reversal
                  symmetry.
"""
function bilbao_bandreps(sgnum::Integer, D::Integer=3;
                         allpaths::Bool=false, spinful::Bool=false,
                         timereversal::Bool=true)
    D ∈ (1,2,3) || throw(DomainError(D, "dimension must be 1, 2, or 3"))
    paths_str = allpaths ? "allpaths" : "maxpaths"
    brtype_str = timereversal ? "elementaryTR" : "elementary"
    filename = joinpath(bandreps_datadir(),
                        "$(D)d/$(brtype_str)/$(paths_str)/$(string(sgnum)).csv")
    open(filename) do io
        dlm2struct(io, sgnum, allpaths, spinful, timereversal)
    end
end

end # module BilbaoBandReps
