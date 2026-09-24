# conversion between delimited text files and array representations
dlm2array(io::IO) = DelimitedFiles.readdlm(io, '|', String, '\n')
dlm2array(str::String) = dlm2array(IOBuffer(str))

# utilities for creation of BandRep and BandRepSet
function dlm2struct(str::Union{String,IO}, sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, 
                    timereversal::Bool=true)
    M = dlm2array(str);
    array2struct(M, sgnum, allpaths, spinful, timereversal)
end

function array2struct(M::Matrix{String}, sgnum::Integer, allpaths::Bool=false, spinful::Bool=false, 
                      timereversal::Bool=true)

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
    # A BandRepSet can either reference single-valued or double-valued irreps, not both; 
    # thus, we "throw out" one of the two here, depending on `spinful`.
    if spinful  # double-valued irreps only (spinful systems)
        delidxs = findall(map(!isspinful, brtags))
    else        # single-valued irreps only (spinless systems)
        delidxs = findall(map(isspinful, brtags))
    end
    for vars in (brtags, wyckpos, sitesym, label, dim)
        deleteat!(vars, delidxs) 
    end
    irlabs, irvecs = get_irrepvecs(brtags)              

    BRs = BandRep.(wyckpos, sitesym, label, dim, map(isspinful, brtags), irvecs,
                   Ref(irlabs))
    
    return BandRepSet(sgnum, BRs, kvs, klabs, irlabs, spinful, timereversal)
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


"""
    bandreps(sgnum::Integer, D::Integer=3; 
             allpaths::Bool=false, spinful::Bool=false, timereversal::Bool=true) --> BandRepSet

Return the elementary band representations (EBRs) as a `BandRepSet` for space group `sgnum`
and dimension `D`.

## Keyword arguments

- `allpaths`: include a minimal sufficient set (`false`, default) or all (`true`) 
              **k**-vectors. 
- `spinful`: single- (`false`, default) or double-valued (`true`) irreps, as appropriate for
             spinless and spinful particles, respectively. Only available for `D=3`.
- `timereversal`: assume presence (`true`, default) or absence (`false`) of time-reversal
                  symmetry.

## References
3D EBRs are obtained from the Bilbao Crystallographic Server's 
[BANDREP program](https://cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl);
please reference the original research papers noted there if used in published work.
"""
function bandreps(sgnum::Integer, D::Integer=3;
                  allpaths::Bool=false, spinful::Bool=false,
                  timereversal::Bool=true)
    D ∉ (1,2,3) && _throw_invalid_dim(D)
    paths_str = allpaths ? "allpaths" : "maxpaths"
    brtype_str = timereversal ? "elementaryTR" : "elementary"
    filename = joinpath(DATA_DIR, 
                        "bandreps/$(D)d/$(brtype_str)/$(paths_str)/$(string(sgnum)).csv")
    open(filename) do io
        brs = dlm2struct(io, sgnum, allpaths, spinful, timereversal)
    end 
end

"""
    basisdim(brs::BandRepSet) --> Int

Return the dimension of the (linearly independent parts) of a band representation set.
This is ``d^{\\text{bs}} = d^{\\text{ai}}`` in the notation of [Po, Watanabe, & Vishwanath,
Nature Commun. **8**, 50 (2017)](https://doi.org/10.1038/s41467-017-00133-2), or 
equivalently, the rank of `stack(brs)` over the ring of integers.
This is the number of linearly independent basis vectors that span the expansions of
a band structure viewed as symmetry data.
""" 
function basisdim(brs::BandRepSet)
    Λ = smith(stack(brs)).SNF
    nnz = count(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix
    return nnz
end

# misc minor utility functions
isspinful(br::AbstractVector{T} where T<:AbstractString) = any(x->occursin(r"\\bar|ˢ", x), br)

function split_paren(str::AbstractString)
    openpar = something(findfirst(==('('), str)) # index of the opening parenthesis
    before_paren = SubString(str, firstindex(str), prevind(str, openpar))
    inside_paren = SubString(str, nextind(str, openpar), prevind(str, lastindex(str)))
    return before_paren, inside_paren
end
