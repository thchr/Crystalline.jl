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
             allpaths::Bool=false, spinful::Bool=false, timereversal::Bool=true)

Returns the elementary band representations (EBRs) as a `BandRepSet` for space group `sgnum`
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
[BANDREP program](http://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl);
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
        BRS = dlm2struct(io, sgnum, allpaths, spinful, timereversal)
    end 
end

"""
$(TYPEDSIGNATURES)

Return the nontrivial (i.e., ≠ {0,1}) elementary factors of an EBR basis, provided as a
`BandRepSet` or `Smith` decomposition.
"""
function nontrivial_factors(F::Smith)
    Λ = F.SNF
    nontriv_idx = findall(is_not_one_or_zero, Λ)
    return Λ[nontriv_idx]
end
function nontrivial_factors(brs::BandRepSet)
    F = smith(stack(brs), inverse=false)
    return nontrivial_factors(F)
end
is_not_one_or_zero(x) = !(isone(x) || iszero(x))

"""
    classification(BRS_or_F::Union{BandRepSet, Smith}) --> String

Return the symmetry indicator group ``X^{\\text{BS}}`` of an EBR basis `F_or_BRS`, provided
as a `BandRepSet` or `Smith` decomposition.

Technically, the calculation answers the question "what direct product of 
``\\mathbb{Z}_n`` groups is the the quotient group
``X^{\\text{BS}} = \\{\\text{BS}\\}/\\{\\text{AI}\\}`` isomorphic to?" (see
[Po, Watanabe, & Vishwanath, Nature Commun. **8**, 50 (2017)](https://doi.org/10.1038/s41467-017-00133-2)
for more information).
"""
function classification(nontriv_Λ::AbstractVector{<:Integer})
    if isempty(nontriv_Λ)
        return "Z₁"
    else
        return "Z"*subscriptify(join(nontriv_Λ, "×Z"))
    end
end
function classification(BRS_or_F::Union{BandRepSet, Smith})
    return classification(nontrivial_factors(BRS_or_F))
end

"""
    basisdim(BRS::BandRepSet) --> Int

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

"""
    wyckbasis(brs::BandRepSet) --> Vector{Vector{Int}}

Computes the (band representation) basis for bands generated by localized orbitals placed at
the Wyckoff positions. Any band representation that can be expanded on this basis with
positive integer coefficients correspond to a trivial insulators (i.e. deformable to atomic
limit).
Conversely, bands that cannot are topological, either fragily (some negative coefficients)
or strongly (fractional coefficients).
"""
function wyckbasis(brs::BandRepSet) 
    # Compute Smith normal form: for an n×m matrix B with integer elements,
    # find matrices S, diagm(Λ), and T (of size n×n, n×m, and m×m, respectively)
    # with integer elements such that B = S*diagm(Λ)*T. Λ is a vector
    # [λ₁, λ₂, ..., λᵣ, 0, 0, ..., 0] with λⱼ₊₁ divisible by λⱼ and r ≤ min(n,m).
    # The matrices T and S have integer-valued pseudo-inverses.
    F = smith(stack(brs)) # Smith normal factorization with λⱼ ≥ 0
    S, S⁻¹, Λ = F.S, F.Sinv, F.SNF
    #T, T⁻¹ = F.T, F.Tinv,

    nnz = count(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix
    nzidxs = OneTo(nnz)

    # If we apply S⁻¹ to a given set of (integer) symmetry data 𝐧, the result 
    # should be the (integer) factors qᵢCᵢ (Cᵢ=Λᵢ here) discussed in Tang, Po,
    # [...], Nature Physics 15, 470 (2019). Conversely, the columns of S gives an integer-
    # coefficient basis for all gapped band structures, while the columns of S*diagm(Λ)
    # generates all atomic insulator band structures (assuming integer coefficients).
    # See also your notes in scripts/derive_sg2_bandrep.jl
    return S[:, nzidxs], diagm(F)[:,nzidxs], S⁻¹[:, nzidxs]
end


# misc minor utility functions
isspinful(br::AbstractVector{T} where T<:AbstractString) = any(x->occursin(r"\\bar|ˢ", x), br)

function split_paren(str::AbstractString)
    openpar = something(findfirst(==('('), str)) # index of the opening parenthesis
    before_paren = SubString(str, firstindex(str), prevind(str, openpar))
    inside_paren = SubString(str, nextind(str, openpar), prevind(str, lastindex(str)))
    return before_paren, inside_paren
end

# TODO: Remove this (unexported method)
"""
    matching_littlegroups(BRS::BandRepSet, ::Val{D}=Val(3))

Finds the matching little groups for each *k*-point referenced in `BRS`. This is mainly a 
a convenience accessor, since e.g. [`littlegroup(::SpaceGroup, ::KVec)`](@ref) could also
return the required little groups. The benefit here is that the resulting operator sorting
of the returned little group is identical to the operator sortings assumed in
[`lgirreps`](@ref) and [`littlegroups`](@ref).

Returns a `Vector{LittleGroup{D}}` (unlike [`littlegroups`](@ref), which returns a
`Dict{String, LittleGroup{D}}`).

## Note 1

Unlike the operations returned by [`spacegroup`](@ref), the returned little groups do not
include copies of operators that would be identical when transformed to a primitive basis.
The operators are, however, still given in a conventional basis.
"""
function matching_littlegroups(BRS::BandRepSet, ::Val{D}=Val(3)) where D
    lgs = littlegroups(num(BRS), Val(D)) # TODO: generalize to D≠3

    klabs_in_brs = klabels(BRS) # find all k-point labels in BandRepSet
    if !issubset(klabs_in_brs, keys(lgs))
        throw(DomainError(klabs_in_brs, "Could not locate all LittleGroups from BandRep"))
    end

    return getindex.(Ref(lgs), klabs_in_brs)
end


function matching_lgirreps(BRS::BandRepSet)
    lgirsd = lgirreps(num(BRS), Val(3))
    # create "physical/real" irreps if `BRS` assumes time-reversal symmetry
    if BRS.timereversal 
        for (klab, lgirs) in lgirsd
            lgirsd[klab] = realify(lgirs)
        end
    end
    # all lgirreps from ISOTROPY as a flat vector; note that sorting is arbitrary
    lgirs = collect(Iterators.flatten(values(lgirsd))) # TODO: generalize to D≠3

    # find all the irreps in lgirs that feature in the BandRepSet, and 
    # sort them according to BandRepSet's sorting
    lgirlabs = label.(lgirs)
    brlabs   = normalizesubsup.(irreplabels(BRS))

    find_and_sort_idxs = Vector{Int}(undef, length(brlabs))
    @inbounds for (idx, brlab) in enumerate(brlabs)
        matchidx = findfirst(==(brlab), lgirlabs)
        if matchidx !== nothing
            find_and_sort_idxs[idx] = matchidx
        else
            throw(DomainError(brlab, "could not be found in ISOTROPY dataset"))
        end
    end

    # find all the irreps _not_ in the BandRepSet; keep existing sorting
    not_in_BRS_idxs = filter(idx -> idx∉find_and_sort_idxs, eachindex(lgirlabs))

    # return: 1st element = lgirs ∈ BandRepSet (matching)
    #         2nd element = lgirs ∉ BandRepSet (not matching)
    return (@view lgirs[find_and_sort_idxs]), (@view lgirs[not_in_BRS_idxs])
end