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
    label, dim = getindex.(temp, 1), parse.(Int64, getindex.(temp, 2)) # label of bandrep

    decomposable = parse.(Bool, vec(@view M[3,2:end])) # whether bandrep can be further decomposed

    brtags = collect(eachcol(@view M[4:end, 2:end])) # set of irreps that jointly make up the bandrep
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
    for vars in (brtags, wyckpos, sitesym, label, dim, decomposable)
        deleteat!(vars, delidxs) 
    end
    irlabs, irvecs = get_irrepvecs(brtags)              

    BRs = BandRep.(wyckpos, sitesym, label, dim, decomposable, map(isspinful, brtags), 
                   irvecs, Ref(irlabs))
    
    return BandRepSet(sgnum, BRs, kvs, klabs, irlabs, allpaths, spinful, timereversal)
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

    irvecs = [zeros(Int64, length(irlabs)) for _ in OneTo(length(brtags))]
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

Returns the `BandRepSet` for space group `sgnum` (provided by the [Bilbao Crystallographic
Server](http://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl)) and dimension `D`.

Keyword arguments:

- `allpaths`: include a minimal sufficient set (`false`, default) or all (`true`) k-vectors. 
- `spinful`: single- (`false`, default) or double-valued (`true`) irreps, as appropriate for
             spinless and spinful particles, respectively.
- `timereversal`: assume presence (`true`, default) or absence (`false`) of time-reversal
                  symmetry.
"""
function bandreps(sgnum::Integer, D::Integer=3;
                  allpaths::Bool=false, spinful::Bool=false,
                  timereversal::Bool=true)
    D ∉ (1,2,3) && _throw_invaliddim(D)
    paths_str = allpaths ? "allpaths" : "maxpaths"
    brtype_str = timereversal ? "elementaryTR" : "elementary"
    filename = (@__DIR__)*"/../data/bandreps/$(D)d/$(brtype_str)/$(paths_str)/$(string(sgnum)).csv"
    open(filename) do io
        BRS = dlm2struct(io, sgnum, allpaths, spinful, timereversal)
    end 
end

"""
    classification(BRS::BandRepSet) --> String

Calculate the symmetry indicator classification of a band representation set, meaning the 
index-classification inferrable on the basis of symmetry alone.

Technically, the calculation answers a question like "what direct product of 
``\\mathbb{Z}_n`` groups is the the quotient group
``X^{\\text{bs}} = \\{\\text{BS}\\}/\\{\\text{AI}\\}`` isomorphic to?".

See e.g. Po, Watanabe, & Vishwanath, Nature Commun. **8**, 50 (2017) for more information.
"""
function classification(BRS::BandRepSet)
    # get the diagonal components of the Smith normal decomposition (≥ 0)
    Λ = smith(matrix(BRS)).SNF
    @assert all(≥(0), Λ)
    nontriv_idx = findall(x-> !(isone(x) || iszero(x)), Λ)
    if isempty(nontriv_idx)
        return "Z₁"
    else
        return ("Z"*join(subscriptify.(string.(sort(@view Λ[nontriv_idx]))), "×Z"))
    end
end

"""
    basisdim(BRS::BandRepSet) --> Int64

Computes the dimension of the (linearly independent parts) of a band representation set.
This is ``d^{\\text{bs}} = d^{\\text{ai}}`` in the notation of Po, Watanabe, & Vishwanath,
Nature Commun. **8**, 50 (2017). In other words, this is the number of linearly independent
basis vectors that span the expansions of a band structure or atomic insulator viewed as
symmetry-data.
""" 
function basisdim(BRS::BandRepSet)
    Λ = smith(matrix(BRS)).SNF
    nnz = count(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix
    return nnz
end

"""
    wyckbasis(BRS::BandRepSet) --> Vector{Vector{Int64}}

Computes the (band representation) basis for bands generated by localized orbitals placed at
the Wyckoff positions. Any band representation that can be expanded on this basis with
positive integer coefficients correspond to a trivial insulators (i.e. deformable to atomic
limit).
Conversely, bands that cannot are topological, either fragily (some negative coefficients)
or strongly (fractional coefficients).
"""
function wyckbasis(BRS::BandRepSet) 
    # Compute Smith normal form: for an n×m matrix B with integer elements,
    # find matrices S, diagm(Λ), and T (of size n×n, n×m, and m×m, respectively)
    # with integer elements such that B = S*diagm(Λ)*T. Λ is a vector
    # [λ₁, λ₂, ..., λᵣ, 0, 0, ..., 0] with λⱼ₊₁ divisible by λⱼ and r ≤ min(n,m).
    # The matrices T and S have integer-valued pseudo-inverses.
    F = smith(matrix(BRS)) # Smith normal factorization with λⱼ ≥ 0
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
# ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
# ==============================   Some context   ==============================
# The thing to remember here is that if B ≡ matrix(BRS), then B is in general a 
# linearly dependent basis for all atomic insulators. Suppose the symmetry data 
# vector is 𝐧, then this means that there exists a coefficient vector ("into" the
# cols of B) 𝐜, such that B𝐜 = 𝐧, where 𝐜 has strictly positive integer coefficients 
# (corresponding to "stacking" of atomic bands). To get a linearly independent basis,
# we seek the column space of B: this generates all possible 𝐧; the price we pay, 
# though, is that in general the column space of B would also allow negative 
# coefficients ("subtraction" of atomic bands). We could of course just get the 
# column space from the SVD, but that would give us non-integral coefficients in 
# general. Instead, we can use the Smith normal form noted above. Specifically,
# the column space of B is S[:,1:r] with r number of nonzero λⱼ. This column 
# space in fact generates both all atomic insulators and all fragile topological 
# insulators (effectively, additions and subtractions of columns of B).
# The trouble with distinguishing fragile and atomic insulators is that the columns of 
# B are linearly dependent; and the span of B obtained from B[:,1:r] may already 
# have included several subtractions. What we want to know is whether, for given 
# 𝐧, there exists a positive-integer coefficient solution to B𝐜 = 𝐧 or not: if we
# just find *some* negative-integer coefficient solution, then that doesn't mean
# that there couldn't also be a positive-integer coefficient solution as well, simply
# because the columns of B are linearly dependent.
# One approach could be to find some solution 𝐜′, such that B𝐜′=𝐧, which may contain
# negative coefficients. From this particular solution, we can generate all possible
# solutions by adding elements from the null space B to 𝐜′: if we can identify an 
# element from the null space to add into 𝐜′ such that the result is positive, 
# we have found solution (= an atomic insulator).
# ==============================================================================
# In searching on this topic, I came across the following KEY WORDS:
#   - Linear system of Diophantine equations: any system Ax=b with integer A and b.
#       See e.g. https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine
#   - Integer linear programming: concerned with solving integer-coefficient equations
#       subject to conditions. We might phrase our problem as a "standard form linear
#       programming" problem:
#           minimize    𝐭ᵀ𝐜
#           subject to  B𝐜 = 𝐧
#                       cᵢ ≥ 0
#                       𝐜 ∈ 𝐙ⁿ
#       where we would choose 𝐭 = (1,1,1,1, ...) with the idea of finding the solution
#       that involves the fewest band representations. Not sure this choice is actually 
#       meaningful; we already have to constrain the total number of bands involved
#       (by adding a (column to A) = dim.(BRS) and a (row to 𝐧) = the number of included
#       bands).
#       It seems like this can be done in JuMP with GLPKSolverLP(). (or via GLPK.jl)
#   - Semirings: technically, the problem is that the set of natural numbers ℕ = 0,1,...
#       is not a ring, but a semiring (no additive inverse). There may be some packages
#       out there that deal specifically with semirings? Actually, since it has an identity
#       element under addition (0) it is a monoid.
#   - The monoid paper from Bernevig: https://arxiv.org/pdf/1905.03262.pdf
#       This actually seems to be a very worthwhile starting point; they're attacking 
#       exactly this point.
# ==============================================================================



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


"""
    matching_littlegroups(BRS::BandRepSet)

Finds the matching little groups for each *k*-point referenced in `BRS`. This is mainly a 
a convenience accessor, since e.g. `littlegroup(::SpaceGroup, ::KVec)` could already give
the required little groups. The benefit here is that the resulting operator sorting of
the returned little group is identical ISOTROPY's, so we can rely on that later on.

Returns a `Vector{<:LittleGroup}` with ordering identical to that of the k-point labels in 
`BRS`. Note that this is different from the return type of `get_littlegroups` which returns
an unordered `Dict`.

## Note 1

The little groups from ISOTROPY do not include copies of operators that would be identical
when transformed to a primitive basis. The operators are, however, still given in a
conventional basis.

## Note 2

An error is thrown if a referenced little group cannot be found (currently, this can happen
for certain k-points in ``Φ-Ω``, see src/special_representation_domain_kpoints.jl)
"""
function matching_littlegroups(BRS::BandRepSet, ::Val{D}=Val(3)) where D
    lgs = get_littlegroups(num(BRS), Val(D)) # TODO: generalize to D≠3

    klabs_in_brs = klabels(BRS) # find all k-point labels in BandRepSet
    if !issubset(klabs_in_brs, keys(lgs))
        throw(DomainError(klabs_in_brs, "Could not locate all LittleGroups from BandRep"))
    end

    return getindex.(Ref(lgs), klabs_in_brs)
end


function matching_lgirreps(BRS::BandRepSet)
    lgirsd = get_lgirreps(num(BRS), Val(3))
    # create "physical/real" irreps if `BRS` assumes time-reversal symmetry
    if BRS.timeinvar 
        for (klab, lgirs) in keys(lgirsd)
            lgirsd[klab] = realify(lgirs)
        end
    end
    # all lgirreps from ISOTROPY as a flat vector; note that sorting is arbitrary
    lgirs = collect(Iterators.flatten(values(lgirsd))) # TODO: generalize to D≠3

    # find all the irreps in lgirs that feature in the BandRepSet, and 
    # sort them according to BandRepSet's sorting
    lgirlabs = label.(lgirs)
    brlabs   = normalizesubsup.(irreplabels(BRS))

    find_and_sort_idxs = Vector{Int64}(undef, length(brlabs))
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