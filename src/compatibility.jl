using SGOps

"""
    subduction_count(Dᴳᵢ, Dᴴⱼ[, αβγᴴⱼ]) --> Int64

For two groups G and H, where H is a subgroup of G, i.e. G>H, with associated 
irreducible representations `Dᴳᵢ`(g) and `Dᴴⱼ`(h) for elements g∈G and h∈H<G, compute 
the compatibility relation between the two irreps from the subduction reduction 
formula (or "magic" formula/Schur orthogonality relation): this is essentially
how many times `nᴳᴴᵢⱼ` the subduced representation `Dᴳᵢ`↓H contains the irrep `Dᴴⱼ`; 
in other words, this gives the compatibility between the two irreps.

Optionally, a vector `αβγᴴⱼ` may be provided, to evaluate the characters/irreps 
of Dᴳᵢ at a concrete value of αβγ. This is e.g. meaningful for LGIrreps at non-
special k-vectors. Defaults to `nothing`.

The reduction formula [e.g. Eq. (15) of https://arxiv.org/pdf/1706.09272.pdf] is:

        nᴳᴴᵢⱼ = |H|⁻¹∑₍ₕ₎ χᴳᵢ(h)χᴴⱼ(h)*

As an example, consider space group 135 and the two compatible k-vectors 
Γ (a point) and Σ (a plane):
```
    lgirvec = get_lgirreps(135, Val(3))
    Γ_lgirs = lgirvec[1] # at Γ ≡ [0.0, 0.0, 0.0]
    Σ_lgirs = lgirvec[4] # at Σ ≡ [α, α, 0.0]
```
We can test their compatibility like so:
```
    [[subduction_count(Γi, Σj) for Γi in Γ_lgirs for Σj in Σ_lgirs]
    > # Γ₁ Γ₂ Γ₃ Γ₄ Γ₅
    >  [ 1, 0, 1, 1, 2] # Σ₁
    >  [ 0, 1, 1, 2, 1] # Σ₂
```
This entails the following compatibility relations between irreps at Γ and Σ:

        Γ₁ → Σ₁           degeneracies: 1 → 1
        Γ₂ → Σ₂                         1 → 1
        Γ₃ → Σ₁ + Σ₂                    2 → 1 + 1
        Γ₄ → Σ₁ + 2Σ₂                   3 → 1 + 2
        Γ₅ → 2Σ₁ + Σ₂                   3 → 2 + 1

where, in this case, all the small irreps are one-dimensional.
"""
function subduction_count(Dᴳᵢ::T, Dᴴⱼ::T, 
                          αβγᴴⱼ::Union{Vector{<:Real},Nothing}=nothing) where T<:AbstractIrrep
    # find matching operations between H & G and verify that H<G 
    boolsubgroup, idxsᴳ²ᴴ = _findsubgroup(operations(Dᴳᵢ), operations(Dᴴⱼ))
    !boolsubgroup && throw(DomainError("Provided irreps are not H<G subgroups"))

    # compute characters 
    # TODO: Care should be taken that the irreps 
    # actually can refer to identical k-points; that should be a check 
    # too, and then we should make sure that the characters are actually
    # evaluated at that KVec
    χᴳᵢ = characters(Dᴳᵢ)
    χᴴⱼ = characters(Dᴴⱼ, αβγᴴⱼ)

    # compute number of times that Dᴴⱼ occurs in the reducible 
    # subduced irrep Dᴳᵢ↓H
    s = zero(ComplexF64)
    @inbounds for (idxᴴ, χᴴⱼ′) in enumerate(χᴴⱼ)
        s += χᴳᵢ[idxsᴳ²ᴴ[idxᴴ]]*conj(χᴴⱼ′)
    end
    (abs(imag(s)) > DEFAULT_ATOL) && throw("unexpected finite imaginary part")
    nᴳᴴᵢⱼ_float = real(s)/order(Dᴴⱼ)
    nᴳᴴᵢⱼ = round(Int64, nᴳᴴᵢⱼ_float)
    abs(nᴳᴴᵢⱼ - nᴳᴴᵢⱼ_float) > DEFAULT_ATOL && throw("unexpected non-integral compatibility count")
    
    return nᴳᴴᵢⱼ
end

"""
    find_compatible_kvec(kv::KVec, kvs′::Vector{KVec})
"""
function find_compatible_kvec(kv::KVec, kvs′::Vector{KVec})
    !isspecial(kv) && throw(DomainError(kv, "input kv must be a special k-point"))

    compat_idxs = Vector{Int64}()
    compat_αβγs = Vector{Vector{Float64}}()
    @inbounds for (idx′, kv′) in enumerate(kvs′)
        isspecial(kv′) && continue # must be a line/plane/general point to match a special point kv
        compat_bool, αβγ′ = is_compatible_kvec(kv, kv′)
        if compat_bool
            push!(compat_idxs, idx′)
            push!(compat_αβγs, αβγ′)
        end
    end

    return compat_idxs, compat_αβγs
end

function is_compatible_kvec(kv::KVec, kv′::KVec)
    # TODO: I think we need to do this in the primitive basis! But it is nontrivial, since
    #       if we match k-points across a G-vector, we also need to transform the irrep
    #       with a suitable phase factor.

    # TODO: this cannot treat finding a compatible plane to a line
    k₀, _  = parts(kv) 
    k₀′, kabc′ = parts(kv′)

    # least squares solve via QR factorization; equivalent to pinv(kabc)*(k₀-k₀′) but faster
    αβγ′ = qr(kabc′, Val(true))\(k₀-k₀′)  
    k′ = k₀′ + kabc′*αβγ′
    # check if least squares solution actually is a solution
    compat_bool = isapprox(k₀, k′, atol=DEFAULT_ATOL) 

    return compat_bool, αβγ′
end

"""
    compatibility(lgirvec)
"""
function compatibility(lgirvec::AbstractVector{<:AbstractVector{LGIrrep{D}}}) where D
    kvs   = kvec.(first.(lgirvec))
    klabs = klabel.(first.(lgirvec))
    Nk    = length(kvs)
    
    # prepare a graph for the connections between k-vectors
    kgraph = MetaDiGraph(Nk)
    foreach((i,kv,kl)->set_props!(kgraph, i, Dict(:kvec=>kv, :klab=>kl)), eachindex(kvs), kvs, klabs)

    for (kidxᴳ,lgirs) in enumerate(lgirvec)                 # parent group 
        kvᴳ = kvs[kidxᴳ]
        !isspecial(kvᴳ) && continue # starting point is always a special k-point
        compat_idxs, compat_αβγs = find_compatible_kvec(kvᴳ, kvs)
        for (kidxᴴ, αβγᴴ) in zip(compat_idxs, compat_αβγs)  # subgroup
            add_edge!(kgraph, kidxᴳ, kidxᴴ)
            for (iᴳ, Dᴳᵢ) in enumerate(lgirs)
                for (jᴴ, Dᴴⱼ) in enumerate(lgirvec[kidxᴴ])
                    nᴳᴴᵢⱼ = subduction_count(Dᴳᵢ, Dᴴⱼ, αβγᴴ)
                    if !iszero(nᴳᴴᵢⱼ) # add an edge between irreps Dᴳᵢ and Dᴴⱼ
                        add_edge!()
                    end
                end
            end
        end
    end
    return kgraph
end


"""
    connectivity(lgirvec)
"""
function connectivity(lgirvec::AbstractVector{<:AbstractVector{LGIrrep{D}}}) where D
    kvs   = kvec.(first.(lgirvec))
    klabs = klabel.(first.(lgirvec))
    Nk    = length(kvs)
    
    # prepare a graph for the connections between k-vectors
    kgraph = MetaDiGraph(Nk)
    foreach((i,kv,kl)->set_props!(kgraph, i, Dict(:kvec=>kv, :klab=>kl)), eachindex(kvs), kvs, klabs)

    Nspecial = 0
    @inbounds for (kidxᴳ,lgirs) in enumerate(lgirvec)       # parent group 
        kvᴳ = kvs[kidxᴳ]
        if isspecial(kvᴳ)
            Nspecial += 1
        else
            continue # starting point is always a special k-point
        end
        compat_idxs, compat_αβγs = find_compatible_kvec(kvᴳ, kvs)
        for (kidxᴴ, αβγᴴ) in zip(compat_idxs, compat_αβγs)  # subgroup
            add_edge!(kgraph, kidxᴳ, kidxᴴ)
        end
    end

    cgraph = MetaGraph(Nspecial) # connectivity graph for special k-vecs
    local_kidx¹ = 0
    @inbounds for kidx¹ in eachindex(lgirvec)
        isspecial(kvs[kidx¹]) || continue      # only compare special vectors
        local_kidx¹ += 1
        set_props!(cgraph, local_kidx¹, Dict(:kvec=>kvs[kidx¹], 
                                             :klab=>klabs[kidx¹], 
                                             :kidx=>kidx¹)) 
        local_kidx² = 0
        for kidx² in eachindex(lgirvec)
            isspecial(kvs[kidx²]) || continue  # only compare special vectors
            local_kidx² += 1
            kidx¹≥kidx² && continue            # avoid double & self-comparisons

            nbs = common_neighbors(kgraph, kidx¹, kidx²)
            for (nbidx, nb) in enumerate(nbs)
                # if the neighbor is just the general point Ω≡[α,β,γ], 
                # we don't consider the two vectors connected
                if kvs[nb] == KVec(zeros(D), Matrix{Float64}(I, D, D))
                    deleteat!(nbs, nbidx)
                    break
                end      
            end
            isempty(nbs) && continue # Ω is only connecting edge (trivial case)
            add_edge!(cgraph, local_kidx¹, local_kidx²) 
            set_props!(cgraph, Edge(local_kidx¹, local_kidx²), 
                               Dict(:klabs=>klabs[nbs],
                                    :kvecs=>kvs[nbs],
                                    :kidxs=>nbs)
                      )
        end
    end          
    return cgraph, kgraph
end


function compatibility_matrix(BRS::BandRepSet)
    lgirs_in, lgirs_out = matching_lgirreps(BRS::BandRepSet)
    for (iᴳ, Dᴳᵢ) in enumerate(lgirs_in)         # super groups
        for (jᴴ, Dᴴⱼ) in enumerate(lgirs_out)    # sub groups
            # we ought to only check this on a per-kvec basis instead of 
            # on a per-lgir basis to avoid redunant checks, but can't be asked...
            compat_bool, αβγ′ = is_compatible_kvec(kvec(Dᴳᵢ), kvec(Dᴴⱼ))
            if compat_bool
                nᴳᴴᵢⱼ = subduction_count(Dᴳᵢ, Dᴴⱼ, αβγ′)
                if !iszero(nᴳᴴᵢⱼ)
                    # TODO: more complicated than I thought: have to match across different special lgirreps
                end 
            end
        end
    end
end

"""
    compatibility_matrix(lgirsvec)

Compute all compatibility relations between lines and connectible points in k-space, each
forming a row in a matrix. Also include all filling constraints (i.e. enforcing that a band
must have the same number of states at every k-point). This is essentially the approach 
described in Song, Zhang, & Fang PRX 8, 031069 (2018)

TODO: Unfortunately, this doesn't appear to work generally yet, as can be verified by 
      comparing with the basis obtained from wyckbasis(bandrep(...)), e.g. as in 
    
    for sgnum in 1:230
        v = (size(SGOps.wyckbasis(bandreps(sgnum))[1],1), 
             size(SGOps.compatibility_basis(get_lgirreps(sgnum)),2))
        println(sgnum, ": ", v[1] == v[2], " ", v)
    end
"""
function compatibility_matrix(lgirsvec)
    kvs    = kvec.(first.(lgirsvec))
    klabs  = klabel.(first.(lgirsvec))
    irlabs = [label(lgirs) for lgirvec in lgirsvec for lgirs in lgirvec]
    #println(irlabs)
    #pop!(irlabs) # remove Γ point (always last)
    #println(irlabs)
    Nk     = length(kvs)
    Nirr   = length(irlabs)
    buffer = zeros(Int, Nirr) # preallocated constraint relation buffer
    Crels  = Vector{Vector{Int}}()
    for (kidxᴴ, Dᴴᵢs) in enumerate(lgirsvec)     # subgroup 𝐤 (lower symmmetry)
        kvᴴ = kvs[kidxᴴ]
        nᵅᵝᵞᴴ = nfreeparams(kvᴴ)
        # we let kvᴴ denote a "lower symmetry" k-line, -plane, or -volume, and want
        # to initially determine how many "higher symmetry" k-point (kvᴳ) it connects to 
        nᵅᵝᵞᴴ ≠ 0 || continue
        for (kidxᴳ, Dᴳᵢs) in enumerate(lgirsvec) # supergroup 𝐤 (higher symmmetry)
            kidxᴳ == kidxᴴ && continue
            kvᴳ = kvs[kidxᴳ]
            nᵅᵝᵞᴳ = nfreeparams(kvᴳ) 
            nᵅᵝᵞᴳ == 0 || continue # we restrict kvᴳ to high symmetry _points_
            # specifically, not all planes (of little group H) that can be connected to a 
            # line (of group G) actually have the necessary relationship H<G
            
            compat_bool, αβγᴴ = is_compatible_kvec(kvᴳ, kvᴴ)
            compat_bool || continue # kvᴴ and kvᴳ must be compatible
            
            subducts = zeros(Int, length(Dᴳᵢs)) # subduction mappings from Dᴳʲ to Dᴴ
            Dᴳᵢidxs = (:)(findfirst(==(label(first(Dᴳᵢs))), irlabs),
                          findfirst(==(label(last(Dᴳᵢs))),  irlabs))
            for Dᴴᵢ in Dᴴᵢs                                     # subgroup irrep
                Dᴴᵢidx = findfirst(==(label(Dᴴᵢ)), irlabs)
                for (j, Dᴳᵢ) in enumerate(Dᴳᵢs)                 # supergroup irrep
                    # how many times does Dᴳ subduce into Dᴴ
                    subducts[j] = subduction_count(Dᴳᵢ, Dᴴᵢ, αβγᴴ)
                end
                # `subducts` now give the a different ways that Dᴳs can subduce into Dᴴ; to 
                # get a conserved quantity, we look at the "sum" all these "options". we
                # have to take special care to allow cases where subduction counts differ:
                # As an example, consider the fictious scenario where irreps at Γ subduce
                # into a particular line irrep Λ₁ 
                #   Γ₁ → Λ₁,   Γ₂ → 2Λ₁,   Γ₃ → 3Λ₁   ⇐ (equiv. to `subducts = [1 2 3]`)
                # In this case, the overall "conservation law" is 
                #   6n(Γ₁) + 3n(Γ₂) + 2n(Γ₃) = 6n(Λ₁) 
                # with n(Dᵏ) denoting the number of times irrep Dᵏ occurs at k. We find the 
                # appropriate prefactors by taking the least common multiple of `subducts`
                max_subduct = lcm(filter(!iszero, subducts))
                for (j,s) in enumerate(subducts)
                    buffer[Dᴳᵢidxs[j]] = (s == 0 ? 0 : Int(max_subduct/s))
                end
                buffer[Dᴴᵢidx]  = -max_subduct

                push!(Crels, copy(buffer))  # write new compatibility relation
                buffer .= 0                 # reset buffer before next irrep from Dᴴᵢs
            end
        end
    end
    
    # TODO: Don't really need to add this equality constraint if the points are connectible
    for (kidxᴳ, Dᴳᵢs) in enumerate(lgirsvec)         # supergroup (higher symmmetry)
        # for all-kpoints that are not connectible, we add the constraints that there must 
        # be equally many states at each k-point
        Dᴳᵢidx_min = findfirst(==(label(first(Dᴳᵢs))), irlabs)
        Dᴳᵢidx_max = findfirst(==(label(last(Dᴳᵢs))), irlabs)
        dimsᴳᵢ = irdim.(Dᴳᵢs)
        nstatesᴳ = lcm(dimsᴳᵢ)
        nirrᴳ = Int64.(nstatesᴳ./dimsᴳᵢ) # same number of states (=nstatesᴳ) across irreps at kᴳ      
        for kidxᴴ in kidxᴳ+1:Nk # no need to include redunant combinations
            buffer .= 0 # reset buffer: new constraint 
            buffer[Dᴳᵢidx_min:Dᴳᵢidx_max] .= nirrᴳ

            #isspecial(kvs[kidxᴴ]) && continue

            Dᴴᵢs = lgirsvec[kidxᴴ]
            Dᴴᵢidx_min = findfirst(==(label(first(Dᴴᵢs))), irlabs)
            Dᴴᵢidx_max = findfirst(==(label(last(Dᴴᵢs))), irlabs)

            dimsᴴᵢ = irdim.(Dᴴᵢs)
            nstatesᴴ = lcm(dimsᴴᵢ)
            nirrᴴ = Int64.(nstatesᴴ./dimsᴴᵢ) # same number of states (=nstatesᴴ) across irreps at kᴴ
            buffer[Dᴴᵢidx_min:Dᴴᵢidx_max] .-= nirrᴴ

            if nstatesᴳ ≠ nstatesᴴ # rebalance if there's an unequal number of states in kᴳ and kᴴ, 
                fᴳᴴ = lcm(nstatesᴴ, nstatesᴳ)
                buffer[Dᴳᵢidx_min:Dᴳᵢidx_max] .*= fᴳᴴ/nstatesᴳ
                buffer[Dᴴᵢidx_min:Dᴴᵢidx_max] .*= fᴳᴴ/nstatesᴴ
            end
            push!(Crels, copy(buffer))
        end
    end
    #pretty_table(vcat(Crels'...), irlabs)

    return vcat(Crels'...) # the action of C≡vcat(Crels'...) on a valid symmetry vector is a zero vector.
end

"""
    compatibility_basis(lgirsvec)

Compute a basis for the space of band structures allowed by compatibility relations and
"filling" constraints. Return a matrix whose columns give an integer-span of all physically 
realizable band structures with symmetry content in `lgirsvec` (a vector of vector of 
`LGIrrep`s, indexed across k-points and distinct irreps). 
The rows of the basis are the distinct irrep labels, ordered [Γ₁, Γ₂, ..., M₁, ..., Ω].

This basis should span the same space as `wyckbasis(bandreps(num(lgirsvec)))[1]'` (provided
that `lgirsvec` feature all the necessary k-points). Note that the irrep-content may not be 
the same in the two, so a comparison should make sure to project out non-shared irreps.

TODO: This doesn't seem to agree with Song, Zhang, & Fang PRX 8, 031069 (2018) Table II,
      nor with SGOps.wyckbasis(bandreps(10))[1]'. Test for SG 10.
"""
function compatibility_basis(lgirsvec)
    C = compatibility_matrix(lgirsvec)

    # Get nullspace of C with integer coefficients; i.e. find the nullspace in a field of 
    # integers ℤ. 
    # Using the Smith Normal Form: For A = SΛT, we can obtain the nullspace of A from the
    # last n columns of T⁻¹ with n denoting the number of zeros in Λ [i.e. n=nullity(A); 
    # contrast this with r=rank(A). See e.g. https://core.ac.uk/download/pdf/82343294.pdf 
    # regarding the Smith normal form and its application to null spaces. See also 
    # scripts/derive_sg2_bandrep.jl
    F = SmithNormalForm.smith(C) # Smith Normal Form
    T⁻¹, Λ = F.Tinv, F.SNF
    r = sum(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix = rank(C)
    zidxs  = r+1:length(Λ)
    basis = T⁻¹[:, zidxs] # the columns of T⁻¹ are the new basis

    return basis
end