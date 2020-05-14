using Crystalline

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
                    # TODO (more complicated than I thought: have to match across different special lgirreps)
                end 
            end
        end
    end
end