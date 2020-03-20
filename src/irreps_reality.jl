const TEST_αβγ = [0.123,0.456,0.789] # arbitrary test numbers for KVecs
# TODO: This implementation should follow the discussion on p. 650-652 in Bradley 
#       & Cracknell's book (there's some discussion in 622-626 as well, but that's 
#       for point groups). Their discussion is for magnetic groups but is generally 
#       applicable, and is by far the most clear and thorough discussion that I've 
#       found so far.
#       Cornwell also does a good job of explicating this.
#       Inui on p. 296-299 also discuss it, but is less clear overall.
function realify(irs::AbstractVector{LGIrrep{D}}, verbose::Bool=false) where D
    Nirr = length(irs)
    kv = kvec(first(irs)) # must be the same for all irreps in list
    kv_αβγ = kv(TEST_αβγ)
    sgnum = num(first(irs))
    lgops = operations(first(irs))
    Nops = order(first(irs)) # order of little group (= # of operations)

    cntr = centering(sgnum, D)
    sgops = operations(spacegroup(sgnum, D))

    verbose && print(klabel(first(irs)), " │ ")

    # Check if -𝐤 is in the star of 𝐤, or if 𝐤 is equivalent to -𝐤: 
    # if so, TR is an element of the little group; if not, it isn't 
    # ║ 𝐑𝐞𝐚𝐬𝐨𝐧: if there is an element g of the (unitary) 𝑠𝑝𝑎𝑐𝑒 group G   
    # ║   that takes 𝐤 to -𝐤 mod 𝐆, then (denoting the TR element by Θ, 
    # ║   acting as θ𝐤 = -𝐤) the antiunitary element θg will take 𝐤 to  
    # ║   𝐤 mod 𝐆, i.e. θg will be an element of the little group of 𝐤
    # ║   M(k) associated with the **gray** space group M ≡ G + θG.
    # ║   Conversely, if no such element g exists, there can be no anti-
    # ║   unitary elements in the little group derived from M; as a result, 
    # ║   TR is not part of the little group and so does not modify its 
    # ║   small irreps (called "co-reps" for magnetic groups).
    # ║   There can then only be type 'x' degeneracy (between 𝐤 and -𝐤)
    # ║   but TR will not change the degeneracy at 𝐤 itself.
    if !isapproxin(-kv, kstar(sgops, kv, cntr), cntr; atol=DEFAULT_ATOL)
        corep_idxs = [[i] for i in Base.OneTo(Nirr)] # TR ∉ M(k) ⇒ smalls irrep (... small co-reps) not modified by TR
        verbose && println(klabel(first(irs)), "ᵢ ∀i (type x) ⇒  no additional degeneracy (star{k} ∌ -k)")

    else
        # Test if 𝐤 is equivalent to -𝐤, i.e. if 𝐤 = -𝐤 + 𝐆
        k_equiv_kv₋ = isapprox(-kv, kv, cntr; atol=DEFAULT_ATOL)

        # Find an element in G that takes 𝐤 → -𝐤 (if 𝐤 is equivalent to -𝐤, 
        # then this is just the unit-element I (if `sgops` is sorted conven-
        # tionally, with I first, this is indeed what the `findfirst(...)`  
        # bits below will find)
        if !k_equiv_kv₋
            g₋ = sgops[findfirst(g-> isapprox(g∘kv, -kv, cntr; atol=DEFAULT_ATOL), sgops)]
        else
            # This is a bit silly: if k_equiv_kv₋ = true, we will never use g₋; but I'm not sure if 
            # the compiler will figure that out, or if it will needlessly guard against missing g₋?
            g₋ = SymOperation{D}(hcat(I, zeros(D))) # ... the unit element I
        end

        # -𝐤 is part of star{𝐤}; we infer reality of irrep from ISOTROPY's data (could also 
        # be done using `herring(...)`). ⇒ deduce new small irreps (... small co-reps).
        corep_idxs = Vector{Vector{Int64}}()
        skiplist = Vector{Int64}()
        for (i, ir) in enumerate(irs)
            if i ∈ skiplist; continue; end # already matched to this irrep previously; i.e. already included now
            verbose && i ≠ 1 && print("  │ ")

            if type(ir) == 1     # real
                push!(corep_idxs, [i])
                if verbose
                    println(formatirreplabel(label(ir)), " (real) ⇒  no additional degeneracy")
                end

            elseif type(ir) == 2 # pseudo-real
                # doubles irrep on its own
                push!(corep_idxs, [i, i])
                if verbose
                    println(formatirreplabel(label(ir)^2), " (pseudo-real) ⇒  doubles degeneracy"); 
                end

            elseif type(ir) == 3 # complex
                # In this case, there must exist a "partner" irrep (say, Dⱼ) which is 
                # equal to the complex conjugate of the current irrep (say, Dᵢ); we 
                # next search for this equivalence.
                # When we check for equivalence between irreps Dᵢ* and Dⱼ we must
                # account for the possibility of a 𝐤-dependence in the matrix-form
                # of the irreps; specifically, for an element g, its small irrep is
                #     Dᵢ[g] = exp(2πik⋅τᵢ[g])Pᵢ[g],
                # where, crucially, for symmetry lines, planes, and general points
                # 𝐤 depends on (one, two, and three) free parameters (α,β,γ).
                # Thus, for equivalence of irreps Dᵢ* and Dⱼ we require that
                #     Dᵢ*[g] ~ Dⱼ[g]       ∀g ∈ G(k)
                #  ⇔ exp(-2πik⋅τᵢ[g])Pᵢ*[g] ~ exp(2πik⋅τⱼ[g])Pⱼ[g]
                # It seems rather tedious to prove that this is the case for all 𝐤s
                # along a line/plane (α,β,γ). Rather than attempt this, we simply test
                # against an arbitrary value of (α,β,γ) [superfluous entires are ignored]
                # that is non-special (i.e. not ={0,0.5,1}); this is `TEST_αβγ`.

                # Characters of the conjugate of Dᵢ, i.e. tr(Dᵢ*) = tr(Dᵢ)*
                θχᵢ = conj.(tr.(irreps(ir, TEST_αβγ))) 
                
                # Find matching complex partner
                partner = 0
                for j = i+1:Nirr
                    if j ∉ skiplist && type(irs[j]) == 3 # only check if j has not previously matched; 
                                                         # similarly, only check if the jth irrep is complex.

                        # Note that we require only equivalence of Dᵢ* and Dⱼ; not equality. 
                        # Cornwell describes (p. 152-153 & 188) a neat trick for checking this 
                        # efficiently: specifically, Dᵢ* and Dⱼ are equivalent irreps if 
                        #     χⁱ(g)* = χʲ(g₋⁻¹gg₋) ∀g ∈ G(k)
                        # with g₋ an element of G that takes 𝐤 to -𝐤, and where χⁱ (χʲ) denotes
                        # the characters of the respective irreps.
                        χⱼ = tr.(irreps(irs[j], TEST_αβγ))
                        match = true
                        for n in Base.OneTo(Nops)
                            if k_equiv_kv₋ # 𝐤 = -𝐤 + 𝐆 ⇒ g₋ = I (the unit element), s.t. g₋⁻¹gg₋ = I⁻¹gI = g
                                χⱼ_g₋⁻¹gg₋ = χⱼ[n]
                            else           # 𝐤 not equivalent to -𝐤, i.e. 𝐤 ≠ -𝐤 + 𝐆
                                g₋⁻¹gg₋ = compose(compose(inv(g₋), lgops[n], false), g₋, false)
                                n′, Δw = findequiv(g₋⁻¹gg₋, lgops, cntr)
                                χⱼ_g₋⁻¹gg₋ = cis(2π*dot(kv_αβγ, Δw)) .* χⱼ[n′] # cis(x) = exp(ix)
                            end
                            
                            match = isapprox(θχᵢ[n], χⱼ_g₋⁻¹gg₋; atol=DEFAULT_ATOL)
                            if !match # ⇒ not a match
                                break
                            end
                        end

                        if match # ⇒ a match
                            partner = j
                            if verbose; 
                                println(formatirreplabel(label(ir)*label(irs[j])), " (complex) ⇒  doubles degeneracy")
                            end
                        end
                    end
                end
                partner === 0 && throw(ErrorException("Didn't find a matching complex partner for $(label(ir))"))
                push!(skiplist, partner)

                push!(corep_idxs, [i, partner])
                
            else
                throw(ArgumentError("Invalid real/pseudo-real/complex type = $(type(ir))"))
            end
        end
    end

    Ncoreps = length(corep_idxs)

    # New small co-rep labels (composite)
    newlabs = Tuple(join(label(irs[i]) for i in corep_idxs[i′]) for i′ in Base.OneTo(Ncoreps))

    # TODO: New small irreps (small co-reps)
    #=
    for i′ in Base.OneTo(Ncoreps)
        idxs = coreps_idxs[i′]
        if length(idxs) == 1      # real or type x
            # same as before
        elseif idxs[1] == idxs[2] # pseudoreal 
            # doubles self
        else                      # complex
            # doubles with complex conjugate
            # what to do about exp(ikτ) dependence? Need new type, different from LGIrrep?
            # maybe the τ values are the same? Could just check...
        end
    end
    =#
    return corep_idxs, newlabs
end


"""
    herring(ir::LGIrrep, sgops::AbstractVector{SymOperation{D}},
            αβγ::Union{Vector{<:Real},Nothing}=nothing)        --> Tuple{Int, Int}

Computes the Herring criterion for a little group irrep `ir`, from 

    [∑ χ({β|b}²)]/[g₀/M(k)] 

over symmetry operations {β,b} that take k → -k. Here g₀ is the order of the point group
of the space group and M(k) is the order of the star(k) [both in a primitive basis].

The returned value, [∑ χ({β|b}²)]/[g₀/M(k)], is one of three integers in {1,-1,0} 
corresponding to {real, pseudoreal, complex} reality. We remind that ISOTROPY's indication
of the same reality types i {1,2,3}.

The provided space group operations `sgops` **must** be the set reduced by 
primitive translation vectors; i.e. using `spacegroup(...)` directly is **not** 
allowable in general. Using the operations from the Γ point of ISOTROPY's 
dataset is, however, fine.

As a sanity check, a value of `αβγ` can be provided to check for invariance
along a symmetry line/plane/general point in k-space. Obviously, the reality 
type should invariant to this choice.

**Implementation:** 
See e.g. Inui's Eq. (13.48), Dresselhaus, p. 618, and 
and Herring's original paper at https://doi.org/10.1103/PhysRev.52.361.
We mainly followed Cornwell, p. 150-152 & 187-188.
"""
function herring(ir::LGIrrep, sgops::AbstractVector{SymOperation{D}}, αβγ::Union{Vector{<:Real},Nothing}=nothing) where D

    lgops = operations(ir)
    kv = kvec(ir)
    kv₋ = -kv
    cntr = centering(num(ir), D)
    Ds = irreps(ir, αβγ) # irrep matrices
    kv_αβγ = kv(αβγ)

    s = zero(ComplexF64)
    for op in sgops
        if isapprox(op∘kv, kv₋, cntr, atol=DEFAULT_ATOL) # check if op∘k == -k; if so, include in sum
            op² = compose(op, op, false) # this is op∘op, _including_ trivial lattice translation parts
            # find the equivalent of `op²` in `lgops`; this may differ by a number of 
            # primitive lattice vectors `w_op²`; the difference must be included when 
            # we calculate the trace of the irrep 𝐃: the irrep matrix 𝐃 is ∝exp(2πi𝐤⋅𝐭)
            idx_of_op²_in_lgops, Δw_op² = findequiv(op², lgops, cntr)
            ϕ_op² = cis(2π*dot(kv_αβγ, Δw_op²)) # phase accumulated by "trivial" lattice translation parts [cis(x) = exp(ix)]
            χ_op² = ϕ_op²*tr(Ds[idx_of_op²_in_lgops]) # χ(op²)

            s += χ_op²
        end
    end

    pgops = pointgroup(sgops) # point group assoc. w/ space group
    g₀ = length(pgops) # order of pgops (denoted h, or macroscopic order, in Bradley & Cracknell)
    Mk = length(kstar(pgops, kv, cntr)) # order of star of k (denoted qₖ in Bradley & Cracknell)
    normalization = round(Int, g₀/Mk) # order of G₀ᵏ; the point group derived from the little group Gᵏ (denoted b in Bradley & Cracknell; [𝐤] in Inui)
    if !isapprox(normalization, g₀/Mk)
        throw(ErrorException("The little group is not factored by its point group and"*
                             "star{k}: this should never happen"))
    end

    # check that output is a real integer and then convert to that for output...
    if norm(imag(s)) < DEFAULT_ATOL 
        sInt = round(Int,real(s)); 
    else 
        throw(error("Herring criterion should yield a real value; obtained complex s=$(s)")) 
    end
    if norm(sInt-real(s)) > DEFAULT_ATOL 
        throw(error("Herring criterion should yield an integer; obtained s=$(s)"))
    end

    # sInt = ∑ χ({β|b}²) and normalization = g₀/M(k) in Cornwell's Eq. (7.18) notation
    herring_type = Int64(sInt/normalization)
    if herring_type ∉ (0,1,-1)
        throw(DomainError(herring_type, "Calculation of the Herring criterion incorrectly"*
                                        "produced a value ∉ (0,1,-1)"))
    end

    return Int64(sInt/normalization) # return [∑ χ({β|b}²)]/[g₀/M(k)]
end