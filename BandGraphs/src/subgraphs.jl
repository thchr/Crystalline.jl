function build_subgraphs(
            n::AbstractSymmetryVector{D},
            subductions,
            lgirsd::Dict{String, <:AbstractVector{LGIrrep{D}}}
            ) where D

    seen_irs = 0 # current number of irreps aggregated across all partitions
    kidx = 0 # number of distinct k-points we've seen

    # maximal k-point partitions
    partitions_max = Partition{D}[]
    for (klab, nsₖ, lgirsₖ) in zip(klabels(n), multiplicities(n), irreps(n))
        kidx += 1
        p_max = build_partition(klab, nsₖ, lgirsₖ, kidx, seen_irs)
        push!(partitions_max, p_max)
        seen_irs = last(p_max.iridxs)
    end

    # monodromy additions (same irrep population as their "non-monodromy"/"original" irrep
    # variants, but connect differently to non-maximal k-manifolds); PRE 96, 023310 (2017)
    monodromy_pairs = Tuple{Int, Int}[]
    for s in subductions
        s.monodromy || continue
        klab_max′ = string(s.c.kᴳ.label) # monodromy k-label (with a '′')
        klab_max = rstrip(klab_max′, '′') # "original" k-label (no '′')

        # skip if we already added info to `partitions_max`
        any(p->p.klab==klab_max′, partitions_max) && continue
        
        # now we have to find the entries of k′ ("monodromy" point) in k ("original" point)
        idx = findfirst(==(klab_max), klabels(n))
        if isnothing(idx)
            # user might have intentionally left out parts of a symmetry vector; then this
            # could fail to find the associated "original" irrep in which case it makes
            # best sense to leave it out; still, warn just in case since this should be very
            # rare
            @warn "failed to add monodromy irrep data at $klab_max′ since $klab_max is not included in symmetry data"
            continue
        else
            kidx += 1
        end
        @assert replace.(s.irlabsᴳ, Ref('′'=>"")) == label.(irreps(n)[idx]) # verify assumption used below
        
        lgirs_max′ = map(enumerate(irreps(n)[idx])) do (i,lgir)
            LGIrrep{D}(s.irlabsᴳ[i], lgir.g, lgir.matrices, lgir.translations, lgir.reality,
                       lgir.iscorep)
        end
        p_max = build_partition(klab_max′, multiplicities(n)[idx], lgirs_max′, kidx, seen_irs)
        push!(partitions_max, p_max)
        seen_irs = last(p_max.iridxs)

        kidx_original = something(findfirst(p->p.klab==klab_max, partitions_max))
        push!(monodromy_pairs, (kidx_original, kidx))
    end

    # non-maximal k-point partitions
    partitions_nonmax = Partition{D}[]
    for p_max in partitions_max # row index (maximal k)
        klab_max = p_max.klab
        idxs = findall(s -> string(s.c.kᴳ.label) == klab_max, subductions)
        subductions_from_max = @view subductions[idxs]
        for sᴳᴴ in subductions_from_max # subductions from maximal to non-maximal manifolds
            klab_nonmax = string(sᴳᴴ.c.kᴴ.label)
            if klab_nonmax ∈ getfield.(partitions_nonmax, Ref(:klab))
                continue # only need to add partition once
            else
                kidx += 1
            end

            # first, find and aggregate all the non-maximal irrep-labels as nodes, without
            # taking any care about the order of them in `irlabs`
            lgirs = LGIrrep{D}[]
            for lgir_max in p_max.lgirs
                irlab_max = label(lgir_max)
                idx_max = something(findfirst(==(irlab_max), sᴳᴴ.irlabsᴳ)) # in subduction irrep sorting
                for (nᴴ, irlabᴴ) in zip(sᴳᴴ.table[idx_max,:], sᴳᴴ.irlabsᴴ)
                    if nᴴ == 0
                        continue
                    else
                        lgirs_nonmax = lgirsd[klab_nonmax]
                        iridxᴴ = findfirst(lgir->label(lgir)==irlabᴴ, lgirs_nonmax)
                        if isnothing(iridxᴴ)
                            _throw_failed_to_find_irrep(irlabᴴ, lgirs_nonmax)
                        else
                            append!(lgirs,
                                    Iterators.repeated(lgirs_nonmax[something(iridxᴴ)], nᴴ))
                        end
                    end
                end
            end
            # now, the elements in `lgirs` are not sorted: restore sorting and compute
            # `multiples` and `global_iridxs` as well
            sort!(lgirs, by = label)
            local_iridx = 1
            multiples = UnitRange{Int}[]
            while local_iridx ≤ length(lgirs)
                local_irlab = label(lgirs[local_iridx])
                nₖᵢ = count(lgir -> label(lgir)==local_irlab, lgirs)
                push!(multiples, local_iridx:local_iridx+nₖᵢ-1)
                local_iridx += nₖᵢ
            end
            global_iridxs = seen_irs .+ (1:length(lgirs)) # indices of irrep nodes in "global sorting"
            seen_irs = last(global_iridxs)

            # now we have all the node info about the non-maximal irreps at `klab_nonmax`
            p_nonmax = Partition(klab_nonmax, lgirs, multiples, false, kidx, global_iridxs)
            push!(partitions_nonmax, p_nonmax)
        end
    end
   
    # TODO: sort `partitions_max` and `partitions_nonmax` so that we benefit maximally from 
    #       from the pinning order when we consider permutations: basically, want to sort
    #       in such a way that subgraphs with most permutations get pinned
    
    # build all the subgraphs that make up the full graph
    subgraphs = SubGraph{D}[]
    for p_max in partitions_max # maximal k-manifolds
        # i = p_max.kidx
        klab_max  = p_max.klab
        lgirs_max = p_max.lgirs
        for p_nonmax in partitions_nonmax # non-maximal k-manifolds
            # j = p_nonmax.kidx
            klab_nonmax  = p_nonmax.klab
            lgirs_nonmax = p_nonmax.lgirs

            idx = findfirst(subductions) do s
                string(s.c.kᴳ.label) == klab_max && string(s.c.kᴴ.label) == klab_nonmax
            end
            if !isnothing(idx)
                Aᵢⱼ = zeros(Int, length(lgirs_max), length(lgirs_nonmax))
                # pick the "most diagonal" connection as our canonical, unpermuted matrix
                s = subductions[idx]
                for (q, lgir_max) in enumerate(lgirs_max)
                    q′ = something(findfirst(==(label(lgir_max)), s.irlabsᴳ))
                    for (r′,v) in enumerate(@view s.table[q′,:])
                        irlab_nonmax = s.irlabsᴴ[r′]
                        r = 0
                        for _ in 1:v
                            while true
                                r = something(findnext(lgir->label(lgir)==irlab_nonmax,
                                              lgirs_nonmax, r+1))
                                if iszero(Aᵢⱼ[:,r])
                                    Aᵢⱼ[q,r] = Crystalline.irdim(lgirs_nonmax[r])
                                    break
                                end
                            end
                        end
                    end
                end
                push!(subgraphs, SubGraph(p_max, p_nonmax, Aᵢⱼ))
            end
        end
    end

    # if we had any monodromy additions, we also add a non-maximal partition associated
    # with the generic point Ω and its trivial irrep Ω₁, in order to tie together every
    # pair of monodromy-related irreps via Ω₁: we do this to ensure that band graphs
    # which are connected in energy are also connected as graphs (via Ω₁, at least)
    μ = occupation(n)
    if !isempty(monodromy_pairs)
        lgir_Ω = only(lgirsd["Ω"])
        lg_Ω = group(lgir_Ω)
        for (original_kidx, monodromy_kidx) in monodromy_pairs
            original_p_max  = partitions_max[original_kidx]
            monodromy_p_max = partitions_max[monodromy_kidx]

            # prepare to add new Ω-partition to graph; but one that is specific to the
            # connection between the two monodromy points
            klab_Ω′ = "Ω" * klabel(original_p_max.klab)
            irlab_Ω′ = klab_Ω′ * "₁"
            
            kidx += 1
            global_iridxs = seen_irs .+ (1:occupation(n))
            seen_irs = last(global_iridxs)

            lg_Ω′ = LittleGroup{D}(lg_Ω.num, lg_Ω.kv, klab_Ω′, lg_Ω.operations)
            lgir_Ω′ = LGIrrep{D}(irlab_Ω′, lg_Ω′, lgir_Ω.matrices, lgir_Ω.translations,
                                 lgir_Ω.reality, lgir_Ω.iscorep)
            p_Ω′ = Partition(klab_Ω′, fill(lgir_Ω′, μ), [1:μ], false, kidx, global_iridxs)
            push!(partitions_nonmax, p_Ω′)
        
            A = zeros(Int, length(original_p_max.lgirs), μ)
            col_start = 1
            for (row, lgir) in enumerate(original_p_max.lgirs)
                cols = col_start:(col_start+irdim(lgir)-1)
                A[row,cols] .= 1
                col_start = last(cols)+1
            end
            push!(subgraphs, SubGraph(original_p_max,  p_Ω′, A, #=pinned=# true))
            push!(subgraphs, SubGraph(monodromy_p_max, p_Ω′, A, #=pinned=# true))
        end
    end

    partitions = vcat(partitions_max, partitions_nonmax)

    return BandGraph(subgraphs, partitions)
end

 function build_partition(klab, nsₖ, lgirsₖ :: AbstractVector{LGIrrep{D}}, kidx, seen_irs) where D
    lgirs = LGIrrep{D}[]
    multiples = UnitRange{Int}[]
    max_local_iridx = 0
    for (nₖᵢ, lgirₖᵢ) in zip(nsₖ, lgirsₖ)
        nₖᵢ == 0 && continue # irrep not included, nothing to add
        append!(lgirs, Iterators.repeated(lgirₖᵢ, nₖᵢ))
        push!(multiples, (max_local_iridx+1):length(lgirs))
        max_local_iridx = length(lgirs)
    end
    global_iridxs = seen_irs .+ (1:max_local_iridx) # indices of irrep nodes in "global sorting"
    return Partition(klab, lgirs, multiples, true, kidx, global_iridxs)
end