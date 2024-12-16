function set_pinned_subgraphs!(subgraphs_ps; with_weyl_filter::Bool = true)
    # pin additional subgraphs via filter-criteria
    foreach(multiple_irrep_filter!, subgraphs_ps)
    foreach(solitary_subduction_path_filter!, subgraphs_ps)
    with_weyl_filter && fake_weyl_filter!(subgraphs_ps)
    # TODO: Think longer about whether it is okay to pin these ahead of the "free-choice"
    #       pins below.

    # first, we exploit our freedom to pin one sub-block for each nonmax manifold; for each
    # nonmax manifold, we identify all sub-blocks connected to it and then pin whichever
    # of these subblocks has the most permutations
    seen_kidxs_nonmax = Set{Int}()
    for subgraph_ps in subgraphs_ps
        subgraph_ps.subgraph.monodromy_tie_via_Ω && continue # skip; unpermuted anyway
        subgraph_ps.is_chained && continue # skip; tied to something else anyway
        kidx_nonmax = subgraph_ps.subgraph.p_nonmax.kidx
        kidx_nonmax ∈ seen_kidxs_nonmax && continue # already processed
        push!(seen_kidxs_nonmax, kidx_nonmax)

        running_maxperm = 0
        pidx = 0 # pinning index (into `subgraphs_ps`)
        for (idx, subgraph_ps′) in enumerate(subgraphs_ps)
            subgraph′ = subgraph_ps′.subgraph
            subgraph_ps′.is_chained && continue # skip; tied to something else anyway
            # only consider subgraphs connected to `kidx_nonmax`
            kidx_nonmax == subgraph′.p_nonmax.kidx || continue
            permutations′ = subgraph_ps′.permutations
            nperms = isnothing(permutations′) ? 1 : length(permutations′)
            if nperms > running_maxperm
                running_maxperm = nperms
                pidx = idx
            end
        end
        pidx == 0 && error("unexpectedly found no subgraphs connected to kidx")

        # remove permutations from `subgraphs_ps[pidx]` & pin it
        subgraphs_ps[pidx].permutations = nothing
        subgraphs_ps[pidx].cols_or_rows = UNPERMUTED
        subgraphs_ps[pidx].pinned = true
    end
    
    return subgraphs_ps
end

function multiple_irrep_filter!(subgraph_ps)
    subgraph = subgraph_ps.subgraph
    # if there's only a single type of irrep included in the considered maximal k-manifold,
    # there is no point in considering variations from it as they are indistinguishable;
    # the same idea IS NOT VALID if there is only a single type of nonmaximal k-manifold
    if length(subgraph.p_max.multiples) == 1
        subgraph_ps.permutations = nothing
        subgraph_ps.cols_or_rows = UNPERMUTED
        subgraph_ps.pinned = true
    end
    return subgraph_ps
end
function solitary_subduction_path_filter!(subgraph_ps)
    # if there is only ever one irrep to subduce into, then there will be no permutations
    # and we might as well pin it for clarity; this is the case when 
    # `m = subgraph.p_nonmax.multiples` consists of solely single-element ranges s.t.
    # `all(perm -> length(perm) == 1, permutations.(m))`; an example is e.g., a subgraph
    # `s` = [S₁S₂ ↓ Q₂Q₂, S₃S₄ ↓ Q₁Q₁] which has no column-permutations for `s.A`;
    # this filter does not actually reduce the number of overall permutations, but better
    # shows that some subgraphs are effectively "locked" anyway (i.e., have no permutations)
    if all(m -> length(m) == 1, subgraph_ps.subgraph.p_nonmax.multiples)
        subgraph_ps.permutations = nothing
        subgraph_ps.cols_or_rows = UNPERMUTED
        subgraph_ps.pinned = true
    end
    return subgraph_ps
end

function fake_weyl_filter!(subgraphs_ps)
    # TODO
    seen_subgraph_ps = zeros(Bool, length(subgraphs_ps))
    chained_max_kidxs = Set{Int}()
    for (i,subgraph_ps) in enumerate(subgraphs_ps)
        seen_subgraph_ps[i] == true && continue # already processed
        if subgraph_ps.subgraph.monodromy_tie_via_Ω
            seen_subgraph_ps[i] = true
            continue # skip; unpermuted anyway
        end
        kidx_nonmax = subgraph_ps.subgraph.p_nonmax.kidx
        kidx_max    = subgraph_ps.subgraph.p_max.kidx
        #row_sorted_A = sortslices(subgraph_ps.subgraph.A, dims=1, rev=true)
        for (j, subgraph′_ps) in enumerate(subgraphs_ps)
            j > i || continue # we will have already processed these
            seen_subgraph_ps[j] == true && continue # also already processed
            kidx_nonmax′ = subgraph′_ps.subgraph.p_nonmax.kidx
            kidx_nonmax′ == kidx_nonmax || continue # not same nonmax; cannot compare
            kidx_max′    = subgraph′_ps.subgraph.p_max.kidx
            kidx_max′ ∈ chained_max_kidxs && continue # already chained; cannot chain again
            subgraph′_ps.pinned && continue # already pinned by other means; no point in chaining
            #row_sorted_A′ = sortslices(subgraph′_ps.subgraph.A, dims=1, rev=true)
            #if row_sorted_A == row_sorted_A′
            if is_weyl_block_related(subgraph_ps.subgraph, subgraph′_ps.subgraph)
                if (subgraph_ps.permutations != subgraph′_ps.permutations &&
                    (!isnothing(subgraph_ps.permutations) && 
                     !isnothing(subgraph′_ps.permutations)))
                    error("unexpectedly found two subgraphs with same A but different, non-nothing permutations")
                end
                # the Weyl filter applies! the two subgraphs must permute in tandem; all
                # independent permutations lead to fake Weyl crossings, which gap out
                if isnothing(subgraph′_ps.permutations)
                    # then prefer to chain `i` to `j`; we'd rather chain to something that
                    # doesn't have any permutations (`j`) than something which does (`i`)
                    push!(@something(subgraph′_ps.chains_others,
                                     subgraph′_ps.chains_others=Int[]), i)
                    subgraph_ps.is_chained = true
                    seen_subgraph_ps[i] = true
                    break # cannot chain anything else to `i` now; already chained
                else
                    # otherwise, just chain `j` to `i`
                    push!(@something(subgraph_ps.chains_others,
                                     subgraph_ps.chains_others=Int[]), j)
                    subgraph′_ps.is_chained = true # chain `j` to `i`
                    seen_subgraph_ps[j] = true
                end
                push!(chained_max_kidxs, kidx_max′)
                # TODO: Should not just blindly chain to first possible option - should look
                #       which chaining options chain the subgraphs with most permutations
            end
        end
    end
    return subgraphs_ps
end

function is_weyl_block_related(subgraph::SubGraph{D}, subgraph′::SubGraph{D}) where D
    # for each nonmax irrep set, extract the associated nonvanishing block (i.e., trim zero
    # rows in each nonmax-same-irrep-column set); if these such blocks for `subgraph` and
    # `subgraph′` are identical for each nonmax-same-irrep subset, any permutations would
    # result in fake-weyl-equivalent connections
    s, s′ = subgraph, subgraph′
    p_nonmax, p_nonmax′ = s.p_nonmax, s′.p_nonmax
    multiples, multiples′ = p_nonmax.multiples, p_nonmax′.multiples
    if !(p_nonmax.klab == p_nonmax′.klab && p_nonmax.iridxs == p_nonmax′.iridxs && 
         multiples == multiples′)
        error("unexpectedly was passed subgraphs with unrelated nonmax irrep partitions")
    end
    A, A′ = s.A, s′.A
    for m in multiples
        Aₘ, Aₘ′ = (@view A[:, m]), (@view A′[:, m])
        rowidxs = findall(!iszero, eachrow(Aₘ))
        rowidxs′ = findall(!iszero, eachrow(Aₘ′))
        nonzero_block_Aₘ = @view Aₘ[rowidxs, :]
        nonzero_block_Aₘ′ = @view Aₘ′[rowidxs′, :]
        if nonzero_block_Aₘ != nonzero_block_Aₘ′
            return false
        end
    end
    return true
end

function subgraph_permutations(
        subgraph :: SubGraph{D};
        max_subblock_permutations :: Union{Nothing, <:Real} = 1e5
        ) where D
    # if the subgraph is a Ω-trivial connection between two monodromy-related manifolds, we
    # deliberately do not want any permutations, since we are using these connections to
    # enforce energy-connectedness between same-irrep nodes across monodromy-pairs; detect
    # and return unpermuted subgraph in this case
    if subgraph.monodromy_tie_via_Ω
        return SubGraphPermutations{D}(subgraph, nothing, UNPERMUTED, #=pinned=# false)
    end

    # for each set of same-irrep nodes, with indices `ms` from `multiples`, we can perform
    # all possible permutations of the columns in the adjacency matrix associated with those
    # nodes; if more than one element in `multiples`, we must perform all product 
    # permutations across each `ms = multiples[i]`.
    # if within each such subset, i.e., each `ms=multiples[i]`, there are identical columns
    # of `subgraph.A`, a permutation would not change the adjacency matrix - i.e., can be 
    # omitted (any two permutations `p1` & `p2` that permute the identical columns of `A`
    # have `A[:,p1] == A[:,p2]`).
    # rather than test which permutations do this by brute force, which is slow, we think
    # about it as a multiset permutation problem: what we need are not all permutations of
    # `ms`, but the multiset permutations accounting for equal column redundancies. The
    # equal-column situation happens e.g., for situations like Γ₁ ↓ Λ₁ + Λ₁ (i.e., subducing
    # into the same nonmax irrep multiple times). Below, we determine the associated
    # equal-column multiset for each `ms`, exploiting the "diagonal" structure of `A` to
    # identify equal columns speedily (they must be consecutive, cf. construction of 
    # the canonical subgraphs `A`)
    multiples = subgraph.p_nonmax.multiples
    A = subgraph.A
    pss = [Vector{Int}[] for _ in 1:length(multiples)]
    for (same_ir_idx, ms) in enumerate(multiples)
        # find frequency of `A`'s columns, exploiting that equal columns must be consecutive
        Aₘ = @view A[:, ms]
        n = 0
        row = firstindex(A, 1) - 1
        f = Vector{Int}()
        for aᵢ in eachcol(Aₘ)
            if row < firstindex(A, 1) || iszero(aᵢ[row]) # exploit `A`'s structure
                row = something(findfirst(!iszero, aᵢ))
                n += 1
                push!(f, 1)
            else #!iszero(aᵢ[row])
                f[n] += 1
            end
        end
        ps′ = multiset_permutation_indices(f) # local indices into `ms`
        
        # in some awful cases, even `ps′` can have too many elements to ever iterate
        # through below; if we exceed `max_subblock_permutations`, bail out
        !isnothing(max_subblock_permutations) && length(ps′)>max_subblock_permutations && return nothing

        offset = first(ms) - 1 # `ms` is consecutive, but may not start at 1
        ps = [(p .= p .+ offset) for p in ps′]
        pss[same_ir_idx] = ps
    end

    # finally, collect the distinct subgraph adjacency matrices for these permutations
    col_permutations = VectorProductIterator(pss)
    if length(col_permutations) ≠ 1
        return SubGraphPermutations{D}(subgraph, col_permutations, COLS, #=pinned=# false)
    else
        return SubGraphPermutations{D}(subgraph, nothing, UNPERMUTED, #=pinned=# false)
    end
end

function permute_subgraphs(
            subgraphs :: AbstractVector{SubGraph{D}};
            set_pinned :: Bool = true,
            max_subblock_permutations :: Union{Nothing, <:Real} = nothing,
            with_weyl_filter :: Bool = true
        ) where D

    # generate all subgraph permutations of each subgraph
    subgraphs_ps = Vector{SubGraphPermutations{D}}(undef, length(subgraphs))
    for (i, subgraph) in enumerate(subgraphs)
        subgraph_ps = subgraph_permutations(subgraph; max_subblock_permutations)
        if !isnothing(max_subblock_permutations) && isnothing(subgraph_ps)
            return nothing # bail out; sub-block permutations > `max_subblock_permutations`
        end
        subgraphs_ps[i] = subgraph_ps::SubGraphPermutations{D}
    end

    # set pinned subgraphs: for pinned subgraphs, we fix the subgraph and do not permute
    set_pinned && set_pinned_subgraphs!(subgraphs_ps; with_weyl_filter)

    return subgraphs_ps
end

# TODO: some entries in `SUBDUCTIONSD[n]` contain trivial connections via the general point
#       Ω: we should probably just remove any nodes associated with Ω. In fact, there are
#       other examples of this where Bilbao/ISOTROPY include "special" k-points which
#       nevertheless have a trivial little group (presumably because it has the composition
#       of inversion and time-reversal?)