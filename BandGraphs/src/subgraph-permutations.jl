function set_pinned_subgraphs!(subgraphs_ps)
    # pin additional subgraphs via filter-criteria
    foreach(multiple_irrep_filter!, subgraphs_ps)
    foreach(solitary_subduction_path_filter!, subgraphs_ps)
    foreach(fake_weyl_filter!, subgraphs_ps)
    # TODO: Think longer about whether it is okay to pin these ahead of the "free-choice"
    #       pins below.

    # first, we exploit our freedom to pin one sub-block for each nonmax manifold; for each
    # nonmax manifold, we identify all sub-blocks connected to it and then pin whichever
    # of these subblocks has the most permutations
    seen_kidxs_nonmax = Set{Int}()
    for subgraph_ps in subgraphs_ps
        subgraph_ps.subgraph.monodromy_tie_via_Ω && continue # skip; unpermuted anyway
        kidx_nonmax = subgraph_ps.subgraph.p_nonmax.kidx
        kidx_nonmax ∈ seen_kidxs_nonmax && continue # already processed
        push!(seen_kidxs_nonmax, kidx_nonmax)

        running_maxperm = 0
        pidx = 0 # pinning index (into `subgraphs_ps`)
        for (idx, subgraph_ps′) in enumerate(subgraphs_ps)
            subgraph′ = subgraph_ps′.subgraph
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
    # similarly if there is only a single type of nonmaximal k-manifold
    if length(subgraph.p_max.multiples) == 1 || length(subgraph.p_nonmax.multiples) == 1
        subgraph_ps.permutations = nothing
        subgraph_ps.cols_or_rows = UNPERMUTED
        subgraph_ps.pinned = true
    end
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
end

function fake_weyl_filter!(subgraph_ps)
    # TODO

function subgraph_permutations(
        subgraph :: SubGraph{D};
        max_subblock_permutations :: Union{Nothing, <:Real} = nothing
        ) where D
    # if the subgraph is a Ω-trivial connection between two monodromy-related manifolds, we
    # deliberately do not want any permutations, since we are using these connections to
    # enforce energy-connectedness between same-irrep nodes across monodromy-pairs; detect
    # and return unpermuted subgraph in this case
    if subgraph.monodromy_tie_via_Ω
        return SubGraphPermutations{D}(subgraph, nothing, UNPERMUTED, #=pinned=# false)
    end

    # for each set of same-irrep nodes, we can perform all possible permutations of the
    # columns in the adjacency matrix associated with those nodes; if there are multiple
    # same-irrep node sets, we must perform all permutations across all sets
    multiples = subgraph.p_nonmax.multiples
    pss = permutations.(multiples)

    # in some awful cases, even `pss` can be too long for us to ever iterate through below,
    # so we have the option to bail out if # of permutations > `max_subblock_permutations`
    if !isnothing(max_subblock_permutations) # don't check if `nothing`
        Np = sum(length, pss)
        Np > max_subblock_permutations && return nothing
    end

    # if the permuted columns of the subgraph adjacency matrix are identical, however, there
    # is no point in including them, since the permutation would leave the associated
    # subgraph adjacency matrix unchanged; we filter them out below
    pss′ = [Vector{Int}[] for same_ir_idx in 1:length(pss)]
    for (same_ir_idx, ps) in enumerate(pss)
        As′_part = Set{Matrix{Int}}()
        ps′ = pss′[same_ir_idx]
        for p in ps
            A′_part = subgraph.A[:,p]
            A′_part ∈ As′_part && continue
            push!(As′_part, A′_part)
            push!(ps′, p)
        end
    end

    # finally, collect the distinct subgraph adjacency matrices for these permutations
    col_permutations = VectorProductIterator(pss′)
    subgraph_ps = SubGraphPermutations{D}(subgraph, col_permutations, COLS, #=pinned=# false)
    return subgraph_ps
end

function permute_subgraphs(
            subgraphs :: AbstractVector{SubGraph{D}};
            set_pinned :: Bool = true,
            max_subblock_permutations :: Union{Nothing, <:Real} = nothing
        ) where D

    # generate all subgraph permutations of each subgraph
    subgraphs_ps = Vector{SubGraphPermutations{D}}(undef, length(subgraphs))
    for (i, subgraph) in enumerate(subgraphs)
        subgraph_ps = subgraph_permutations(subgraph; max_subblock_permutations)
        if !isnothing(max_subblock_permutations) && isnothing(subgraph_ps)
            return nothing # bail out; subblock permutations > `max_subblock_permutations`
        end
        subgraph_ps::SubGraphPermutations{D}
        subgraphs_ps[i] = subgraph_ps
    end

    # set pinned subgraphs: for pinned subgraphs, we fix the subgraph and do not permute
    set_pinned && set_pinned_subgraphs!(subgraphs_ps)

    return subgraphs_ps
end

# TODO: some entries in `SUBDUCTIONSD[n]` contain trivial connections via the general point
#       Ω: we should probably just remove any nodes associated with Ω. In fact, there are
#       other examples of this where Bilbao/ISOTROPY include "special" k-points which
#       nevertheless have a trivial little group (presumably because it has the composition
#       of inversion and time-reversal?)