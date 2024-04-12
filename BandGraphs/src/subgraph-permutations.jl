function set_pinned_subgraphs!(subgraphs)
    # we assume and exploit a sorting of subgraphs here, where entries are sorted
    # lexicographically in the sense that
    #       [(g.p_max.kidx, g.p_nonmax.kidx) for g in subgraphs]
    # increases first in its second index (nonmax manifold), then its first (max manifold).

    seen_kidxs_nonmax = Set{Int}()
    for subgraph in subgraphs
        kidx_nonmax = subgraph.p_nonmax.kidx
        if kidx_nonmax ∈ seen_kidxs_nonmax
            # a subgraph has already been pinned with this (nonmax) kidx, so we cannot pin
            # another block in the column associated with this kidx
            subgraph.pinned = false
            continue
        else
            push!(seen_kidxs_nonmax, kidx_nonmax)
            subgraph.pinned = true
        end
    end

    # pin additional subgraphs via filter-criteria
    multiple_irrep_filter!.(subgraphs)
    solitary_subduction_path_filter!.(subgraphs)
    fake_weyl_filter!.(subgraphs)

    return subgraphs
end

function multiple_irrep_filter!(subgraph)
    # if there's only a single type of irrep included in the considered maximal k-manifold,
    # there is no point in considering variations from it as they are indistinguishable;
    # similarly if there is only a single type of nonmaximal k-manifold
    if length(subgraph.p_max.multiples) == 1 || length(subgraph.p_nonmax.multiples) == 1
        subgraph.pinned = true
    end
end
function solitary_subduction_path_filter!(subgraph)
    # if there is only ever one irrep to subduce into, then there will be no permutations
    # and we might as well pin it for clarity; this is the case when 
    # `m = subgraph.p_nonmax_multiples` consists of solely single-element ranges s.t.
    # `all(perm -> length(perm) == 1, permutations.(m))`; an example is e.g., a subgraph
    # `s` = [S₁S₂ ↓ Q₂Q₂, S₃S₄ ↓ Q₁Q₁] which has no column-permutations for `s.A`;
    # this filter does not actually reduce the number of overall permutations, but better
    # shows that some subgraphs are effectively "locked" anyway (i.e., have no permutations)
    if all(m -> length(m) == 1, subgraph.p_nonmax.multiples)
        subgraph.pinned = true
    end
end

function fake_weyl_filter!(subgraph)
    # TODO
end

function subgraph_permutations(subgraph)
    # if the subgraph is pinned, we do not generate any permutations
    if subgraph.pinned
        As = [subgraph.A]
        return SubGraphPermutations(subgraph.p_max, subgraph.p_nonmax, As, #=pinned=# true)
    end
    
    # for each set of same-irrep nodes, we can perform all possible permutations of the
    # columns in the adjacency matrix associated with those nodes; if there are multiple
    # same-irrep node sets, we must perform all permutations across all sets
    multiples = subgraph.p_nonmax.multiples   
    pss = collect.(permutations.(multiples))
    Np = prod(length, pss) # aggregate number of same-irrep permutations, across irreps
    As = Vector{Matrix{Int}}(undef, Np)
    for (j,is) in enumerate(CartesianIndices(Tuple(Base.OneTo.(length.(pss)))))
        # TODO: do this type-stably and without `collect`ing `pss`
        cols = reduce(vcat, [pss[k][i] for (k,i) in enumerate(Tuple(is))])
        As[j] = subgraph.A[:,cols]
    end
    return SubGraphPermutations(subgraph.p_max, subgraph.p_nonmax, As, #=pinned=# false)
end

function permute_subgraphs(subgraphs; set_pinned::Bool=true)
    # find pinned subgraphs: only include a single, canonical, trivial permutation then
    set_pinned && set_pinned_subgraphs!(subgraphs)

    # generate all subgraph permutations of non-pinned subgraphs
    return subgraph_permutations.(subgraphs) # :: Vector{SubGraphPermutations}
end

# TODO: some entries in `SUBDUCTIONSD[n]` contain trivial connections via the general point
#       Ω: we should probably just remove any nodes associated with Ω. In fact, there are
#       other examples of this where Bilbao/ISOTROPY include "special" k-points which
#       nevertheless have a trivial little group (presumably because it has the composition
#       of inversion and time-reversal?)