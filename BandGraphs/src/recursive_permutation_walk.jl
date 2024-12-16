function findall_separable_vertices_treecutting(
    criterion,
    n :: SymmetryVector{D},
    subts :: AbstractVector{SubductionTable{D}},
    lgirsd :: AbstractDict{String, Collection{LGIrrep{D}}};
    separable_degree :: Union{Nothing, <:Integer} = nothing,
    kws...
) where D
    T = Tuple{SeparabilityState, Tuple{LGIrrep{D}, Int}, Union{Nothing, NTuple{2, BandGraph{D}}}}
    bandg = build_subgraphs(n, subts, lgirsd)

    # find which vertices fulfil criterion; i.e., which vertices to check seperability for
    vs = findall_vertices(lgir -> isspecial(lgir) && (criterion(lgir) :: Bool), bandg)
    results = Vector{T}(undef, length(vs))
    isempty(vs) && return results

    # fast-path "necessary conditions" check which does not require permutation construction
    unfeasible_vs_idxs = fast_path_separability_irdim_infeasibility(bandg, vs, separable_degree)

    # full check, exploring permutations via recursive tree-cutting
    for (i,v) in enumerate(vs)
        if i in unfeasible_vs_idxs
            results[i] = (SeparabilityState(INSEPARABLE_CANNOT_GROUP), v, nothing)
            continue
        else
            r = findall_separable_vertices_treecutting(bandg, v; separable_degree, kws...)
            results[i] = (r[1], v, r[2])
        end
    end

    return results
end

function findall_separable_vertices_treecutting(
    bandg :: BandGraph{D},
    v :: Tuple{LGIrrep{D}, Int};
    with_weyl_filter = false,
    max_subblock_permutations = 1e5,
    kws...
) where D
    bandgp = BandGraphPermutations(bandg; with_weyl_filter, max_subblock_permutations)
    if isnothing(bandgp) # more subblock permutations than `max_subblock_permutations`
        s = SeparabilityState(BandGraphs.INFEASIBLE_TOO_MANY_SUBBLOCK_PERMUTATIONS)
        return s, nothing
    end
    bandgp :: BandGraphPermutations{D}
    return findall_separable_vertices_treecutting(bandgp, v; kws...)
end

function findall_separable_vertices_treecutting(
    bandgp :: BandGraphPermutations{D},
    v :: Tuple{LGIrrep{D}, Int};
    verbose :: Bool = false,
    kws...
) where D
    lgir = first(v)
    v_klab = klabel(lgir)
    subgraphs_ps = bandgp.subgraphs_ps
    
    # pick a first set of subgraphs that are connected to `v_klab` and have "trivial"
    # permutations only (i.e., length 1)
    idxs, klabs = add_unpermuted_max2nonmax2max_connections!(Int[], [v_klab], subgraphs_ps)
    
    # start the recursive search
    iter_count, perm_skips, sep, split_info = recursive_separability_search(bandgp, v, idxs, klabs; verbose, kws...)
    if verbose && !is_feasible(sep)
        println(sep, ": ", perm_skips, " vs. ", length(bandgp), " (", 
                    round(perm_skips/length(bandgp), digits=2), "%)")
    end
    return sep, split_info
end

function recursive_separability_search(
    bandgp  :: BandGraphPermutations{D},
    v       :: Tuple{LGIrrep{D}, Int}, 
    idxs    :: Vector{Int}, 
    klabs   :: Vector{String},
    recur_depth :: Int = 0, 
    iter_count  :: Int = 0, 
    perm_skips  :: BigInt = Base.GMP.ZERO;
    separable_degree :: Union{Nothing, Int} = nothing,
    max_permutations :: Real = 1e4,
    verbose :: Bool = true
) where D
    subgraphs_ps = bandgp.subgraphs_ps

    # check if we have already grown up to the "end" - of so, we can give final answers:
    if idxs == eachindex(subgraphs_ps)
        for bandg in bandgp
            iter_count += 1
            if iter_count > max_permutations
                sep = SeparabilityState(BandGraphs.INFEASIBLE_TOO_MANY_BANDGRAPH_PERMUTATIONS)
                return iter_count, perm_skips, sep, nothing
            end
            may_be_articulation = only(may_be_articulation_vertices(bandg, [v]))
            may_be_articulation || continue
            split_bandg = is_separable_at_vertex_check_all_splits(bandg, v;
                                                                  verbose, separable_degree)
            if !isnothing(split_bandg) # success: found a separable solution!
                sep = SeparabilityState(BandGraphs.SEPARABLE)
                return iter_count, perm_skips, sep, (bandg, split_bandg)
            end
        end
        # `bandgp` does not contain any separable solutions at `v`
        sep = SeparabilityState(BandGraphs.INSEPARABLE_CANNOT_SPLIT)
        return iter_count, perm_skips, sep, nothing
    end

    # include any new no-permutation max subgraphs connected to the above nonmax klabs (but
    # is not already in `idxs`)
    Nidxs = length(idxs)
    idxs, klabs = add_permuted_max2nonmax2max_connection!(idxs, klabs, subgraphs_ps)
    length(idxs) == Nidxs && error("failed to grow induced graph; infinite loop imminent")

    # we have now added something to `idxs` & `klabs`: there could now be more unpermuted
    # subgraphs connected to `klabs` than before: include them if they exist
    idxs, klabs = add_unpermuted_max2nonmax2max_connections!(idxs, klabs, subgraphs_ps)

    # build the induced (permuted) band graph corresponding to `bandgp.subgraphs_ps[idxs]`
    induced_bandgp = BandGraphs.subgraph_subset_induced_subgraph(bandgp, idxs)

    # how many permutations remain in the left-out bits of `bandg`
    remaining_permutations = count_remaining_permutations(idxs, subgraphs_ps)
    
    # iterate over permutations of the induced graph
    for (idx, induced_bandg) in enumerate(induced_bandgp)
        iter_count += 1
        if iter_count > max_permutations
            sep = SeparabilityState(BandGraphs.INFEASIBLE_TOO_MANY_BANDGRAPH_PERMUTATIONS)
            return iter_count, perm_skips, sep, nothing
        end
        may_be_articulation = only(may_be_articulation_vertices(induced_bandg, [v]))
        if !may_be_articulation
            perm_skips += remaining_permutations
            continue
        end
        
        split_bandg = is_separable_at_vertex_check_all_splits(induced_bandg, v;
                                                              verbose, separable_degree)
        if !isnothing(split_bandg)
            # lock permutations involved so far, and proceed to grow the induced subgraph
            # recursively; if it is a dead end, we continue in the loop
            partially_locked_bandgp = _copy_aliased(bandgp)
            lock_permutations!(partially_locked_bandgp, induced_bandgp, idx)

            iter_count, perm_skips, sep, split_info = recursive_separability_search(
                    partially_locked_bandgp, v, copy(idxs), copy(klabs),
                    recur_depth+1, iter_count, perm_skips;
                    separable_degree, max_permutations, verbose)
            if !isnothing(split_info)
                # split_info = (bandg, split_bandg)
                return iter_count, perm_skips, sep, split_info
            end
        else
            perm_skips += remaining_permutations
        end
    end
    sep = SeparabilityState(BandGraphs.INSEPARABLE_CANNOT_SPLIT)
    return iter_count, perm_skips, sep, nothing
end

function add_unpermuted_max2nonmax2max_connections!(
    idxs::Vector{Int},     # mutated
    klabs::Vector{String}, # mutated
    subgraphs_ps::AbstractVector{<:BandGraphs.SubGraphPermutations}
)   
    # idxs of subgraps connected to `klabs`, but not yet in `idxs`
    idxs_from_klabs = Int[]
    for (i, s_ps) in enumerate(subgraphs_ps)
        s_ps.subgraph.p_max.klab ∈ klabs || continue
        insorted(i, idxs) && continue
        if isnothing(s_ps.permutations) || isone(length(s_ps.permutations)) ||
            s_ps.is_chained
            push!(idxs_from_klabs, i)
        end
    end
    isempty(idxs_from_klabs) && return idxs, klabs

    # idxs of subgraphs that have max-to-nonmax-to-max connections that are pairwise trivial
    new_idxs = Int[]
    for i in idxs_from_klabs
        nonmax_klab = subgraphs_ps[i].subgraph.p_nonmax.klab
        for (j, s_ps) in enumerate(subgraphs_ps)
            insorted(j, idxs_from_klabs) && continue
            insorted(j, idxs) && continue
            s_ps.subgraph.p_nonmax.klab == nonmax_klab || continue
            if isnothing(s_ps.permutations) || isone(length(s_ps.permutations))
                i ∉ new_idxs && push!(new_idxs, i)
                j ∉ new_idxs && push!(new_idxs, j)
            elseif s_ps.is_chained
                chains_others_i = subgraphs_ps[i].chains_others
                if !isnothing(chains_others_i) && j ∈ chains_others_i
                    for k in chains_others_i
                        k ∉ new_idxs && !insorted(k, idxs) && push!(new_idxs, k)
                    end
                end
            end
        end
    end
    isempty(new_idxs) && return idxs, klabs

    # make sure we also included every chained subgraph; including things that were chained
    # away "multiple" max points
    for i in new_idxs
        s_ps_i = subgraphs_ps[i]
        if !isnothing(s_ps_i.chains_others)
            for j in s_ps_i.chains_others
                j ∉ new_idxs && !insorted(j, idxs) && push!(new_idxs, j)
            end
        end
    end

    # add new k-labels to `klabs` and append & re-sort `new_idxs` to `idxs`
    for i in new_idxs
        max_klabᵢ = subgraphs_ps[i].subgraph.p_max.klab
        max_klabᵢ ∉ klabs && push!(klabs, max_klabᵢ)
        nonmax_klabᵢ = subgraphs_ps[i].subgraph.p_nonmax.klab
        nonmax_klabᵢ ∉ klabs && push!(klabs, nonmax_klabᵢ)
    end
    #println("idxs_at_start: ", idxs)
    sort!(append!(idxs, new_idxs))
    
    # at this point, we have included new "trivially permuted" subgraphs; which means we now
    # how more things to add to `klabs`; in turn, the new klabs could be connected to MORE
    # "trivially permuted" subgraphs; we recurse to make sure we get them all
    return add_unpermuted_max2nonmax2max_connections!(idxs, klabs, subgraphs_ps)
end

function add_permuted_max2nonmax2max_connection!(
    idxs::Vector{Int},     # mutated
    klabs::Vector{String}, # mutated
    subgraphs_ps::AbstractVector{<:BandGraphs.SubGraphPermutations}
)
    # TODO: Will need to avoid/deal with picking/scanning over `chained` subgraphs - cannot
    #       add them without adding their "chainer" also (but they should always be
    #       connected by a nonmax partition - so maybe okay to include both?)

    # find the subgraph with the most connections connected to `klabs`, which is not already
    # in `idxs` and is not chained
    minp = typemax(Int)
    idxA = 0
    for (i, s_ps) in enumerate(subgraphs_ps)
        i ∈ idxs && continue
        s = s_ps.subgraph
        (s.p_max.klab ∈ klabs || s.p_nonmax.klab ∈ klabs) || continue
        s_ps.is_chained && continue
        if minp > length(s_ps)
            idxA, minp = i, length(s_ps)
        end
    end

    # if we did not find any permuted connections attached to `klabs`, we will look for
    # connections not attached to `klabs`; the resulting graph will not be connected.
    if minp == typemax(Int) # no new connection identified
        for (i, s_ps) in enumerate(subgraphs_ps)
            i ∈ idxs && continue
            s_ps.is_chained && continue
            if minp > length(s_ps)
                idxA, minp = i, length(s_ps)
            end
        end
    end
    idxA == typemax(Int) && error("failed to find a connection to add to the graph; `idxs` may equal `eachindex(subgraphs_ps)`?")

    # two cases:
    #   1. if the connection is from a nonmax partition in `idxs` to a max partition (not
    #      already in `idxs`, by construction), it is enough to add that & move on
    #   2. if the connection is from a max partition in `idxs` to a nonmax partition (again,
    #      not already in `idxs`), we need to add another connection to this new nonmax
    #      partition from some max partition, to avoid "dangling" nonmax partitions
    
    # branch on case 1 vs. case 2
    s_ps_A = subgraphs_ps[idxA]
    if s_ps_A.subgraph.p_nonmax.klab ∈ klabs
        # case 1 (new klab is a max partition)
        push!(idxs, idxA) # NB: could do `searchsortedfirst` & insert to preserve sorting
        # if idxA chained anything, we add its chainees
        if !isnothing(s_ps_A.chains_others)
            for i in s_ps_A.chains_others
                push!(idxs, i)
                klabᵢ = subgraphs_ps[i].subgraph.p_max.klab
                klabᵢ ∉ klabs && push!(klabs, klabᵢ)
            end
        end
    else
        # case 2 (new klab is a nonmax partition)
        idxB = 0
        minp = typemax(Int)
        nonmax_klab = s_ps_A.subgraph.p_nonmax.klab
        for (i, s_ps) in enumerate(subgraphs_ps) # find least-permuted new `nonmax_klab`-connection
            i == idxA && continue
            s_ps.subgraph.p_nonmax.klab == nonmax_klab || continue
            i ∈ idxs && continue # impossible?
            if minp > length(s_ps)
                idxB, minp = i, length(s_ps)
            end
        end
        idxB == typemax(Int) && error(lazy"unexpectedly could not find an unexplored connection via $nonmax_klab")

        push!(idxs, idxA, idxB) # NB: could do `searchsortedfirst` & insert to preserve sorting
        s_ps_B = subgraphs_ps[idxB]
        klab_A = s_ps_A.subgraph.p_nonmax.klab
        klab_B = s_ps_B.subgraph.p_max.klab
        klab_A ∈ klabs || push!(klabs, klab_A)
        klab_B ∈ klabs || push!(klabs, klab_B)

        # if idxA or idxB chained anything, we add their chainees
        if !isnothing(s_ps_A.chains_others)
            for i in s_ps_A.chains_others
                i == idxB && continue # already discovered; don't take twice
                push!(idxs, i)
                klabᵢ = subgraphs_ps[i].subgraph.p_max.klab
                klabᵢ ∉ klabs && push!(klabs, klabᵢ)
            end
        end
        if !isnothing(s_ps_B.chains_others)
            for i in s_ps_B.chains_others
                i == idxA && continue # already discovered; don't take twice (impossible?)
                push!(idxs, i)
                klabᵢ = subgraphs_ps[i].subgraph.p_max.klab
                klabᵢ ∉ klabs && push!(klabs, klabᵢ)
            end
        end
    end

    # re-sort resulting indices and check assumptions
    sort!(idxs)
    @assert allunique(idxs)
    @assert allunique(klabs)

    return idxs, klabs
end

function lock_permutations!(
        bandgp :: BandGraphPermutations{D},
        induced_bandgp :: BandGraphPermutations{D},
        idx :: Integer
    ) where D

    N = length(induced_bandgp)
    @boundscheck 1 ≤ idx ≤ N || throw(BoundsError(induced_bandgp, idx))
    idx = idx - 1
    for s_ps in induced_bandgp.subgraphs_ps
        s_ps.is_chained && continue # will be treated by the chaining subgraph
        isnothing(s_ps.permutations) && continue # no change needed

        r1 = length(s_ps)
        idx′ = div(idx, r1)
        sub_idx_i = idx - r1*idx′ + 1 # `i`th subscript in the linear-to-subscript problem
        idx = idx′
        locked_permutation = s_ps.permutations[sub_idx_i]

        # now the relevant permutation is `locked_permutation`; but we want to lock it in
        # `bandgp`, not `induced_bandgp`, so we need to find the associated subgraph in
        # `bandgp` first - then change its permutations to be just `locked_permutation`
        s = s_ps.subgraph
        j = findfirst(bandgp.subgraphs_ps) do s_ps′
            s′ = s_ps′.subgraph
            s′.p_max.klab == s.p_max.klab && s′.p_nonmax.klab == s.p_nonmax.klab
        end
        @assert length(bandgp.subgraphs_ps[j].permutations) == length(s_ps.permutations)
        bandgp.subgraphs_ps[j].permutations = BandGraphs.VectorProductIterator([[locked_permutation]])
    end
    return bandgp
end

function _copy_aliased(bandgp :: BandGraphPermutations{D}) where D
    subgraphs_ps′ = map(bandgp.subgraphs_ps) do s_ps
        BandGraphs.SubGraphPermutations{D}(
            s_ps.subgraph,
            s_ps.permutations, # don't need to unalias/copy; will change reference not contents
            s_ps.cols_or_rows, s_ps.pinned, s_ps.is_chained, s_ps.chains_others)
    end
    return BandGraphPermutations(bandgp.partitions, subgraphs_ps′)
end
#=
function recur_search(bandgp, klabs)
    qp = partition_induced_subgraph(bandgp, klabs)

    # check if separable
    #   - if no:  the permutation so far (involving `klabs` partitions only) can be abandoned; it's no good.
    #   - if yes: we can continue, and should recur, adding a new klab to klabs.
    #       - the permutation so far is still acceptable and worth exploring ensuing permutations of
end
=#

function count_remaining_permutations(idxs, subgraphs_ps)
    # count the number of permutations that are not "locked" by `idxs`
    c = Base.GMP.ONE # :: BigInt(1), without allocation
    for (i, s_ps) in enumerate(subgraphs_ps)
        # since `c` is a `BigInt`, we skip trivial multiplications (i.e., by 1) for perf.
        s_ps.is_chained || insorted(i, idxs) && continue
        c *= length(s_ps)
    end
    return c
end