"""
    complete_split(bandg::BandGraph, v::Tuple{String, Int})

Compute an iterator over the allowable complete splits of a vertex `v` in a band
graph `bandg`.

`v` must be a vertex in the maximal partition of the band graph; this vertex is split into
`irdim(lgir)` one-dimensional "surrogate" irreps, with `lgir` denoting the irrep associated
with `v`. The split is "complete" in the sense that it splits the irrep fully into
one-dimensional surrogate split-off irreps.
It must be specified in terms of its irrep label and its index in the multiples of that 
irrep. E.g., if `bandg` contains Γ₁ twice, and we want to split the first of such irrep
contained in `bandg`, we have `v = ("Γ₁", 1)`.
"""
function complete_split(
    bandg :: BandGraph{D}, 
    v::Tuple{String, Int}; 
    separable_degree :: Union{Nothing, Int} = nothing
    ) where D
    
    partitions = bandg.partitions
    irlab, irmul = v

    # first, we find out which partition `v` belongs to - and where in that partition, i.e.,
    # which `iridxs` index corresponds to `v`
    tmp = find_vertex_in_partitions(partitions, irlab, irmul)
    isnothing(tmp) && error(lazy"failed to find vertex `v = `$v in `bandg` partitions")
    kidx, iridx, midx = tmp
    p = partitions[kidx]
    p.maximal || error(lazy"input vertex $(v[1]) is not associated with a maximal k-point")
    irdim_v = irdim(p.lgirs[iridx])
    irdim_v == 1 && error("cannot split a one-dimensional irrep")
    isnothing(separable_degree) && (separable_degree = irdim_v)

    # now that we have some information regarding the vertex `v`, we proceed to split the
    # associated irrep & find all the possible admissible splits of its subgraphs
    splits_bandg = _complete_split(bandg, irmul, kidx, iridx, midx, irdim_v,
                                   separable_degree)
    if splits_bandg === nothing
        # unable to split; likely connected to a nonmaximal irrep that is not 1-dimensional
        return nothing
    end

    # check if vertex has a monodromy-related partner; if so, it must be split with `v`
    klab = klabel(irlab)
    irlab′ = replace(irlab, klab=>klab * "′")
    irmul′ = irmul
    tmp′ = find_vertex_in_partitions(partitions, irlab′, irmul′)
    if !isnothing(tmp′)
        kidx′, iridx′, midx′ = tmp′
        # we must also split up the monodromy-related vertex - and in general, we can
        # just explore all subgraphs for the affected subgraphs as well
        splits_bandg = _complete_split(splits_bandg, irmul′, kidx′, iridx′, midx′, irdim_v,
                                       separable_degree, #=monodromy=# true)
    end

    return splits_bandg
end
function complete_split(
    bandg :: BandGraph{D},
    v::Tuple{LGIrrep{D}, Int};
    kws...
    ) where D
    lgir, irmul = v
    return complete_split(bandg, (label(lgir), irmul); kws...)
end

function _complete_split(
    bandg :: Union{BandGraph{D}, BandGraphPermutations{D}},
    irmul :: Int,
    kidx :: Int,
    iridx :: Int,
    midx :: Int,
    irdim :: Int,
    separable_degree :: Int = irdim,
    monodromy :: Bool = false,
    ) where D

    partitions = bandg.partitions

    # modify the specific partition that `v = (kidx, iridx, midx)` belongs to by removing
    # it and adding two new "split-off" surrogate vertices in its place
    split_p = build_split_partition(partitions, kidx, iridx, midx, irdim)

    # all subsequent partitions in `bandg` have their `iridxs` incremented by `irdim-1`
    # since the degree of the graph is now `irdim-1` larger
    split_partitions = Vector{Partition{D}}(undef, length(partitions))
    for (kidx′, p′) in enumerate(partitions)
        if kidx′ < kidx
            split_partitions[kidx′] = p′
        elseif kidx′ == kidx
            split_partitions[kidx′] = split_p
        else # kidx′ > kidx
            split_partitions[kidx′] = Partition(p′.klab, p′.lgirs, p′.multiples, p′.maximal, 
                                                p′.kidx, p′.iridxs .+ (irdim-1))
        end
    end

    # `split_partitions` is now the new "store" of vertex data for the graph; but the actual
    # graph structure is specified in `subgraphs`, so this must be updated accordingly;
    # this is the real work of the splitting - and where we must create the different kinds
    # of admissible vertex splits by reconnecting the edges in a meaningful and
    # comprehensive way
    subgraphs = if bandg isa BandGraphPermutations
        getfield.(bandg.subgraphs_ps, Ref(:subgraph))
    else
        bandg.subgraphs
    end :: Vector{SubGraph{D}}

    # first, find the parts of `subgraphs` that are affected by the split
    affected_subgraphs_idxs = findall(s -> s.p_max.kidx == kidx, subgraphs)
    affected_subgraphs = subgraphs[affected_subgraphs_idxs]

    # check that we are not trying to split a subgraph that already has permutations; this
    # is annoying to implement and not currently done (TODO)
    if bandg isa BandGraphPermutations
        for idx in affected_subgraphs_idxs
            subgraph_ps = bandg.subgraphs_ps[idx]
            if !isnothing(subgraph_ps.permutations)
                error("cannot split a subgraph which already has multiple permutations \
                      (not currently implemented)")
            end
        end
    end            

    # now, get started on finding the permutations of the splitting, in terms of a new 
    # subgraph matrix
    row_in_A = partitions[kidx].multiples[midx][irmul]
    affected_subgraphs_new_split_A = Vector{Matrix{Int}}(undef, length(affected_subgraphs))
    affected_subgraphs_row_permutations = Union{Vector{Vector{Int}}, Nothing}[
                                        Vector{Int}[] for _ in 1:length(affected_subgraphs)]
    for (i, s) in enumerate(affected_subgraphs)
        A = s.A

        cols_in_A = findall(!iszero, A[row_in_A,:])
        A_part = A[row_in_A, cols_in_A]
        # At this point, `A_part` will usually be a vector of `1`s, indicating that the
        # maximal irrep at `v` (corresponding to `row_in_A`) is connected to nonmaximal
        # irreps that are all 1-dimensional. In this case, the split is relatively straight-
        # forward conceptually, as it involves the row-permutations of an identity matrix.
        # E.g., for a 2D maximal irrep connected to two 1D nonmaximal irreps, we would have
        # `A_part = [1,1]` and the splits would splan [1 0; 0 1] and [0 1; 1 0]. The row
        # permutations correspond to relabelling the "split-off" vertices.
        # However, more generally, the maximal irrep could be connected to nonmaximal irreps
        # that are not 1-dimensional. In this case, the splits are somewhat more
        # complicated - we explain the procedure below:
        # Assuming `A_part = [n₁, n₂, ... nₘ]` where `nᵢ` is the dimensionality of the
        # `i`-th nonmaximal irrep connected to the maximal irrep, we consider the unique
        # row-permutations of the following column-block-identity-like matrix:
        #    A = [1_n₁ 0_n₁ ... 0_n_1; 0_n₂ 1_n₂ ... 0_n₂; ...; 0_nₘ 0_nₘ ... 1_nₘ]
        # with 1_nᵢ and 0_nᵢ denoting `nᵢ`-dimensional unit and zero vectors. We want to
        # find all row-permutations of this matrix, corresponding to all possible
        # relabellings of the split-off vertices. Not all permutations of the rows will give
        # a new matrix, however: in particular, permutations within an `nᵢ`-block will not.
        # Generally, what we seek here are the indices of the so-called multiset
        # permutations of a multiset {n₁ ⋅ 1, n₂ ⋅ 2, ... nₘ ⋅ m}. We implement the tooling
        # to find these in `multiset_permutation_indices.jl`.
        # Technically, the new split-off subgraph adjacency matrix probably doesn't fulfil
        # all the sum-rules we have written up for the subblock adjacency matrices of band
        # graphs, but it doesn't matter as those that might be broken are not used anywhere.
        sum(A_part) == irdim || error("unexpected sum of nonmaximal irrep dimensions")
        # fast-path check
        if length(A_part) < separable_degree
            # If `separable_degree` is bigger than the number of nonmaximal irreps
            # connected to `v` in any given partition, there is no hope for a split to lead
            # to a separable configuration; in this case, we can skip further consideration
            # to speed things up
            return nothing
        end
        
        # now, we replace A_part by the row permutations of the above-mentioned matrix.
        new_A_part = column_block_identity_matrix(A_part)
        # if this is the first of the affected subgraphs, we only include the unpermuted
        # `new_A_part` matrix, since the inclusion of additional permutations would give
        # graphs that are isomorphic under edge contraction with respect to the split-off
        # surrogate vertices; if it's a monodromy-related vertex, we just permute in any
        # case to be safe
        top = @view A[1:row_in_A-1, :]
        bottom = @view A[row_in_A+1:end, :]
        new_A_part_full_col_size = zeros(Int, irdim, size(A, 2))
        new_A_part_full_col_size[:, cols_in_A] .= new_A_part
        affected_subgraphs_new_split_A[i] = vcat(top, bottom, new_A_part_full_col_size)
        for new_rows_in_A in multiset_permutation_indices(A_part)
            if first(s.p_nonmax.klab) == 'Ω'
                # we want specifically to retain a mapping relationship between monodromy-
                # paired irreps that are tied together via Ω<...> since these edges signify
                # that we physically require monodromy-paired irreps to have the same energy
                affected_subgraphs_row_permutations[i] = nothing
                break
            elseif !monodromy && i == 1
                # we pin this split-subgraph; i.e., it provides the reference point for
                # other subgraphs split permutations
                affected_subgraphs_row_permutations[i] = nothing
                break
            end
            new_rows_in_A .+= size(A,1) - 1 # go from 1-based to referencing the last row 
            # in `A`, which is where the new split vertices will be placed:
            # `new_rows_in_A` covers a contiguous set of the rows of the new split-subgraph
            # adjacency block, but not necessarily consecutively. I cannot wrap my head
            # around whether this could be expressed using a `VectorProductIterator` as we
            # do for `SubGraphPermutations` but it doesn't matter: there can be at most 
            # 720 permutations (`factorial(max(irdim))` = 6!) of this kind, so we can just
            # list them out explicitly
            row_permutation = collect(1:(size(A,1)+irdim-1))
            row_permutation[size(A,1):size(A,1)+irdim-1] .= new_rows_in_A
            push!(affected_subgraphs_row_permutations[i], row_permutation)
        end
    end

    # build a new `split_subgraphs_ps`, using `split_partitions`, `affected_subgraphs_idxs`,
    # `affected_subgraphs_new_split_A`, & `affected_subgraphs_row_permutations`
    i′ = 0
    split_subgraphs_ps = Vector{SubGraphPermutations{D}}(undef, length(subgraphs))
    for (i, s) in enumerate(subgraphs)
        new_split_A = if i ∉ affected_subgraphs_idxs
            s.A
        else
            affected_subgraphs_new_split_A[i′+=1]
        end
        p_max    = split_partitions[s.p_max.kidx]
        p_nonmax = split_partitions[s.p_nonmax.kidx]
        split_subgraph = SubGraph(p_max, p_nonmax, new_split_A)
        if i ∈ affected_subgraphs_idxs
            # we already checked previously that subgraphs[i] was not already permuted, so
            # we can go ahead and use the split-permutations without worry that we're
            # overwriting existing permutations
            _permutations = affected_subgraphs_row_permutations[i′]
            permutations, cols_or_rows, pinned = if isnothing(_permutations)
                nothing, UNPERMUTED, true
            else
                # `BandGraphPermutations` expects a `VectorProductIterator` for the
                # permutations; here, we actually only want to iterate over `_permutations`,
                # so we can just "wrap" `_permutations` in `[]` and then
                # `VectorProductIterator`, since `VectorProductIterator([a])` iterates like
                # `a::Vector` itself.
                VectorProductIterator([something(_permutations)]), ROWS, false
            end
            split_subgraphs_ps[i] = SubGraphPermutations(
                                        split_subgraph, permutations, cols_or_rows, pinned)
        else
            # subgraphs[i] might have represented an already-permuted subgraph; in this case
            # it would usually represent the graph having already been passed through
            # `_complete_split` once before); i.e., this is the `monodromy`-update and we're
            # here copying over the permutations from the previous split
            permutations, cols_or_rows, pinned = if bandg isa BandGraphPermutations
                subgraph_ps = bandg.subgraphs_ps[i]
                (subgraph_ps.permutations, subgraph_ps.cols_or_rows, subgraph_ps.pinned)
            else
                nothing, UNPERMUTED, false
            end

            split_subgraphs_ps[i] = SubGraphPermutations(
                                        split_subgraph, permutations, cols_or_rows, pinned)
        end
    end

    split_bandg = BandGraphPermutations(split_partitions, split_subgraphs_ps)
    return split_bandg
end

# given input `nv = [n₁, n₂, ... nₘ]`, return the column-block-identity-like matrix:
#    A = [1_n₁ 0_n₁ ... 0_n_1; 0_n₂ 1_n₂ ... 0_n₂; ...; 0_nₘ 0_nₘ ... 1_nₘ]
# with 1_nᵢ and 0_nᵢ denoting `nᵢ`-dimensional unit and zero vectors.
function column_block_identity_matrix(nv :: Vector{Int})
    cols = length(nv) # "m" in the above
    rows = sum(nv)
    A = zeros(Int, rows, cols)
    colidx = 1
    for (i, nᵢ) in enumerate(nv)
        nᵢ > 0 || error("unexpected zero or negative dimensionality")
        A[colidx:colidx + nᵢ - 1, i] .= 1
        colidx += nᵢ
    end
    return A
end

function find_vertex_in_partitions(
            partitions :: Vector{<:Partition}, irlab :: AbstractString, irmul :: Int)
    klab = klabel(irlab)
    for (kidx, p) in enumerate(partitions)
        klab == p.klab || continue
        kidx == p.kidx || error("unexpected misalignment of partition `kidx` and its index in the `bandg.partitions` vector")

        idx = 0 # below is equivalent to idx=findall(...)[irmul], but without allocations
        for _ in 1:irmul 
            idx = findnext(lgir -> label(lgir) == irlab, p.lgirs, idx + 1)
            if isnothing(idx)
                error(lazy"could not find irrep label $irlab with multiplicity $irmul at $klab ($p.lgirs)")
            end
        end
        iridx = idx :: Int
        midx = something(findfirst(∋(iridx), p.multiples))

        # return indices corresponding to `v = (irlab, irmul)` into:
        #   - kidx:  partitions[kidx]
        #   - iridx: partitions[kidx].lgirs[iridx]
        #   - midx:  partitions[kidx].multiples[midx]
        return (kidx, iridx, midx)
    end
    return nothing
end

function build_split_partition(
        partitions :: AbstractVector{Partition{D}}, kidx, iridx, midx, irdim) where D
          
    # we want to build a new, augmented partition, where the vertex corresponding to
    # `(kidx, iridx, midx)` is split-up into two new "split-surrogate" vertices; this
    # new augmented partition differs from the original partition in `lgirs`, `multiples`,
    # and `iridxs`

    # original partition
    p₀ = partitions[kidx]

    # augmented `multiples`
    multiples₀ = p₀.multiples
    complete_removal = length(multiples₀[midx]) == 1
    # ↑ if `true`, there will be no more irreps of the removed kind in the partition after
    # removal; if `false`, one or more irreps of the removed kind will remain in the
    # partition, after removal

    split_multiples = Vector{UnitRange{Int}}(undef, length(multiples₀) + !complete_removal)
    split_multiples[1:midx-1] = multiples₀[1:midx-1] # entries before `midx`
    for i in (midx+!complete_removal):length(split_multiples)-1 # entries after `midx`
        split_multiples[i] = multiples₀[i+complete_removal] .- 1
    end
    # insert augmented-irrep multiples
    split_multiples[end] = (last(multiples₀[end])) : (last(multiples₀[end]) + (irdim-1))
    if !complete_removal # retain some of removed multiples
        split_multiples[midx] = first(multiples₀[midx]) : (last(multiples₀[midx]) - 1)
    end
    
    # augmented irreps
    lgirs₀ = p₀.lgirs
    split_lgirs = [lgir for (i, lgir) in enumerate(lgirs₀) if i != iridx]

    # insert two dummy irreps corresponding to the split-off irreps; we give them "trivial"
    # irrep matrices, since all we care about is their labels - for the labels, we give
    # them the original label, appended with an "ˣ" suffix
    lgir₀ = lgirs₀[iridx]
    split_irlab = label(lgir₀) * 'ˣ'
    split_lgir = LGIrrep(
        split_irlab, 
        group(lgir₀),
        [ones(ComplexF64, 1, 1) for _ in 1:length(group(lgir₀))],
        [zeros(Float64, D) for _ in 1:length(group(lgir₀))],
        REAL, false)
    append!(split_lgirs, fill(split_lgir, irdim)) # insert `irdim` times

    # augmented indices into irreps
    iridxs₀ = p₀.iridxs
    split_iridxs = first(iridxs₀) : (last(iridxs₀) + (irdim-1))

    return Partition(p₀.klab, split_lgirs, split_multiples, p₀.maximal, p₀.kidx,
                     split_iridxs)
end


# a split can be invalid if it involves a monodromy-paired vertex: in particular, if the 
# split happens to disconnect the monodromy-paired vertices, then the split is invalid since
# these vertices must always be connected (i.e., live together in a connected component);
# equivalently, all monodromy-related vertices v and v′ must be connected by a path in the
# graph; this function checks that and is intended to be used to "filter out" invalid 
# splits that might be generated by the combinatorial splitting in `complete_split`
function is_valid_split(
        bandg :: BandGraph,
        g :: Graph = assemble_simple_graph(bandg),
        next_work::Vector{Int} = Vector{Int}(undef, nv(g)), # work buffer
        seen_work::Vector{Bool} = zeros(Bool, nv(g))        # work buffer
        )
    partitions = bandg.partitions
    any(p->last(p.klab) == '′', partitions) || return true

    for p′ in Iterators.reverse(partitions)
        klab′ = p′.klab
        last(klab′) == '′' || continue
        
        klab = SubString(klab′, firstindex(klab′), prevind(klab′, lastindex(klab′))) # faster `chop(klab′)` 
        i = findfirst(p -> p.klab == klab, partitions)
        if isnothing(i)
            error(lazy"could not find a monodromy-related partition for $klab′ (trying to match to $klab); available partitions are $([p.klab for p in partitions])")
        end
        p = partitions[i] # "original" partition (not monodromy-translated)

        iridxs′ = p′.iridxs
        iridxs  = p.iridxs

        # now check whether there is a path in `g` from each `iridx` to each `iridx′`
        for (iridx, iridx′) in zip(iridxs, iridxs′)
            # we use a modified version of Graphs.jl's `has_path` which allows us to reuse
            # the work arrays and hence is faster for repeated use
            if !has_path_with_workarrays!(g, iridx, iridx′, next_work, seen_work)
                return false
            end
        end
    end
    return true
end

# a variant of the `has_path` implementation from https://github.com/JuliaGraphs/Graphs.jl/pull/406,
# which additionally takes work buffers for `next` and `seen` arrays, allowing them to be
# reused on repeated calls. Checks if there is a path from vertex `u` to `v` in graph `g`.
function has_path_with_workarrays!(
    g :: AbstractGraph{T},
    u :: Integer,
    v :: Integer,
    next :: Vector{T},    # work buffer (mutated)
    seen :: Vector{Bool}; # work buffer (mutated)
    reset_seen :: Bool = true
    ) where T

    u == v && return true # cannot be separated
    reset_seen && fill!(seen, false)

    front = back = 1
    next[front] = u
    seen[u] = true
    while front <= back
        src = next[front]
        front += 1 # dequeue; pop from queue
        for vertex in outneighbors(g, src)
            vertex == v && return true
            if !seen[vertex]
                back += 1
                next[back] = vertex # enqueue; push onto queue
                seen[vertex] = true
            end
        end
    end
    return false
end


function nv_after_split(bandg, v)
    lgir = v[1]
    irdim_v = irdim(lgir)
    # after split, there will be `irdim_v - 1` extra vertices in the graph - except if the
    # vertex `v` has a monodromy partner, then there will be `2(irdim_v - 1)` extra vertices
    # check if vertex has a monodromy-related partner; if so, it must be split with `v`
    klab = klabel(lgir)
    irlab′ = klab * "′"
    has_monodromy_partner = any(bandg.partitions) do p
        p.maximal && p.klab == irlab′
    end
    extra_nv = has_monodromy_partner ? 2(irdim_v - 1) : irdim_v - 1
    return nv(bandg) + extra_nv
end