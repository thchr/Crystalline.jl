using Crystalline: irdim
using BandGraphs: SubGraphPermutations, BandGraphPermutations, assemble_simple_graph

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
function complete_split(bandg :: BandGraph{D}, v::Tuple{String, Int}) where D
    
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

    # now that we have some information regarding the vertex `v`, we proceed to split the
    # associated irrep & find all the possible admissible splits of its subgraphs
    splits_bandg = _complete_split(bandg, irmul, kidx, iridx, midx, irdim_v)
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
                                        #=monodromy=# true)
    end

    return splits_bandg
end
function complete_split(bandg :: BandGraph{D}, v::Tuple{LGIrrep{D}, Int}) where D
    lgir, irmul = v
    return complete_split(bandg, (label(lgir), irmul))
end

function _complete_split(
    bandg :: Union{BandGraph{D}, BandGraphPermutations{D}},
    irmul :: Int,
    kidx :: Int,
    iridx :: Int,
    midx :: Int,
    irdim :: Int,
    monodromy :: Bool = false
    ) where D

    partitions = bandg.partitions

    # modify the specific partition that `v = (kidx, iridx, midx)` belongs to by removing
    # it and adding two new "split-off" surrogate vertices in its place
    split_p = build_split_partition(partitions, kidx, iridx, midx, irdim)

    # all subsequent partitions in `bandg` have their `iridxs` incremented by 1 since the
    # degree of the graph is now one larger
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
    subgraphs = bandg isa BandGraphPermutations ? bandg.subgraphs_ps : bandg.subgraphs

    # first, find the parts of `subgraphs` that are affected by the split
    affected_subgraphs_idxs = findall(s -> s.p_max.kidx == kidx, subgraphs)
    affected_subgraphs = subgraphs[affected_subgraphs_idxs]

    row_in_A = partitions[kidx].multiples[midx][irmul]
    affected_subgraphs_As = [Matrix{Int}[] for _ in 1:length(affected_subgraphs)]
    for (i, s) in enumerate(affected_subgraphs)
        A = if bandg isa BandGraphPermutations
            As = s.As
            length(As) == 1 || error("cannot split a subgraph which already has multiple permutations")
            As[1]
        else
            s.A
        end

        cols_in_A = findall(!iszero, A[row_in_A,:])
        A_part = A[row_in_A, cols_in_A]
        if !(length(A_part) == irdim && all(isone, A_part))
            # unexpected nonzero parts of subgraph A: expected `fill(1, irdim)`; vertex is 
            # not fully splittable at `v`; this likely indicates that `v` is connected to
            # a nonmaximal irrep which is not 1-dimensional; return `nothing` as sentinel
            return nothing
            # FIXME: Consider if this really is meaningful behavior in the context of
            #        keyword argument `separable_degree` to `findall_separable_vertices`:
            #        it doesn't seem that meaningful - it's not clear that this is really
            #        even necessary: maybe the vertex should just not be split completely
            #        in this case. Or we should, synthetically, allow a block adjacency 
            #        matrix like [1 0; 0 1; 0 1] for a case where `A_part` = [1,2] and
            #        we are splitting a 3-dimensional irrep. I wonder if that breaks any
            #        of our "sum rules" for band graphs
        end
        
        # now, we replace A = [1,1,...] by a permutation of the identity matrix (e.g., for
        # irdim=2, we replace A = [1,1] by either [1 0; 0 1] or [0 1; 1 0]); if this is the
        # first of the affected subgraphs, we only include the identity matrix, since the
        # inclusion of additional permutations would give graphs that are isomorphic under
        # edge contraction with respect to the split-off surrogate vertices; if it's a
        # monodromy-related vertex, we just permute in any case to be safe
        for (q, perm) in enumerate(permutations(1:irdim))
            if !monodromy
                (i == 1 && q ≠ 1) && continue # skip [0 1; 1 0] cf. above
            end
            if first(s.p_nonmax.klab) == 'Ω'
                # we want specifically to retain a mapping relationship between monodromy-
                # paired irreps that are tied together via Ω<...> since these edges signify
                # that we physically require monodromy-paired irreps to have the same energy
                q ≠ 1 && continue # skip [0 1; 1 0]
            end
            
            top = A[1:row_in_A-1, :]
            bottom = A[row_in_A+1:end, :]
            A_part = I(irdim)[:, perm] # TODO: optimize to avoid explicit creation of A_part
            vcat_to_bottom = zeros(Int, irdim, size(A, 2))
            vcat_to_bottom[:, cols_in_A] .= A_part
            split_A = vcat(top, bottom, vcat_to_bottom)

            push!(affected_subgraphs_As[i], split_A)
        end
    end

    # build a new `split_subgraphs_ps`, using `split_partitions` & `affected_subgraphs_As`
    i′ = 0
    split_subgraphs_ps = Vector{SubGraphPermutations{D}}(undef, length(subgraphs))
    for (i, s) in enumerate(subgraphs)
        As = if i ∉ affected_subgraphs_idxs
            bandg isa BandGraphPermutations ? s.As : [s.A]
        else
            affected_subgraphs_As[i′+=1]
        end
        p_max    = split_partitions[s.p_max.kidx]
        p_nonmax = split_partitions[s.p_nonmax.kidx]
        split_s = SubGraphPermutations(p_max, p_nonmax, As, s.pinned)
        split_subgraphs_ps[i] = split_s
    end

    split_bandg = BandGraphPermutations(split_partitions, split_subgraphs_ps)
    return split_bandg
end

function find_vertex_in_partitions(
            partitions :: Vector{<:Partition}, irlab :: AbstractString, irmul :: Int)
    klab = klabel(irlab)
    for (kidx, p) in enumerate(partitions)
        klab == p.klab || continue

        idxs = findall(lgir -> label(lgir) == irlab, p.lgirs)

        isempty(idxs) && error(lazy"could not find irrep label $irlab at $klab ($p.lgirs)")
        irmul > length(idxs) && error(lazy"irrep index of $irmul exceeds cardinality of irrep $irlab")
        kidx == p.kidx || error("unexpected misalignment of partition `kidx` and its index in the `bandg.partitions` vector")      
        
        iridx = idxs[irmul]
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
function is_valid_split(bandg :: BandGraph)
    partitions = bandg.partitions
    any(p->last(p.klab) == '′', partitions) || return true

    g = assemble_simple_graph(bandg)
    for p′ in Iterators.reverse(partitions)
        klab′ = p′.klab
        last(klab′) == '′' || continue
        
        klab = SubString(klab′, firstindex(klab′), prevind(klab′, lastindex(klab′))) # faster `chop(klab′)` 
        i = findfirst(p -> p.klab == klab, partitions)
        if isnothing(i)
            klabs = [p.klab for p in partitions]
            error(lazy"could not find a monodromy-related partition for $klab′ (trying to match to $klab); available partitions are $klabs")
        end
        p = partitions[i] # "original" partition (not monodromy-translated)

        iridxs′ = p′.iridxs
        iridxs  = p.iridxs

        # now check that there is a path in `g` from `iridx` to `iridx′`
        bool = all(zip(iridxs, iridxs′)) do (iridx, iridx′)
            has_path(g, iridx, iridx′)
        end
        if !bool
            return false
        end
    end
    return true
end