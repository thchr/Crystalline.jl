@enum Separability begin
    SEPARABLE
    INSEPARABLE_NO_MET_IRREP_CRITERIA
    INSEPARABLE_CANNOT_SPLIT
    INSEPARABLE_CANNOT_GROUP
    INFEASIBLE_TOO_MANY_BANDGRAPH_PERMUTATIONS
    INFEASIBLE_TOO_MANY_SPLIT_PERMUTATOINS
    INFEASIBLE_TOO_MANY_SUBBLOCK_PERMUTATIONS
end

is_separable(s::Separability) = s == SEPARABLE
function is_feasible(s::Separability)
    return (s != INFEASIBLE_TOO_MANY_BANDGRAPH_PERMUTATIONS &&
            s != INFEASIBLE_TOO_MANY_SPLIT_PERMUTATOINS)
end

struct SeparabilityState
    s :: Separability
    permutations :: Union{Nothing, BigInt}
end
SeparabilityState(s::Separability)  = SeparabilityState(s, nothing)
is_separable(sepstate::SeparabilityState) = is_separable(sepstate.s)
is_feasible(sepstate::SeparabilityState)  = is_feasible(sepstate.s)

"""
Find all the vertices in a band graph `bandg` whose associated irrep `lgir` for which
`criterion(lgir)` is true. 

Return a vector of 2-tuples, whose first element is a valid `lgir`, and whose second element
is a linear index into the multiplicity of this irrep in `bandg`.
"""
function findall_vertices(criterion::F, bandg :: BandGraph{D}) where {F,D}
    vs = Vector{Tuple{LGIrrep{D}, Int}}()
    for p in bandg.partitions
        # only include non-monodromy partitions, since monodromy irreps will be split
        # with their parents anyway; i.e., don't double count
        is_monodromy = endswith(p.klab, "′")
        is_monodromy && continue

        # only consider splitting at maximal vertices
        p.maximal || continue
            
        for ms in p.multiples
            for (i, m) in enumerate(ms)
                lgir = p.lgirs[m]
                if criterion(lgir) :: Bool
                    push!(vs, (lgir, i))
                end
            end
        end
    end
    return vs
end

"""
    has_vertex(criterion, n; ignore_unit_occupation=true) --> Bool

Test if a symmetry vector `n` includes an irrep that fulfils the criterion `criterion`,
returning `true` if it does, and `false` otherwise. `criterion` is a function or functor
that takes an `LGIrrep` and returns a boolean.

If `ignore_unit_occupation` is `true` (default), then the function will return `false` if
the occupation of the symmetry vector is 1, regardless of `criterion`.
"""
function has_vertex(
        criterion :: F,
        n :: AbstractSymmetryVector{D};
        ignore_unit_occupation :: Bool = true
    ) where {F,D}
    # usually nothing interesting in 1-band cases; skip unless `ignore_unit_occupation = 0`
    ignore_unit_occupation && occupation(n) == 1 && return false

    for (lgirs, mults) in zip(irreps(n), multiplicities(n))
        for (lgir, m) in zip(lgirs, mults)
            m ≠ 0 && criterion(lgir) && return true
        end
    end
    return false
end

function findall_separable_vertices(
            criterion::F, n::AbstractSymmetryVector{D}, subts, lgirsd; kws...) where {F,D}
    bandg = build_subgraphs(n, subts, lgirsd)
    return findall_separable_vertices(criterion, bandg; kws...)
end
function findall_separable_vertices(
            criterion :: F,
            bandg :: BandGraph{D};
            max_permutations :: Union{Nothing, <:Real} = 1e5,
            separable_degree :: Union{Nothing, <: Integer} = nothing,
            max_subblock_permutations :: Union{Nothing, <:Real} = 1e6
        ) where {F,D}
    
    separable = Tuple{BandGraph{D}, BandGraph{D}, LGIrrep{D}}[]

    # find which vertices fulfil criterion; i.e., which vertices to check seperability for
    vs = findall_vertices(lgir -> isspecial(lgir) && (criterion(lgir) :: Bool), bandg)
    isempty(vs) && return SeparabilityState(INSEPARABLE_NO_MET_IRREP_CRITERIA), nothing

    # fast-path "necessary conditions" check which does not require permutation construction
    unfeasible_vs_idxs = fast_path_separability_irdim_infeasibility(bandg, vs, separable_degree)
    deleteat!(vs, unfeasible_vs_idxs)
    isempty(vs) && return SeparabilityState(INSEPARABLE_CANNOT_GROUP), nothing
    
    # move on to additional conditions which require us to look at individual permutations
    subgraphs_ps = permute_subgraphs(bandg.subgraphs; max_subblock_permutations)
    if isnothing(subgraphs_ps)
        return SeparabilityState(INFEASIBLE_TOO_MANY_SUBBLOCK_PERMUTATIONS), nothing
    end
    bandgp = BandGraphPermutations(bandg.partitions, something(subgraphs_ps))

    # if there are too many permutations to consider, we bail out
    Nbandgp = safe_length(bandgp)
    if !isnothing(max_permutations) && Nbandgp > max_permutations
        sepstate = SeparabilityState(INFEASIBLE_TOO_MANY_BANDGRAPH_PERMUTATIONS, Nbandgp)
        return sepstate, nothing
    end

    # okay, now we actually have to do the "full" checks for any remaining vertices by
    # explicitly considering each possible permutation one by one
    success_vs = zeros(Bool, length(vs))
    work_g = Graph(nv(bandg))
    work_split_g_for_vs = Graph.((nv(bandg) - 1) .+ irdim.(first.(vs))) # one for each `vs` (in case they have different `irdim`)
    for bandg′ in bandgp # iterate over graph permutations
        may_be_articulation_vs = may_be_articulation_vertices(bandg′, vs, work_g)
        for (i, v) in enumerate(vs)
            # constrain ourselves to return only a single representative seperable band
            # permutation if any such permutation exists
            success_vs[i] && continue
            # v must be an articulation point to be a separable point (necessary condition)
            may_be_articulation_vs[i] || continue

            lgir, irmul = v
            split_bandg = is_separable_at_vertex_check_all_splits(
                            bandg′, (lgir, irmul), work_g, work_split_g_for_vs[i];
                            verbose=true, separable_degree)
            if !isnothing(split_bandg)
                push!(separable, (bandg′, split_bandg, lgir))
                success_vs[i] = true
            end
        end
        if all(success_vs)
            # all vertices have been split; no need to look at more graph permutations
            break
        end
    end

    if isempty(separable)
        return SeparabilityState(INSEPARABLE_CANNOT_SPLIT), nothing
    else
        return SeparabilityState(SEPARABLE), separable
    end
end

function fast_path_separability_irdim_infeasibility(
    bandg::BandGraph{D},
    vs::Vector{Tuple{LGIrrep{D}, Int}},
    separable_degree :: Union{Nothing, <:Integer} = nothing
    ) where D
    # - for `v` to be splittable, it must be possible to separate the irreps of every
    #   partition into groups that contain a single split-off surrogate vertex `vˣ`, with
    #   each such group having a consistent number of bands; this is precisely the job
    #   of `solve_subset_sum_variant` to decide if this is possible (but there, across
    #   bands over subsets of partitions, not over individual partitions, as here).
    # - this is equivalent to ignoring the additional constraints from compatibility
    #   relations; but since they are _additional_ this is a still a necessary condition.
    # - if `v` has a monodromy partner, we simply remove that partition from the
    #   consideration, since the monodromy partner would be split up in exactly the same
    #   way
    # - we return a list of indices `infeasible_idxs` into `vs`, with each index signifying
    #   a vertex that failed this test

    partitions_irdims = [irdim.(p.lgirs) for p in bandg.partitions]
    infeasible_idxs = Int[]
    for (i, v) in enumerate(vs)
        lgir, irmul = v
        irlab = label(lgir)
        kidx, iridx, _ = something(find_vertex_in_partitions(bandg.partitions, irlab, irmul))
        v_irdim = partitions_irdims[kidx][iridx]

        # check if there are any monodromy-related partners to `v`; if so, don't include the
        # monodromy-related partition in the set of "exterior" partitions
        klab′ = klabel(irlab) * '′'
        kidx′ = findfirst(p->p.klab == klab′, bandg.partitions)

        # now see if there is a possible irrep-grouping across partitions exists that allow
        # for a single split-off surrogate vertex in each group; collect irrep dimensions 
        # "outside" the one belong to `v` in `Sʲs`, and those belong to `v` in `T` and `R`
        # with `T` denoting the irrep dimensions of the split-off surrogate vertices `vˣ`
        # and `R` denoting the remaining irreps in the partition where `v` lives
        Sʲs = if isnothing(kidx′) # no monodromy partner
            [Sʲ for (j, Sʲ) in enumerate(partitions_irdims) if j ≠ kidx]
        else
            [Sʲ for (j, Sʲ) in enumerate(partitions_irdims) if j ∉ (kidx, kidx′)]
        end
        T = fill(1, v_irdim) # TODO: maybe someday a more general thing based on non-maximal irrep multiplicities incident on `v`?
        R = Int[v′_irdim for (iridx′, v′_irdim) in enumerate(partitions_irdims[kidx]) if iridx′≠iridx]

        feasible, _ = if isnothing(separable_degree) || separable_degree == irdim(lgir)
            solve_subset_sum_variant(Sʲs, T, R)
        else
            solve_subset_sum_variant_flexibleT(Sʲs, T, R, separable_degree)
        end
        if !feasible
            push!(infeasible_idxs, i)
        end
    end

    infeasible_idxs
end

function may_be_articulation_vertices(
        bandg :: BandGraph{D}, 
        vs :: Vector{Tuple{LGIrrep{D}, Int}},
        g :: Graph = Graph(nv(bandg)) # work-graph (mutated)
    ) where D
    # --- fast path check ---
    # "fast-path" articulation point check: a necessary (but insufficient) condition is that
    # a vertex must be an articulation point in order to be a separable point; it is 
    # advantageous to check this immediately, since we then don't need to explore any of the
    # (possibly many) splits of `bandg`, only `bandg` itself. However, we can only use this
    # if `v` doesn't have a monodromy partner (if it does, the monodromy partner would
    # prevent it from being a genuine articulation point (it would require removal of both))
    may_be_articulation_vs = zeros(Bool, length(vs))
    for (i, v) in enumerate(vs)
        lgir = first(v)
        klab′ = klabel(lgir) * '′'
        if any(p->p.klab == klab′, bandg.partitions)
            # monodromy partner case: we can't tell, so it _may_ effectively, still be an
            # articulation point (in the sense that removal of both vertices would cut)
            may_be_articulation_vs[i] = true
        end
    end
    # if all are monodromy-paired; no point in checking further
    all(may_be_articulation_vs) && return may_be_articulation_vs

    # okay; some non-monodromy-paired vertices; get actual articulation vertices of `bandg`
    g = assemble_simple_graph!(g, bandg; reset_edges = true)
    articulation_vertex_idxs = articulation(g)

    # check which `vs` are among articulation vertices of `bandg`
    for (i, v) in enumerate(vs)
        may_be_articulation_vs[i] && continue # monodromy-paired; skip further checks

        lgir, irmul = v
        irlab = label(lgir)
        
        # not monodromy-paired; can do the check
        kidx, iridx, _ = something(find_vertex_in_partitions(bandg.partitions, irlab, irmul))

        p = bandg.partitions[kidx]
        vertex_idx = p.iridxs[iridx] # "global" vertex index into `g` below
        may_be_articulation_vs[i] = vertex_idx ∈ articulation_vertex_idxs
    end

    return may_be_articulation_vs
end

function is_separable_at_vertex_check_all_splits(
        bandg :: BandGraph{D}, 
        v :: Tuple{LGIrrep{D}, Int},
        g :: Graph = Graph(nv(bandg)),                         # `bandg` work graph (mutated)
        split_g :: Graph = Graph(nv(bandg) + irdim(v[1]) - 1); # `split_bandg` work graph (mutated)
        verbose :: Bool = true,
        separable_degree :: Union{Nothing, <:Integer} = nothing
    ) where D

    # The function determines whether a graph G is separable at a vertex v, as described in 
    # our Overleaf document; before doing "heavy-lifting", we do a fast-path "necessary
    # conditions" check, to potentially avoid actually creating any explicit band splits
    # if the vertex is clearly not separable

    split_lgir = first(v)

    # --- "slow-path" checks ---
    # TODO: Return and pass through `INFEASIBLE_TOO_MANY_SPLIT_PERMUTATIONS` if there are 
    #       too many splits    
    splits_bandg = complete_split(bandg, v; separable_degree)
    isnothing(splits_bandg) && return nothing # no valid splits at `v`

    # set how many components the band graph should be separated into after the split, in
    # order for us to classify it is a "succesful" split
    split_irdim = irdim(split_lgir)
    if isnothing(separable_degree)
        # "full" separability: into as many distinct band graphs as the split-irrep's dimen.
        separable_degree = split_irdim
    else
        # "custom": any number less than or equal to `split_irdim`
        if split_irdim < separable_degree
            error(lazy"`separable_degree` (=$separable_degree) cannot exceed the splitted irrep's dimension (=$split_irdim)")
        end
    end

    # obtain the connected components of the split-up graph, across both partitions and
    # energy-separation: `bandg_cs[partitions-grouping][energy-grouping]`
    g = assemble_simple_graph!(g, bandg; reset_edges = true)
    bandg_cs = group_partition_connected_components(bandg, g) # components

    # initialize `is_valid_split` work buffers
    work_valid_next = Vector{Int}(undef, nv(split_g))
    work_valid_seen = zeros(Bool, nv(split_g))
    if isnothing(bandg_cs)
        # CASE 1: original band graph is fully connected - easy to check separability
        #         (all distinct components after the split will cover the same set of 
        #          partitions)
        # initialize `count_connected_components` work buffers
        work_countc_search_queue = Int[]
        work_countc_label = zeros(Int, nv(split_g))
        # iterate over split permutations
        for split_bandg in splits_bandg 
            split_g = assemble_simple_graph!(split_g, split_bandg; reset_edges = true)
            if !is_valid_split(split_bandg, split_g, work_valid_next, work_valid_seen)
                #verbose && printstyled("   ... skipping an invalid split\n"; color=:yellow)
                verbose && error("encountered an invalid split - that's new!")
                continue
            end
            # get conventional graph for `split_bandg` & get its number of connected
            # components; in fact, we do not need the edge weights, since connectivity is
            # decided only by presence/absence of an edge, so we just build & use the simple
            # graph (instead of the weighted and labelled `assemble_graph`) for efficiency
            split_N = count_connected_components(
                split_g, work_countc_label, work_countc_search_queue; reset_label=true)
            if split_N == separable_degree
                return split_bandg
            end
        end
    else
        # CASE 2: original band graph is not fully connected - harder to check separability
        #         (the distinct components may involve different partitions)
        split_surrogate_irlab = label(split_lgir) * "ˣ" # label for irrep after split
        for split_bandg in splits_bandg # iterate over split permutations
            split_g = assemble_simple_graph!(split_g, split_bandg; reset_edges = true)
            if !is_valid_split(split_bandg, split_g, work_valid_next, work_valid_seen)
                #verbose && printstyled("   ... skipping an invalid split\n"; color=:yellow)
                verbose && error("encountered an invalid split - that's new!")
                continue
            end
            split_bandg_cs = group_partition_connected_components(split_bandg, split_g)
            ijs = findall_vertices_in_components(split_bandg_cs, split_surrogate_irlab)

            # after split, there must be at least `separable_degree` components that
            # each contain a "split surrogate" vertex; if not, no need to look further
            length(ijs) == separable_degree || continue

            # now we must check if it is possible to group the "split-off" components into
            # isolated band graphs, each with a consistent occupation number
            i = first(first(ijs)) # index into `partition-group` in `split_bandg_cs`
            @assert all(ij -> ij[1] == i, ijs)
            T = [occupation(split_bandg_cs[i][j]) for (i,j) in ijs]
            i_components = split_bandg_cs[i]
            js = getindex.(ijs, 2)
            R = [occupation(i_components[j′]) for j′ in 1:length(i_components) if j′ ∉ js]
            Sʲs = [[occupation(split_bandg_c) for split_bandg_c in split_bandg_cs[i′]] 
                                                for i′ in 1:length(split_bandg_cs) if i′≠i]
            feasible, _ = solve_subset_sum_variant(Sʲs, T, R; verbose=false)
            if feasible
                # it's feasible to group the split-off vertices into their own distinct and
                # isolated band graphs
                return split_bandg
            end
        end
    end   

    return nothing
end

function Graphs.induced_subgraph(bandg :: BandGraph{D}, vlist) where D
    partitions = bandg.partitions
    issorted(vlist) || (vlist = sort(vlist))

    # determine which partitions are correspond to the vertex indices in `vlist`
    # (this aligns with the vertex-indexing ordering from `assemble_graph`)
    kidx = 1
    irmax = 0
    code = 0
    partitions′ = empty(partitions)
    vmap = Int[]
    lgirs_idxs = Int[]
    for p in partitions
        for j in eachindex(p.lgirs)
            code += 1
            if insorted(code, vlist)
                push!(lgirs_idxs, j)
                push!(vmap, code)
            end
        end
        isempty(lgirs_idxs) && continue
        lgirs = p.lgirs[lgirs_idxs]
        multiples = _identify_multiples(lgirs, first(lgirs))
        p′ = Partition{D}(p.klab, lgirs, multiples, p.maximal, kidx, irmax+1:irmax+length(lgirs))
        push!(partitions′, p′)

        kidx += 1
        irmax = last(p′.iridxs)
        empty!(lgirs_idxs)
    end
    
    # determine the subgraphs connected by elements in vlist (exploit that this is a
    # bipartite graph)
    subgraphs = bandg.subgraphs
    subgraphs′ = empty(subgraphs)
    for subgraph in subgraphs
        klab_max    = subgraph.p_max.klab
        klab_nonmax = subgraph.p_nonmax.klab
        idx_max = findfirst(p-> p.klab==klab_max, partitions′)
        isnothing(idx_max) && continue
        idx_nonmax = findfirst(p -> p.klab==klab_nonmax, partitions′)
        if isnothing(idx_nonmax)
            error("could not find expected nonmaximal partition from k-point $klab_nonmax")
        end

        # determine include parts of the subgraph adjacency matrix
        p_max′    = partitions′[idx_max]
        p_nonmax′ = partitions′[idx_nonmax]
        A = subgraph.A
        A′ = Matrix{Int}(undef, length(p_max′.lgirs), length(p_nonmax′.lgirs))
        for (local_idx_max′, iridx_max′) in enumerate(p_max′.iridxs)
            for (local_idx_nonmax′, iridx_nonmax′) in enumerate(p_nonmax′.iridxs)
                iridx_max = vmap[iridx_max′]
                iridx_nonmax = vmap[iridx_nonmax′]
                local_idx_max = something(findfirst(==(iridx_max), subgraph.p_max.iridxs))
                local_idx_nonmax = something(findfirst(==(iridx_nonmax), subgraph.p_nonmax.iridxs))
                A′[local_idx_max′, local_idx_nonmax′] = A[local_idx_max, local_idx_nonmax]
            end
        end
        subgraph′ = SubGraph(p_max′, p_nonmax′, A′)
        push!(subgraphs′, subgraph′)
    end
    
    bandg′ = BandGraph(subgraphs′, partitions′)
    return bandg′, vmap
end

function _identify_multiples!(multiples :: Vector{UnitRange{Int}}, lgirs, lgir_target, start)
    # we assume that `lgirs` is sorted in the sense that identical `lgirs` are adjacent
    stop = start
    for j in start+1:length(lgirs)
        lgir = lgirs[j]
        if label(lgir) == label(lgir_target)
            stop = j
        else
            push!(multiples, start:stop)
            return _identify_multiples!(multiples, lgirs, lgirs[j], j)
        end
    end
    push!(multiples, start:stop)
    return multiples
end
function _identify_multiples(lgirs, lgir_target, start::Int=1)
    _identify_multiples!(UnitRange{Int}[], lgirs, lgir_target, start)
end

# return the connected components of `bandg`; returned components are a 
# `Vector{Vector{BandGraph}}` with indexing `bandg_cs[partitions-grouping][energy-grouping]`
# if there is only one connected component (fully connected), returns `nothing` as sentinel
function group_partition_connected_components(
        bandg :: BandGraph{D},
        g :: Graph = assemble_simple_graph(bandg)
        ) where D
    cs = connected_components(g)
    if length(cs) == 1
        # `nothing` as sentinel for being fully connected (i.e., 1 connected component)
        return nothing
    end

    klabs_cs = Set{String}[]
    bandg_cs = Vector{BandGraph{D}}[]
    for c in cs
        bandg_c, vmap = induced_subgraph(bandg, c)
        klabs_c = Set(p.klab for p in bandg_c.partitions)
        idx = findfirst(==(klabs_c), klabs_cs)
        if isnothing(idx)
            push!(klabs_cs, klabs_c)
            push!(bandg_cs, [bandg_c])
        else
            push!(bandg_cs[idx], bandg_c)
        end
    end

    return bandg_cs
end

function findall_vertices_in_components(
            bandg_cs :: Vector{Vector{BandGraph{D}}}, irlab :: String) where D
    klab = klabel(irlab)
    contains_vertices_idxs = Tuple{Int, Int}[]
    for (i, bandgs) in enumerate(bandg_cs) # over each k-connected partition
        for (j, bandg_c) in enumerate(bandgs) # over (energy-separated) components of a k-connected partition
            is_matching_partition_group = false
            partitions_c = bandg_c.partitions
            for partition_c in partitions_c
                if partition_c.klab == klab
                    is_matching_partition_group = true
                    if any(lgir->label(lgir)==(irlab), partition_c.lgirs)
                        push!(contains_vertices_idxs, (i,j))
                    end
                    break # found what we looked for; go to next graph
                else
                    continue # keep looking
                end
            end
            is_matching_partition_group || break # wrong partitions; skip related graphs
        end
    end
    return contains_vertices_idxs
end