using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie
using KdotP

timereversal = true
weyl_irs_d = Dict{Int, Vector{LGIrrep{3}}}()
for sgnum in 1:230
    sgnum == 1 && empty!(weyl_irs_d)
    lgirsd = lgirreps(sgnum)
    timereversal && (lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirsd))
    
    for lgirs in values(lgirsd)
        isspecial(first(lgirs)) || continue
        for lgir in lgirs
            if KdotP.isweyl(lgir; timereversal)
                # TODO: This appears to miss several space groups that supports Weyls and 
                #       are listed in Table S1 of https://doi.org/10.1016/j.scib.2021.10.023
                haskey(weyl_irs_d, sgnum) || (weyl_irs_d[sgnum] = Vector{String}())
                push!(weyl_irs_d[sgnum], lgir)
            end
        end
    end

    # now we go on to check individual Hilbert base vectors, if space group allows Weyls
    haskey(weyl_irs_d, sgnum) || continue
    irs = weyl_irs_d[sgnum]
    irlabs = label.(irs)

    print("space group: ", sgnum, " (")
    join(stdout, irlabs, ", ")
    println(")")

    subts = subduction_tables(sgnum; timereversal)
    sb, _ = compatibility_basis(sgnum; timereversal)
    checked_vectors = 0
    n_permutations = BigInt(0)
    for _n in sb
        n = SymVector(_n, sb.irlabs, lgirsd)
        any(irlab -> contains(string(n), irlab), irlabs) || continue
        checked_vectors += 1

        bandg = build_subgraphs(n, subts, lgirsd);
        subgraphs_ps = permute_subgraphs(bandg.subgraphs);
        bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps);
        n_permutations += prod(BigInt.(BandGraphs.permutation_counts(bandgp)))
    end

    # subgroup check
    gr = maximal_subgroups(sgnum)
    previously_seen_idxs = findall(num -> num ∈ keys(weyl_irs_d),   @view gr.nums[2:end])
    previously_seen = (@view gr.nums[2:end])[previously_seen_idxs]
    if false
    # TODO: Now check if the Weyl-supporting irrep of supergroup G subduces into one of the
    #       Weyl-supporting irreps of subgroup H; if not, we the group-relation is not 
    #       useful
    previously_seen_kill_idxs = Int[]
    for (i, sgnumᴴ) in enumerate(previously_seen)
        @show sgnumᴴ
        @show getfield.(gr[sgnum].children, Ref(:num))
        idxᴴ = something(findfirst(c->c.num==sgnumᴴ, gr[sgnum].children))
        # FIXME: This fails for G = 210 and H = 24 because the group relation is not maximal
        #        and goes "through" an intermediary subgroup 98. So, one needs to
        #        sequentially apply the transformations. Tedious
        conjcls = gr[sgnum].children[idxᴴ].classes
        length(conjcls)>1 && @warn "multiple conjugacy classes; presently unhandled"
        conjrel = conjcls[1]
        P, p = conjrel.P, conjrel.p
        
        irsᴴ = weyl_irs_d[sgnumᴴ]
        candidate = false
        for irᴴ in irsᴴ
            # do quite a dance in order to make sure the two irreps share coordinates
            # TODO: should probably also check that k-vectors of H and G are equivalent
            #       in transformed setting

            lgᴴ′ = primitivize(group(irᴴ))
            kᴴ′ = position(lgᴴ′)
            irᴴ′ = LGIrrep{3}(label(irᴴ), lgᴴ′, irᴴ.matrices, irᴴ.translations, irᴴ.reality, irᴴ.iscorep)
            for irᴳ in irs
                println("   comparing ", label(irᴳ), " vs. ", label(irᴴ), " (⋕", sgnumᴴ, ")")
                lgᴳ = group(irᴳ)
                kᴳ  = position(lgᴳ)
                kᴳ′ = transform(primitivize(position(lgᴳ), centering(lgᴳ)), P)

                kᴳ′ == kᴴ′ || (println("    k-vectors do not match! skipping"); continue)
                lgᴳ′_ops = primitivize.(transform.(lgᴳ, Ref(P), Ref(p)), centering(group(irᴳ)))
                lgᴳ′ = LittleGroup{3}(sgnum, kᴳ′, klabel(lgᴳ), lgᴳ′_ops)
                irᴳ′ = LGIrrep{3}(label(irᴳ), lgᴳ′, irᴳ.matrices, irᴳ.translations, irᴳ.reality, irᴳ.iscorep)

                if !Crystalline._findsubgroup(lgᴳ′, lgᴴ′)[1] # we do not have H < G
                    println("   little groups do not match! skipping")
                    # FIXME: This seems a little overzealous, at least in combination with 
                    #        the k-vector check: this way, we seem to be missing out on
                    #        some relations, e.g., for G = 178 and H = 152 where there
                    #        should be a match between K₃ (178) and KA₃ (152), I think.
                    continue
                end
                if subduction_count(irᴳ′, irᴴ′) ≠ 0
                    println("    ", label(irᴳ), " ↓ ", label(irᴴ), " (⋕", sgnumᴴ, ")")
                    candidate=true
                end
            end
        end
        if !candidate
            push!(previously_seen_kill_idxs, i)
            println("no relevant irrep subductions in ⋕", sgnumᴴ)
        end
    end
    deleteat!(previously_seen_idxs, previously_seen_kill_idxs)
    end

    # print info
    println("  Hilbert vectors to check   : ", checked_vectors)
    println("  total permutations to check: ", n_permutations)
    if !isempty(previously_seen)
    print("  supergroup of              : ")
        join(stdout, previously_seen, ", ")
        println()
    end
    println()
end

    
