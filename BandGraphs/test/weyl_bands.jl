using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie: plot, save, Figure, Axis, Makie
#using CairoMakie: plot, save, Figure, Axis, Makie
using KdotP
using LinearAlgebra: norm, dot

timereversal = true

weyl_irs_d = Dict{Int, Vector{LGIrrep{3}}}()
checked_vectors_d = Dict{Int, Int}()
n_permutations_d = Dict{Int, BigInt}()
for sgnum in 1:230
    sgnum == 1 && empty!(weyl_irs_d)
    # --- identify space groups with Weyl points ---
    lgirsd = lgirreps(sgnum)
    timereversal && realify!(lgirsd)
    for lgirs in values(lgirsd)
        isspecial(first(lgirs)) || continue
        for lgir in lgirs
            if KdotP.isweyl(lgir; timereversal)
                haskey(weyl_irs_d, sgnum) || (weyl_irs_d[sgnum] = Vector{String}())
                push!(weyl_irs_d[sgnum], lgir)
            end
        end
    end
    haskey(weyl_irs_d, sgnum) || continue # skip further work if Weyl points not supported
    irlabs = label.(weyl_irs_d[sgnum])

    # --- find Hilbert basis vectors w/ Weyl points ---
    subts = subduction_tables(sgnum; timereversal)
    sb, _ = compatibility_basis(sgnum; timereversal)
    checked_vectors = 0
    n_permutations = BigInt(0)
    for _n in sb
        n = SymmetryVector(_n, sb.irlabs, lgirsd)
        any(irlab -> contains(string(n), irlab), irlabs) || continue
        checked_vectors += 1

        bandg = build_subgraphs(n, subts, lgirsd)
        subgraphs_ps = permute_subgraphs(bandg.subgraphs)
        bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps)
        n_permutations += length(bandgp)
    end
    checked_vectors_d[sgnum] = checked_vectors
    n_permutations_d[sgnum] = n_permutations
end


# --- subgroup check ---
# check whether the Weyl-supporting irrep of supergroup G subduces into one of the Weyl-
# supporting irreps of a subgroup H; if not, the group-relation is not useful
weyl_gr = MetaGraph(DiGraph(), Int, Vector{LGIrrep{3}}, Vector{Vector{LGIrrep{3}}})
for (sgnum, lgirs) in weyl_irs_d
   add_vertex!(weyl_gr, sgnum, lgirs)
end

for sgnumᴳ in 1:230
    haskey(weyl_irs_d, sgnumᴳ) || continue
    
    # --- basic check of subgroup compatibility ---
    gr = maximal_subgroups(sgnumᴳ)
    previously_seen_idxs = findall(gr[sgnumᴳ]) do rᴴ
        sgnumᴴ = rᴴ.num
        sgnumᴴ ∈ keys(weyl_irs_d)
    end
    previously_seen = getfield.(gr[sgnumᴳ][previously_seen_idxs], Ref(:num))

    # --- print info ---
    printstyled("G: ", sgnumᴳ, " (", join(label.(weyl_irs_d[sgnumᴳ]), ", "), ")\n";
                bold=true, color=:blue)
    println("   vecs (permutations) : ", 
            checked_vectors_d[sgnumᴳ], " (", n_permutations_d[sgnumᴳ], ")")
    if !isempty(previously_seen)
        print("   supergroup of       : ")
        join(stdout, previously_seen, ", ")
        println()
    else
        printstyled("   base-case\n"; color=:red)
    end

    # --- more careful subgroup/irrep/k-vector check ---
    irsᴳ = weyl_irs_d[sgnumᴳ]
    sgᴳ = reduce_ops(spacegroup(sgnumᴳ, Val(3)), centering(sgnumᴳ, 3))
    for (i, sgnumᴴ) in enumerate(previously_seen)
        printstyled("   G<H: ", sgnumᴳ, "<", sgnumᴴ, "\n"; color=:blue, bold=true)
        
        rels = Crystalline.conjugacy_relations(gr, sgnumᴳ, sgnumᴴ)
        
        irsᴴ = weyl_irs_d[sgnumᴴ]
        candidate = false
        for irᴴ in irsᴴ
            # --- check whether H and G irreps are compatible ---
            # (i.e., that irᴳ could subduce to irᴴ)

            lgᴴ = group(irᴴ) # already has no redundant copies (cf. `littlegroups` contract)
            kᴴ = position(lgᴴ)
            cntrᴴ = centering(lgᴴ)
            for (j, irᴳ) in enumerate(irsᴳ)
                print("      ", label(irᴳ), " vs. ", label(irᴴ), ": ")
                lgᴳ = group(irᴳ)
                kᴳ  = position(lgᴳ)
                kᴳs = sgᴳ .* Ref(kᴳ) # orbit of kᴳ (including _all_ equivalent G-points)
                for t in rels
                    P, p = t.P, t.p
                    
                    # check that k-vectors of H and G are equivalent in transformed setting;
                    # the check must encompass the full orbit of one of the k-vectors,
                    # including equivalent reciprocal points, since we may need to change
                    # the irreps of G accordingly (for nonsymmorphic groups)
                    kᴳᴴs = [transform(_kᴳ, P) for _kᴳ in kᴳs] # kᴳ points referred to conventional H basis (Gsᴴ)
                    idx_kᴳᴴs_eq_kᴴ = findfirst(kᴳᴴ->isapprox(kᴳᴴ, kᴴ, nothing, false), kᴳᴴs)
                    if isnothing(idx_kᴳᴴs_eq_kᴴ)
                        # kᴳ is not equivalent by either symmetry or reciprocal lattice
                        # vectors to kᴴ; we cannot compare the little groups
                        printstyled("÷", color=:red)
                        printstyled(" (k-points)"; color=:light_black)
                        continue
                    end
                    idx_kᴳᴴs_eq_kᴴ = something(idx_kᴳᴴs_eq_kᴴ)
                    q_coset = sgᴳ[idx_kᴳᴴs_eq_kᴴ] # relevant coset representative
                    kᴳ′ᴴ = kᴳᴴs[idx_kᴳᴴs_eq_kᴴ]

                    # since kᴳ might not equal kᴴ, we might need to remap the operations in
                    # irrep according to the coset-representative that connects kᴳ to kᴴ;
                    # but we need to compare the two k-points in a shared setting; we choose
                    # to first map kᴴ to the setting of kᴳ, calling this kᴴᴳ
                    kᴴᴳ = transform(kᴴ, inv(P)) # kᴴ in G basis
                    tmp = Crystalline.remap_lgirreps_to_point_in_kstar([irᴳ], kᴴᴳ, [q_coset])
                    irᴳ′ = tmp[1]
                    lgᴳ′ = group(irᴳ′)

                    # now transform to the basis of H (lgᴳ ops in H setting at kᴴ position)
                    lgᴳ′ᴴ_ops = transform.(lgᴳ′, Ref(P), Ref(p), false)

                    # the operator-ordering of the H and G irreps might differ; find the 
                    # permutation between the two (or, if their little groups are not
                    # subgroups, bail out)
                    boolsub, idxsᴳ²ᴴ = Crystalline._findsubgroup(lgᴳ′ᴴ_ops, operations(lgᴴ), cntrᴴ)
                    if !boolsub # we do not have H < G
                        printstyled("÷", color=:red)
                        printstyled(" (little groups)"; color=:light_black)
                        continue
                    end

                    # even if we now have a subgroup, `lgᴳ′ᴴ` and `lgᴴ` may still differ
                    # in their operations by primitive lattice vectors; and this may change
                    # the irreps of nonsymmorphic groups; so we have to correct for that
                    # by identifying possible translation mismatches between operations
                    # and multiplying by the appropriate Bloch phase
                    Δts = translation.(lgᴴ) .- translation.(lgᴳ′ᴴ_ops[idxsᴳ²ᴴ])
                    translationsᴳ′ = [copy(τ) for τ in irᴳ′.translations]
                    #matricesᴳ′ = [copy(m) for m in irᴳ′.matrices]
                    for (i, Δt) in enumerate(Δts)
                        norm(Δt) > Crystalline.DEFAULT_ATOL || continue
                        iᴳ = idxsᴳ²ᴴ[i]
                        lgᴳ′ᴴ_ops[iᴳ] = SymOperation(Δt) * lgᴳ′ᴴ_ops[iᴳ]
                        translationsᴳ′[iᴳ] += Δt
                        #matricesᴳ′[iᴳ] .*= cispi(2*dot(kᴳ′ᴴ.cnst, Δt))
                        # TODO: do something similar w/ irᴳ.translations & kᴳ′ᴴ.free
                    end
                    lgᴳ′ᴴ = LittleGroup{3}(sgnumᴳ, kᴳ′ᴴ, klabel(lgᴳ), lgᴳ′ᴴ_ops[idxsᴳ²ᴴ])

                    irᴳ′ᴴ = LGIrrep{3}(label(irᴳ), lgᴳ′ᴴ, irᴳ.matrices[idxsᴳ²ᴴ],
                                       translationsᴳ′[idxsᴳ²ᴴ], irᴳ.reality, irᴳ.iscorep)

                    if subduction_count(irᴳ′ᴴ, irᴴ) ≠ 0
                        printstyled("✓"; color=:green)
                        candidate=true

                        if !has_edge(weyl_gr, code_for(weyl_gr, sgnumᴳ), code_for(weyl_gr, sgnumᴴ))
                            edge_data = [Vector{LGIrrep{3}}() for _ in 1:length(irsᴳ)]
                            push!(edge_data[j], irᴴ)
                            add_edge!(weyl_gr, sgnumᴳ, sgnumᴴ, edge_data)
                        else
                            push!(weyl_gr[sgnumᴳ, sgnumᴴ][j], irᴴ)
                        end
                    else
                        printstyled("÷", color=:red)
                        printstyled(" (subduction)"; color=:light_black)
                    end
                end
                println()
            end
        end
        if !candidate
            printstyled("      no relevant irrep subductions\n"; color=:red)
        end
    end
end

## -----------------
function layout_y_position(weyl_gr)
    sgnums = collect(labels(weyl_gr))
    orders = Crystalline.SG_PRIMITIVE_ORDERs[3][sgnums]
    layers′ = CrystallineGraphMakieExt.map_numbers_to_oneto(orders; rev = false)
    maxlayer = maximum(layers′)
    1:length(orders) .=> (maxlayer+1) .- layers′
end

const CrystallineGraphMakieExt = Base.get_extension(Crystalline, :CrystallineGraphMakieExt)
const Zarate = CrystallineGraphMakieExt.Zarate
const solve_positions = CrystallineGraphMakieExt.solve_positions
function zarate_layout(weyl_gr′, force_y_layer)
    # `solve_positions` implicitly assumes `SimpleGraph`/`SimpleDiGraph` type graphs, so we
    # must first convert `gr` to `gr′::SimpleDiGraph`
    gr′ = SimpleDiGraphFromIterator(edges(weyl_gr′))
    xs, ys, _ = solve_positions(Zarate(), gr′; force_layer=force_y_layer)
    # TODO: This is broken in that it does not respect the `force_y_layer` if it is not
    #       consecutive - so cannot be compared across disconnected components
    return GraphMakie.Point2f.(ys, -xs)
end


_tmp = induced_subgraph.(Ref(weyl_gr), connected_components(weyl_gr))
weyl_connected_grs = getindex.(_tmp, 1)
weyl_connected_idxs = getindex.(_tmp, 2)
y_layers = layout_y_position(weyl_gr)
connected_y_layers = [1:length(idxs) .=> last.(y_layers[idxs]) for idxs in weyl_connected_idxs]
f = Figure()
N_connected = length(weyl_connected_grs)
for (n, weyl_gr′) in enumerate(weyl_connected_grs)
    nlabels_text = map(vertices(weyl_gr′)) do i # label for each vertex
        sgnum = label_for(weyl_gr′, i)
        irlabs = join(label.(weyl_irs_d[sgnum]), ", ")
        string(sgnum) * "\n(" * irlabs * ")"
    end
    nlabels_color = map(vertices(weyl_gr′)) do i
        sgnum = label_for(weyl_gr′, i)
        isempty(outneighbors(weyl_gr′, i)) ? :red : :black
    end
    elabels_text = map(edges(weyl_gr′)) do e
        sgnumᴳ, sgnumᴴ = label_for(weyl_gr′, src(e)), label_for(weyl_gr′, dst(e))
        irlabsᴳ = label.(weyl_irs_d[sgnumᴳ])
        irsᴴv = weyl_gr′[sgnumᴳ, sgnumᴴ]
        irlabsᴴ = map(irsᴴv) do irsᴴ
            irlabsᴴ = label.(irsᴴ)
            length(irlabsᴴ) == 1 ? irlabsᴴ[1] : "(" * join(irlabsᴴ, ", ") * ")"
        end
        join([irlabsᴳ[i] * " ↓ " * irlabsᴴ[i] for i in 1:length(irlabsᴳ)], "\n")
    end

    ax = Axis(f[1, n])
    p = graphplot!(ax, weyl_gr′;
        nlabels = nlabels_text,
        nlabels_align = (:center, :bottom),
        nlabels_offset = GraphMakie.Point(0, .025),
        nlabels_color = nlabels_color,
        nlabels_fontsize = 13,
        elabels = elabels_text,
        elabels_offset = GraphMakie.Point(0, .025),
        elabels_fontsize = 13,)

    # adjust appearance of graph
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    xy′ = zarate_layout(weyl_gr′, connected_y_layers[n])
    p[:node_pos][] = xy′ # update layout of vertices
    y_extrema = (-1) .* reverse(extrema(last.(y_layers)))
    y_extrema = y_extrema .+ (0.1 * (-)(y_extrema...)) .* (1, -1)
    Makie.ylims!(y_extrema...)
end
f