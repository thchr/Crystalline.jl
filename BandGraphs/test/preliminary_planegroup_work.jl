using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using SymmetryBases
using GraphMakie
using GLMakie

# BEWARE: There are many choices of `sgnum²ᴰ` for which this does not work, because the 
#         2D cut might not be oriented along z or because the irreps have different labels
#         etc. It does work for a subset of the plane groups though, e.g., plane groups
#         ⋕6 (p2mm), ⋕11 (p4mm), ⋕14 (p3m1; although it has some visualization troubles),
#         ⋕15 (p31m), ⋕17 (p6mm), ⋕12 (p4gm). But seems to fail e.g., for other cases of
#         interest such as ⋕4 (p1g1), ⋕5 (c1m1), and ⋕9 (c2mm) (seems mainly due to
#         different k-vector labels - but could be other things also).

# Example for p2mm (plane group ⋕6)
timereversal = true
sgnum²ᴰ = 12
sgnum³ᴰ = Crystalline.PLANE2SPACE_NUMS[sgnum²ᴰ]
subts³ᴰ = subduction_tables(sgnum³ᴰ; timereversal)
subts²ᴰ = BandGraphs.SubductionTable{2}[]
cts²ᴰ   = BandGraphs.Connection{2}[]
for (i, t) in enumerate(subts³ᴰ)
    if sgnum²ᴰ == 2
        xy = [3,1]
        z = 2
    elseif sgnum²ᴰ in (3, 5)
        xy = [2,1]
        z = 3
    elseif sgnum²ᴰ == 4
        xy = [2, 3]
        z = 1
    else
        xy = [1,2]
        z = 3
    end

    if !iszero(t.c.kᴳ.kv.cnst[z]) || !iszero(t.c.kᴴ.kv.cnst[z]) #|| !iszero(t.c.kᴴ.kv.free[z,:])
        continue # the connection involves k₃≠0
    else
        if !iszero(t.c.kᴴ.kv.free[:,z])
            @warn "Unexpectedly used γ as free parameter: may not be handled as intended..."
        end
    end
    c = BandGraphs.Connection{2}(
            BandGraphs.LabeledKVec{2}(t.c.kᴳ.label, KVec(t.c.kᴳ.kv.cnst[xy])), 
            BandGraphs.LabeledKVec{2}(t.c.kᴴ.label, KVec(t.c.kᴴ.kv.cnst[xy], t.c.kᴴ.kv.free[xy, xy]))
            )
    if sgnum²ᴰ == 6
        swapidxs = [1,2,4,3] # there is an awful Γ₃↔Γ₄ notational swap from 3D to 2D in p2mm
    else
        swapidxs = 1:size(t.table, 1)
    end
    table = t.table[swapidxs, :] 
    push!(subts²ᴰ, 
          BandGraphs.SubductionTable{2}(sgnum²ᴰ, c, t.irlabsᴳ, t.irlabsᴴ, table, t.monodromy)
    )
    push!(cts²ᴰ, c)
end
cts²ᴰ

## --------------------------------------------------------------------------------------- #
# Visualize "interesting" (i.e., more than 1 band) Hilbert vectors

sb, brs = compatibility_basis(sgnum²ᴰ, 2; timereversal)
lgirsd = lgirreps(sgnum²ᴰ, Val(2))
timereversal && realify!(lgirsd)

GLMakie.closeall()
for (j, _n) in enumerate(sb)
    println(j)
    _n[end] == 1 && continue # nothing interesting to draw for one-band configurations

    n = SymmetryVector(_n, brs.irlabs, lgirsd)
    bandg = build_subgraphs(n, subts²ᴰ, lgirsd)

    subgraphs_ps = permute_subgraphs(bandg.subgraphs);
    bandgp = BandGraphs.BandGraphPermutations(bandg.partitions, subgraphs_ps);
    BandGraphs.permutation_info(bandgp)
    length(bandgp) > 10 && @info("Be warned: many permutations ($(length(bandgp)))")

    xys = nothing
    maxplot = 4
    fs = Vector{Figure}(undef, maxplot)
    for (n, bandg′) in enumerate(bandgp)
        n > maxplot && break
        faxp, plot_data = plot_flattened_bandgraph(bandg′; xys=xys)
        fs[n], ax, p = faxp
        ax.title = "Basis vector j = $j"
        display(GLMakie.Screen(), fs[n])
        xys = (plot_data.xs, plot_data.ys)
    end
end