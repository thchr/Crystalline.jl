using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using MetaGraphsNext
using GLMakie
using BandGraphs: BandGraphPermutations
using BandGraphs: findall_separable_vertices
using Crystalline: irdim
## --------------------------------------------------------------------------------------- #

D = 2
timereversal = false
criterion = (lgir) -> irdim(lgir) == 2

contenders = Vector{Pair{Int, SymmetryVector{D}}}()
lgirsdd = Dict{Int, Dict{String, Collection{LGIrrep{D}}}}()
subtsd = Dict{Int, Vector{SubductionTable{D}}}()
for sgnum in 1:MAX_SGNUM[D]
    #sgnum ∈ (4) && continue # TODO: generalize to handle disconnected partitions, but #4 is not separable anyway
    lgirsd = lgirreps(sgnum, D)
    timereversal && realify!(lgirsd)

    sb, brs = compatibility_basis(sgnum, D; timereversal)
    μs = fillings(sb)
    had_contender = false
    for i in eachindex(sb)
        μs[i] == 1 && continue # nothing interesting in 1-band cases

        _n = sb[i]
        n = SymmetryVector(_n, brs.irlabs, lgirsd)
        for (lgirs, mults) in zip(irreps(n), multiplicities(n))
            for (lgir, m) in zip(lgirs, mults)
                if m ≠ 0 && criterion(lgir)
                    push!(contenders, sgnum => n)
                    had_contender = true
                    @goto SYMVECTOR_DONE
                end
            end
        end
        @label SYMVECTOR_DONE
    end
    if had_contender
        lgirsdd[sgnum] = lgirsd
        subtsd[sgnum] = subduction_tables(sgnum, D; timereversal)
    end
end


## --------------------------------------------------------------------------------------- #
separable_summary = Dict{Int, Vector{Tuple{Vector{LGIrrep{D}}, SymmetryVector{D}}}}()
separable_details = Dict{Int, Vector{Tuple{SymmetryVector{D}, Vector{Tuple{BandGraph{D}, BandGraph{D}, LGIrrep{D}}}}}}()
inseparable = Dict{Int, Vector{SymmetryVector{D}}}()

for (i, (sgnum, n)) in enumerate(contenders)
    subts = subtsd[sgnum]
    lgirsd = lgirsdd[sgnum]
    has_split, bandg_splits = findall_separable_vertices(criterion, n, subts, lgirsd)
    if has_split
        lgir = getindex.(bandg_splits, 3)
        push!(get!(separable_summary, sgnum, Vector{SymmetryVector{D}}()), (lgir, n))
        push!(get!(separable_details, sgnum, 
                Vector{Tuple{SymmetryVector{D}, typeof(bandg_splits)}}()), 
                (n, bandg_splits))      
    else
        push!(get!(inseparable, sgnum, Vector{SymmetryVector{D}}()), n)
    end
end
separable_summary

## --------------------------------------------------------------------------------------- #

using KdotP
degeneracies = Dict{Int, Vector{KdotP.HamiltonianExpansion{D}}}()
for (sgnum, info) in separable_summary
    lgirs = unique(reduce(vcat, getindex.(info, 1)))
    degeneracies[sgnum] = kdotp.(lgirs)
end

## --------------------------------------------------------------------------------------- #
