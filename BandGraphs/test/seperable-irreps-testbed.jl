using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using Graphs
using SymmetryBases
using GraphMakie
using GLMakie
using BandGraphs: findall_separable_vertices
using Crystalline: irdim
using KdotP
using Dictionaries

## --------------------------------------------------------------------------------------- #
#logpath = joinpath((@__DIR__), "log.txt")
#redirect_stdio(; stdout=logpath, stderr=logpath, stdin=devnull) do
## --------------------------------------------------------------------------------------- #

D = 3
timereversal = true
criterion = lgir -> isspecial(lgir) && irdim(lgir) == 3 # KdotP.isweyl(lgir; timereversal)
separable_degree = 2

sgnums = 1:MAX_SGNUM[D]
contendersd = Dictionary{Int, Vector{SymVector{D}}}()
lgirsdd = Dictionary{Int, Dict{String, Crystalline.IrrepCollection{LGIrrep{D}}}}()
subtsd = Dictionary{Int, Vector{SubductionTable{D}}}()
for sgnum in sgnums
    lgirsd = lgirreps(sgnum, D)
    timereversal && realify!(lgirsd)
    has_relevant_lgirs = any(lgirsd) do (klab, lgirs)
        any(criterion, lgirs)
    end
    if has_relevant_lgirs
        insert!(contendersd, sgnum, SymVector{D}[])
        insert!(lgirsdd,     sgnum, lgirsd)
        insert!(subtsd,      sgnum, subduction_tables(sgnum, D; timereversal))
    else
        continue
    end
    printstyled("#", sgnum; bold=true)

    sb, brs = compatibility_basis(sgnum, D; timereversal)
    μs = fillings(sb)
    for i in eachindex(sb)
        μs[i] == 1 && continue # nothing interesting in 1-band cases

        n = SymVector(sb[i], brs.irlabs, lgirsd)
        for (lgirs, mults) in zip(n.lgirsv, n.mults)
            for (lgir, m) in zip(lgirs, mults)
                if m ≠ 0 && criterion(lgir)
                    push!(contendersd[sgnum], n)
                    @goto SYMVECTOR_DONE
                end
            end
        end
        @label SYMVECTOR_DONE
    end
    printstyled(" (", length(contendersd[sgnum]), " contenders)\n"; color=:light_black)
    flush(stdout)
    flush(stderr)
end


## --------------------------------------------------------------------------------------- #
separable_summary = Dictionary{Int, Vector{Tuple{Vector{LGIrrep{D}}, SymVector{D}}}}()
separable_details = Dictionary{Int, Vector{Tuple{SymVector{D}, Vector{Tuple{BandGraph{D}, BandGraph{D}, LGIrrep{D}}}}}}()
inseparable = Dictionary{Int, Vector{SymVector{D}}}()

done = Set{Int}()
too_hard = Set{Int}()
for (sgnum, ns) in pairs(contendersd)
    #sgnum == 84 && continue
    #sgnum ∈ (178, 181, 179, 180, 214, 199) && continue # untractably many permutations
    #sgnum < 211 && continue
    #sgnum == 210 && continue # segfault!?
    #sgnum ≠ 214 && continue # know there has to be examples here (gyroid)! segfault!?
    #sgnum ≠ 230 && continue # know there has to be examples here (double gyroid)! segfault!?
    printstyled("#", sgnum, " "; bold=true)
    printstyled("(", length(ns), " contenders)\n", color=:light_black)
    for (i, n) in enumerate(ns)
        #i ≤ 324 && continue
        println(" i=", i, ": ", n)
        subts = subtsd[sgnum]
        lgirsd = lgirsdd[sgnum]
        
        try
            tmp = findall_separable_vertices(
                        criterion, n, subts, lgirsd; 
                        separable_degree = separable_degree,
                        max_permutations = 1e5)
            if ismissing(tmp)
                push!(too_hard, sgnum)
                #printstyled("  Skipping any remaining symmetry vectors in this space group!\n"; 
                #    color=:red)
                continue
                #break # too many permutations
            end
            has_split, bandg_splits = tmp
            if has_split
                lgir = getindex.(bandg_splits, 3)
                push!(get!(separable_summary, sgnum, Vector{SymVector{D}}()), (lgir, n))
                push!(get!(separable_details, sgnum,
                        Vector{Tuple{SymVector{D}, typeof(bandg_splits)}}()), 
                        (n, bandg_splits))
                printstyled("   HAS SEPARABLE CONFIGURATIONS!!!\n", color=:green)    
            else
                push!(get!(inseparable, sgnum, Vector{SymVector{D}}()), n)
            end
        catch e
            printstyled("  ", e, "\n",
                       #"  Skipping any remaining symmetry vectors in this space group!\n"; 
                       color=:red)
            push!(too_hard, sgnum)
            continue
            #break # too many permutations
        end
    end
    if sgnum ∉ too_hard
        push!(done, sgnum)
    end
    flush(stdout)
    flush(stderr)
end
separable_summary
#end
## --------------------------------------------------------------------------------------- #

using KdotP
degeneracies = Dictionary{Int, Vector{KdotP.HamiltonianExpansion{D}}}()
for (sgnum, info) in pairs(separable_summary)
    lgirs = unique(reduce(vcat, getindex.(info, 1)))
    inset!(degeneracies, sgnum, kdotp.(lgirs))
end

## --------------------------------------------------------------------------------------- #
