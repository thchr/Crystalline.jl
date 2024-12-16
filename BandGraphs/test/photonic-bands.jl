using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline, BandGraphs, PhotonicBandConnectivity, SymmetryBases
using GLMakie, GraphMakie
const B = BandGraphs

## --------------------------------------------------------------------------------------- #

criterion = (lgir) -> isspecial(lgir) && irdim(lgir) == 2 && !is_vrep(lgir)
separable_degree = 2

D = 3
timereversal = false


sepsd = Dict{Int,Vector{Any}}()
for sgnum in 1:230#1:MAX_SGNUM[D]
    sgnum ∈ (47, 83, 123) && (println("#", sgnum, " (SKIPPED)"); continue) # slow
    #sgnum ∈ (39, 41, 45, 76, 78, 108, 120) && continue # bugged for 2D irrep
    #sgnum ∈ (118, 139, 144, 145, 209) && (println("#", sgnum, " (SKIPPED)"); continue) # broken
    print("#$sgnum")
    lgirsd = timereversal ? realify!(lgirreps(sgnum, Val(D))) : lgirreps(sgnum, Val(D))
    
    brs = calc_bandreps(sgnum, Val(D); timereversal)
    ns = transverse_symmetry_vectors(brs, lgirsd)
    subts = subduction_tables(sgnum; timereversal)
    lgirsv = irreps(first(ns))
    vrep = something((lgirsv -> begin for lgirs in lgirsv
        i = findfirst(is_vrep, lgirs)
        isnothing(i) || return lgirs[i]
    end end)(lgirsv))
    BandGraphs.add_transverse_vrep!(subts, vrep)
    
    #=
    sb, brs = compatibility_basis(sgnum; timereversal)
    ns = SymmetryVector.(sb, Ref(brs.irlabs), Ref(lgirsd))
    subts = subduction_tables(sgnum; timereversal)
    =#

    println(" (", length(ns), " solutions):")
    for i in eachindex(ns)
        n = ns[i]
        #=
        sep, bandg_splits = findall_separable_vertices(
                                criterion, n, subts, lgirsd; 
                                max_permutations = 1e4,
                                separable_degree = separable_degree,
                                max_subblock_permutations = 1e4,
                                with_weyl_filter = true)

        if is_separable(sep) || !is_feasible(sep)
            println("    ", i, ": ", sep)
        end
        if is_separable(sep)
            #exfiltrate_transverse_vrep(n)
            push!(get!(sepsd, sgnum, Vector{Any}()), (sep, n, bandg_splits))
        end
        =#

        results = findall_separable_vertices_treecutting(
            criterion, n, subts, lgirsd; 
            max_permutations = 1e6,
            separable_degree = separable_degree,
            max_subblock_permutations = 1e7,
            with_weyl_filter = false,
            verbose = false)
        seps         = getindex.(results, 1)
        vs           = getindex.(results, 2)
        bandg_splits = getindex.(results, 3)

        i_print_str, i_shown = sprint(print, "  ", i, ": "), false
        sep_idxs = findall(is_separable, seps)
        if !isempty(sep_idxs)
            sep = first(seps[sep_idxs])
            i_shown || (i_shown = true; print(i_print_str))
            println(sep, " @ ", join(first.(vs[sep_idxs]), ", "))
            push!(get!(sepsd, sgnum, []), (sep, n, bandg_splits[sep_idxs]))
        end
        insep_idxs  = findall(s->is_feasible(s) && !is_separable(s), seps)
        if !isempty(insep_idxs)
            for i in insep_idxs
                i_shown ? print(" "^length(i_print_str)) : (i_shown=true; print(i_print_str))
                println(seps[i], " @ ", vs[i])
            end
        end
        infeas_idxs = findall(!is_feasible, seps)
        if !isempty(infeas_idxs)
            for i in infeas_idxs
                i_shown ? print(" "^length(i_print_str)) : (i_shown=true; print(i_print_str))
                println(seps[i], " @ ", vs[i])
            end
        end
    end
end

#occupation.(ns) .=> seps