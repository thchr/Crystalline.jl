using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)
using Crystalline
using BandGraphs
using SymmetryBases
using GraphMakie
using GLMakie
using Crystalline: irdim
using ProgressMeter

# --------------------------------------------------------------------------------------- #

D = 3
timereversal = true
criterion = (lgir) -> isspecial(lgir) && irdim(lgir) == 4
separable_degree = 2
check_sgnums = 1:230

contenders = Dict{Int, Vector{SymmetryVector{D}}}()
for sgnum in check_sgnums
    #sgnum in (47, 123) && continue # skip to be fast
    lgirsd = lgirreps(sgnum, D)
    timereversal && realify!(lgirsd)
    
    # check if any irreps fulfil `criterion` before computing Hilbert basis (can be costly!) 
    any(any(criterion, lgir for lgir in lgirs) for lgirs in values(lgirsd)) || continue
    print("#$sgnum, ")

    sb, brs = compatibility_basis(sgnum, D; timereversal)
    μs = fillings(sb)
    for i in eachindex(sb)
        μs[i] == 1 && continue # nothing interesting in 1-band cases
        n = SymmetryVector(sb[i], brs.irlabs, lgirsd)
        if has_vertex(criterion, n)
            push!(get!(contenders, sgnum, SymmetryVector{D}[]), n)
        end
    end
end
println("\n")

## --------------------------------------------------------------------------------------- #
# update-printing utilities

function formatted_time(t_s)
    if t_s < 60
        return string(round(Int, t_s)) * " s"
    else # minutes
        return (string(round(Int, div(t_s, 60))) * " min " * 
                string(floor(Int, rem(t_s, 60))) * " s")
    end
end

function print_status(
        io::IO, t₀, N, progress, separable_count, inseparable_count, infeasible_count;
        clear_line=true)
    # print time since start (t₀)
    printstyled(io, 
        clear_line ? "\r\e[K" : "", 
        "   ", formatted_time(time()-t₀), ": "; color=:light_black)
    
    # print iteration progress
    printstyled(io, progress, "/", N, " (", round(Int, progress/N*100), "%)";
                color=:light_black)
    
    # print number of separable, inseparable, and infeasible configurations   
    print(io, " [")
    separable_count > 0 && printstyled(io, separable_count; color=:green)
    separable_count > 0 && (inseparable_count > 0 || infeasible_count > 0) && print(io, " ")
    inseparable_count > 0 && printstyled(io, inseparable_count; color=:yellow)
    inseparable_count > 0 && infeasible_count > 0 && print(io, " ")
    infeasible_count > 0 && printstyled(io, infeasible_count; color=:red)
    print(io, "]")
end

function print_status(t₀, N, progress, separable_count, inseparable_count, infeasible_count)
    print_status(stdout, t₀, N, progress, separable_count, inseparable_count, infeasible_count)
end

## --------------------------------------------------------------------------------------- #
T_result = Tuple{BandGraph{D}, BandGraph{D}, LGIrrep{D}}
separable = Dict{Int, Vector{Tuple{Vector{LGIrrep{D}}, SymmetryVector{D}}}}()
separable_details = Dict{Int, Vector{Tuple{SymmetryVector{D}, Vector{T_result}}}}()
inseparable = Dict{Int, Vector{Tuple{SeparabilityState, SymmetryVector{D}}}}()
infeasible = Dict{Int, Vector{Tuple{SeparabilityState, SymmetryVector{D}}}}()

# using Profile
# Profile.init(; n = UInt(100_000_000), delay=0.01)

sgnums = sort!(collect(keys(contenders)))
for sgnum in sgnums
    sgnum == 130 && continue
    #sgnum ∈ (200, 205) && continue # stuck
    #sgnum > 199 && continue # slow
    #sgnum > 197 && continue
    ns = contenders[sgnum]
    subts = subduction_tables(sgnum, D; timereversal)
    lgirsd = timereversal ? realify!(lgirreps(sgnum, D)) : lgirreps(sgnum, D)
    vs = filter(criterion, collect(Iterators.flatten(values(lgirsd))))
    print("#$sgnum ($(length(ns)) symmetry vectors) [$(join(vs, ", "))]:\n")

    N_ns = length(ns)
    accumulate_results = Vector{Tuple{SeparabilityState, Union{Vector{T_result}, Nothing}}}(undef, N_ns)
    infeasible_count = Threads.Atomic{Int}(0)
    feasible_count = Threads.Atomic{Int}(0)
    separable_count = Threads.Atomic{Int}(0)
    progress_count = Threads.Atomic{Int}(0)

    t₀ = time()
    Threads.@threads for i in eachindex(ns)
        n = ns[i]
        sep, bandg_splits = findall_separable_vertices(
                                criterion, n, subts, lgirsd; 
                                max_permutations = 1e6,
                                separable_degree = separable_degree,
                                max_subblock_permutations = 1e6)
        accumulate_results[i] = (sep, bandg_splits)

        # print progress
        print_status(t₀, N_ns, progress_count[], separable_count[], feasible_count[], infeasible_count[])
        
        # update progress counters
        if is_separable(sep)
            Threads.atomic_add!(separable_count, 1)  # atomic `separable_count += 1`
        elseif is_feasible(sep)
            Threads.atomic_add!(feasible_count, 1)   # atomic `feasible_count += 1`
        else
            Threads.atomic_add!(infeasible_count, 1) # atomic `infeasible_count += 1`
        end
        Threads.atomic_add!(progress_count, 1) # atomic `progress_count += 1`
    end
    print("\r\e[K") # clear last progress line

    # separate out the results in various bins
    for ((sep, bandg_splits), n) in zip(accumulate_results, ns)
        if is_separable(sep)
            lgir = getindex.(bandg_splits, 3)
            push!(get!(separable, sgnum, SymmetryVector{D}[]), (lgir, n))
            push!(get!(separable_details, sgnum, Tuple{SymmetryVector{D}, T_result}[]), 
                (n, bandg_splits))
        elseif is_feasible(sep)
            push!(get!(inseparable, sgnum, Tuple{SeparabilityState, SymmetryVector{D}}[]), 
                (sep, n))
        else # not feasible; calculation failed
            push!(get!(infeasible, sgnum, Tuple{SeparabilityState, SymmetryVector{D}}[]), 
                (sep, n))
        end
    end
    
    if haskey(separable, sgnum)
        printstyled("   $(length(separable[sgnum])) separable configurations\n"; color=:green)
    end
    if haskey(infeasible, sgnum)
        min_perms, max_perms = extrema(something.(getfield.(getindex.(infeasible[sgnum], 1), Ref(:permutations))))
        printstyled("   $(length(infeasible[sgnum])) infeasible configurations ($min_perms to $max_perms permutations)\n"; color=:red)
    end
    if haskey(inseparable, sgnum)
        printstyled("   $(length(inseparable[sgnum])) inseparable configurations\n"; color=:yellow)
    end
    printstyled("   elapsed time: ", formatted_time(time() - t₀) * "\n"; color=:light_black)
end
separable

## --------------------------------------------------------------------------------------- #
#=
using KdotP
degeneracies = Dict{Int, Vector{KdotP.HamiltonianExpansion{D}}}()
for (sgnum, info) in separable
    lgirs = unique(reduce(vcat, getindex.(info, 1)))
    degeneracies[sgnum] = kdotp.(lgirs)
end
=#
## --------------------------------------------------------------------------------------- #
