using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)

using Crystalline
using BandGraphs
using JLD2
using Test

connectionsd = load_object(joinpath((@__DIR__), "..", "data/connections/3d/connections.jld2"))
subductionsd = load_object(joinpath((@__DIR__), "..", "data/connections/3d/subductions.jld2"))

subductionsd′ = Dict{Int, Vector{SubductionTable{3}}}()
failures = Dict{Int, Any}()
errors   = Dict{Int, Any}()
for sgnum in 1:230
    subts = subductionsd[sgnum]
    # ERROR: nonzero reciprocal vector between nonspecial kv and kv′ requires explicit handling: not yet implemented
    sgnum ∈ (26, 27, 28, 29, 30, 31, 32, 33, 34, 39, 40, 41, 100, 101, 102, 103, 104, 106, 116, 117, 118, 226) && continue
    # ERROR: DomainError with Provided irreps are not H<G subgroups:
    #sgnum ∈ (43, 70, 199, 203, 206, 210, 214, 220, 227, 228, 230) && continue
    # ERROR: failed to prove compatibility of X and Y
    sgnum in (121, 122, 203) && continue

    sg = spacegroup(sgnum, Val(3))
    lgirsd = lgirreps(sgnum, Val(3))
    subts′ = SubductionTable{3}[]
    println(sgnum)
    for subt in subts
        c = subt.c
        println("   ", c)
        subt′ = try
            SubductionTable(c, sg, lgirsd)
        catch
            haskey(errors, sgnum) || (errors[sgnum] = [])
            push!(errors[sgnum], c)
            printstyled("error: ", color=:red)
            println("#", sgnum, " & ", c)
            continue
        end
        permuteᴳ_idxs′ = [findfirst(==(irlab), subt′.irlabsᴳ) for irlab in subt.irlabsᴳ]
        permuteᴴ_idxs′ = [findfirst(==(irlab), subt′.irlabsᴴ) for irlab in subt.irlabsᴴ]
        if subt.table != subt′.table[permuteᴳ_idxs′, permuteᴴ_idxs′]
            display(subt)
            display(subt′)
            
            haskey(failures, sgnum) || (failures[sgnum] = [])
            push!(failures[sgnum], (; subt=subt, subt′=subt′))
            printstyled("failure: ", color=:yellow)
            println("#", sgnum, " & ", c)
        end
        push!(subts′, subt)
    end
    #subductionsd′[sgnum] = subts
end