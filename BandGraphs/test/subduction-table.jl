using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)

using Crystalline
using BandGraphs
using JLD2
using Test

timereversal = true
trtag = timereversal ? "-tr" : ""
dir = joinpath((@__DIR__), "..", "data/connections/3d")
connectionsd = load_object(joinpath(dir, "connections$(trtag).jld2"))
subductionsd = load_object(joinpath(dir, "subductions$(trtag).jld2"))

subductionsd′ = Dict{Int, Vector{SubductionTable{3}}}()
failures = Dict{Int, Any}()
errors   = Dict{Int, Any}()
for sgnum in 1:230
    subts = subductionsd[sgnum]

    sg = spacegroup(sgnum, Val(3))
    lgirsd = lgirreps(sgnum, Val(3))
    timereversal && realify!(lgirsd)
    
    subts′ = SubductionTable{3}[]
    println(sgnum)
    for subt in subts
        c = subt.c
        println("   ", c)
        subt′ = try
            SubductionTable(c, sg, lgirsd)
        catch e
            haskey(errors, sgnum) || (errors[sgnum] = [])
            push!(errors[sgnum], (; c=c, error=e))
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