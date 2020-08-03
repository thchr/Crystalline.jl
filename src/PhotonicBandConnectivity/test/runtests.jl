if !isdefined(Main, :(PhotonicBandConnectivity))
    include("../src/PhotonicBandConnectivity.jl")
    using Main.PhotonicBandConnectivity
end

Ms = getindex.(
        minimal_expansion_of_zero_freq_bands.(1:230, timereversal=true, verbose=false),
     2)

# Compare with Watanabe & Lu
Base.ndigits(::Nothing) = 1 # pirate haaaack
include("scripts/watanabelu_results.jl") # loads Watanabe & Lu data (in `Msᵂᴸ`)
Q = [[sg, M, Mbound] for (sg, M, Mbound) ∈ zip(1:230, Ms, getindex.(Msᵂᴸ, 2))]
Q′ = filter(x-> x[2]!==nothing, Q) # filter out those sgs that are not currently implemented (i.e. allow only regular 2T)
issues = map(x->x[2]===nothing ? "─" : (x[2]≥(x[3]) ? " " : "!"), Q)
differences = map(x->x[2]===nothing ? "─" : (x[2]==(x[3]) ? " " : "≠"), Q)

foreach(vcat.(Q, issues, differences)) do x
    println("|", " "^(4-ndigits(x[1])), x[1], " |", " "^(3-ndigits(x[2])),  # SG no.
    x[2] === nothing ? "─" : x[2], " | ",  # our M predictions
    x[3] == 2 ? "=" : "≥", x[3], " | ",    # M-bounds from Watanabe & Lu
    x[4], " | ",                           # bound violations
    x[5], " |"                             # differences from bound?
    )
end