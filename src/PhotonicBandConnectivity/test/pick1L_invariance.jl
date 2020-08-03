if !isdefined(Main, :(PhotonicBandConnectivity))
    includet((@__DIR__)*"/../src/PhotonicBandConnectivity.jl") # TODO: broken for execution in vscode?
    using Main.PhotonicBandConnectivity
end

using Test
#= 
This is a test script to verify that it the choice of 1L mode doesn't matter (so long that
it is among the valid choices that actually match the symmetry constraints) and doesn't 
influence the eventual transverse symmetry representation (of course, it _does_ influence
the 1L+2T representation, but not 2T = [1L+2T]-[1L])
=#

# unshuffled 1L pick (first minimal)
sgnums = 1:230
data = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=true);
cⁱss   = getindex.(data, 1) # coefficients of expansions
νᵀs    = getindex.(data, 2) # fillings for tranverse branch
sbs    = getindex.(data, 3) # symmetry bases
idx¹ᴸs = getindex.(data, 4) # index for chosen 1L branch

# shuffled 1L pick (next element after first minimal)
data′ = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=true, shuffle_1Lpick=true);
cⁱss′   = getindex.(data′, 1)
νᵀs′    = getindex.(data′, 2)
sbs′    = getindex.(data′, 3)
idx¹ᴸs′ = getindex.(data′, 4)

# test that we get the same symmetry bases and the same transverse fillings
@test νᵀs == νᵀs′
@test sbs == sbs′

# associated transverse symmetry vectors nᵀ = n - n¹ᴸ
function transverse_expansion(sb, cⁱ, idx¹ᴸ)
    n = sum_symbases(sb, cⁱ)
    if !isnothing(idx¹ᴸ) && idx¹ᴸ > 0
        n .-= sb[idx¹ᴸ]
    end
    return n
end
nsᵀs  = [unique!(sort(transverse_expansion.(Ref(sb), cⁱs, idx¹ᴸ))) for (sb, cⁱs, idx¹ᴸ) in zip(sbs,  cⁱss,  idx¹ᴸs )]
nsᵀs′ = [unique!(sort(transverse_expansion.(Ref(sb), cⁱs, idx¹ᴸ))) for (sb, cⁱs, idx¹ᴸ) in zip(sbs′, cⁱss′, idx¹ᴸs′)]

# test that we get the same symmetry vectors for T-branches, regardless of our n¹ᴸ pick
for (sgnum, nsᵀ, nsᵀ′) in zip(sgnums, nsᵀs, nsᵀs′)
    @test nsᵀ == nsᵀ′
end