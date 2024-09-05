using Bravais
using LinearAlgebra: norm

# Example from Krivy & Gruber [1](https://doi.org/10.1107/S0567739476000636):
# NB: the example given by Krivy & Gruber first gives values for basis lengths & angles,
# i.e., a,b,c,α,β,γ - and then give associated values for the niggli parameters A, B, C, ξ,
# η, ζ. However, the values they give for the Niggli parameters are given with 0 significant
# digits - and it appears they actually use these zero-significant-digit values in the
# actual calculation they show afterwards. So here, we start from the provided values for 
# A, b, C, ξ, η, ζ, and then calculate associated values for a, b, c, α, β, γ; the latter
# then deviates from those given in [1] due to their rounding, but are consistent with the
# actual calculation that they do (although they appear to have done the calculation with
# integers rather than floating points; and their evaluation of the associated angles is
# ultimately also very imprecise; see large `atol` values below)
# Starting from the a, b, c, α, β, γ values given in [1], rather than the Niggli parameters,
# actually lead to a very different Niggli reduced cell, highlighting that the approach is
# not terribly robust to imprecision / measurement errors.

A, B, C, ξ, η, ζ = 9, 27, 4, -5, -4, -22
a, b, c = sqrt(A), sqrt(B), sqrt(C)
α, β, γ = acos(ξ/(2b*c)), acos(η/(2c*a)), acos(ζ/(2a*b))
Rs = crystal(a, b, c, α, β, γ)
Rs′, P = niggli_reduce(Rs; rtol = 1e-5, max_iterations=100)
abc = norm.(Rs′)
αβγ = [(Bravais.angles(Rs′) .* (180/π))...]

@test abc ≈ [2.0,3.0,3.0]      atol=1e-3 # norm accuracy from [1] is low
@test αβγ ≈ [60.0,75.31,70.32] atol=3e-1 # angle accuracy from [1] is _very_ low
@test Rs′ ≈ transform(Rs, P)