using Crystalline # reexports Bravais' API
using LinearAlgebra: norm
using Test

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

@testset "Niggli reduction" begin
A, B, C, ξ, η, ζ = 9, 27, 4, -5, -4, -22
a, b, c = sqrt(A), sqrt(B), sqrt(C)
α, β, γ = acos(ξ/(2b*c)), acos(η/(2c*a)), acos(ζ/(2a*b))
Rs = crystal(a, b, c, α, β, γ)
Rs′, P = nigglibasis(Rs; rtol = 1e-5, max_iterations=100)
abc = norm.(Rs′)
αβγ = [(Bravais.angles(Rs′) .* (180/π))...]

@test abc ≈ [2.0,3.0,3.0]      atol=1e-3 # norm accuracy from [1] is low
@test αβγ ≈ [60.0,75.31,70.32] atol=3e-1 # angle accuracy from [1] is _very_ low
@test Rs′ ≈ transform(Rs, P)

# the angles and lengths of the Niggli bases extracted from two equivalent lattice bases,
# mutually related by an element of SL(3, ℤ), must be the same (Niggli cell is unique)
P₁ = [1 2 3; 0 1 8; 0 0 1]
P₂ = [1 0 0; -3 1 0; -9 2 1]
P  = P₂ * P₁ # random an element of SL(3, ℤ) (i.e., `det(P) == 1`)
for sgnum in 1:MAX_SGNUM[3]
    Rs  = directbasis(sgnum)
    Rs′ = transform(Rs, P)
    niggli_Rs  = nigglibasis(Rs)[1]
    niggli_Rs′ = nigglibasis(Rs′)[1]
    abc  = norm.(niggli_Rs)
    abc′ = norm.(niggli_Rs′)
    αβγ  = [Bravais.angles(niggli_Rs)...]
    αβγ′ = [Bravais.angles(niggli_Rs′)...]
    @test abc ≈ abc′
    @test αβγ ≈ αβγ′
    @test sort(abc) ≈ abc # sorted by increasing norm
    ϵ_ϕ = 1e-10 # angle tolerance in comparisons below
    @test all(ϕ -> ϕ<π/2 + ϵ_ϕ, αβγ) || all(ϕ -> ϕ>π/2 - ϵ_ϕ, αβγ)
    @test isapprox(collect(Bravais.niggli_parameters(niggli_Rs)),  # equivalent to checking
                   collect(Bravais.niggli_parameters(niggli_Rs′))) # abc ≈ abc′ & αβγ ≈ αβγ′

    # Idempotency of Niggli-reduced Niggli-parameters
    niggli_niggli_Rs  = nigglibasis(niggli_Rs)[1]
    niggli_niggli_Rs′ = nigglibasis(niggli_Rs′)[1]
    @test isapprox(collect(Bravais.niggli_parameters(niggli_Rs)),
                   collect(Bravais.niggli_parameters(niggli_niggli_Rs)))
    @test isapprox(collect(Bravais.niggli_parameters(niggli_Rs)),
                   collect(Bravais.niggli_parameters(niggli_niggli_Rs′)))
    # the following tests are not guaranteed (i.e., the Niggli parameters (lengths & angles)
    # must be idempotent, but cannot guarantee exact choice of basis vectors `(...)Rs`)
    @test niggli_Rs ≈ niggli_niggli_Rs
    @test niggli_Rs′ ≈ niggli_niggli_Rs′
    #@test niggli_Rs ≈ niggli_Rs′         # <-- in particular, this does not hold generally
end

# Niggli reduction of reciprocal lattice
for sgnum in 1:MAX_SGNUM[3]
    Rs = directbasis(sgnum)
    Gs = reciprocalbasis(Rs)

    niggli_Rs, P_from_Rs = nigglibasis(Rs)
    niggli_Gs, P_from_Gs = nigglibasis(Gs)
    @test reciprocalbasis(niggli_Rs) ≈ niggli_Gs
    @test reciprocalbasis(niggli_Gs).vs ≈ niggli_Rs.vs # duality of Niggli reduction
    @test niggli_Gs ≈ transform(Gs, P_from_Gs)
    @test P_from_Rs ≈ P_from_Gs
end

# Niggli reduction of 2D lattice
P₁²ᴰ = [1 4; 0 1]
P₂²ᴰ = [1 0; -3 1]
P²ᴰ  = P₂²ᴰ * P₁²ᴰ # random an element of SL(2, ℤ)
for sgnum in 1:MAX_SGNUM[2]
    Rs = directbasis(sgnum, Val(2))
    niggli_Rs, niggli_P = nigglibasis(Rs)
    @test niggli_Rs ≈ transform(Rs, niggli_P)

    Rs′ = transform(Rs, P²ᴰ)
    niggli_Rs′ = nigglibasis(Rs′)[1]
    ab  = norm.(niggli_Rs)
    ab′ = norm.(niggli_Rs′)
    α  = only(Bravais.angles(niggli_Rs))
    α′ = only(Bravais.angles(niggli_Rs′))
    @test ab ≈ ab′
    @test α ≈ α′
    @test sort(ab) ≈ ab # sorted by increasing norm
    ϵ_ϕ = 1e-10 # angle tolerance in comparisons below
    @test α < π/2 + ϵ_ϕ || α > π/2 - ϵ_ϕ
end

end # @testset "Niggli reduction"