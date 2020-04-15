"""
    find_representation(symvals::AbstractVector{Number}, 
                        lgirs::AbstractVector{<:AbstractIrrep},
                        αβγ::Union{AbstractVector{<:Real},Nothing}=nothing)
                                                                            --> Vector{Int}

From a vector (or vector of vectors) of symmetry eigenvalues `symvals` sampled along all the
operations of a group gᵢ, whose irreps are contained in `irs`, return the multiplicities of
each irrep.

Effectively, this applies the projection operator P⁽ʲ⁾ of each irrep (j = 1 , ..., Nⁱʳʳ) to
the symmetry data sᵢ ≡ `symvals`:
    P⁽ʲ⁾  ≡ (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*gᵢ         [characters χ⁽ʲ⁾(gᵢ), irrep dimension dⱼ]
    P⁽ʲ⁾s = (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*sᵢ = nⱼ,   [number of bands that transform like jth irrep]
returning the irrep multiplicities mⱼ ≡ nⱼ/dⱼ.
"""
function find_representation(symvals::AbstractVector{<:Number}, 
                             lgirs::AbstractVector{<:AbstractIrrep},
                             αβγ::Union{AbstractVector{<:Real},Nothing}=nothing)
    ct = CharacterTable(lgirs, αβγ)
    χs = ct.chartable # character table as matrix (irreps-as-columns & operations-as-rows)

    # METHOD/THEORY: From projection operators, we have
    #   (dⱼ/|g|) ∑ᵢ χᵢⱼ*sᵢ = nⱼ   ∀j                                                    (1)
    # where χᵢⱼ≡χ⁽ʲ⁾(gᵢ) are the characters of the ith operation gᵢ of the jth irrep, sᵢ are 
    # the symmetry eigenvalues (possibly, a sum!)¹ of gᵢ, and nⱼ is the number of bands that
    # transform like the jth irrep; finally, dⱼ is the jth irrep's dimension, and |g| is the
    # number of operations in the group. What we want here is the multiplicity mⱼ, i.e. the
    # number of times each irrep features in the band selection. mⱼ and nⱼ are related by 
    # mⱼ ≡ nⱼ/dⱼ and are both integers since each irrep must feature a whole number of
    # times. Thus, we have
    #   |g|⁻¹ ∑ᵢ χᵢⱼ*sᵢ = mⱼ   ∀j                                                      (2a)
    #   |g|⁻¹ χ†s = m          (matrix notation)                                       (2b)
    # Adrian suggested an equivalent approach, namely that ∑ⱼ χᵢⱼmⱼ = sᵢ ∀i, i.e. that 
    #   χm = s,                                                                         (3)
    # such that s can be decomposed in a basis of irrep characters. This is clearly an
    # appealing perspective. It is easy to see that (2) and (3) are equivalent; but only if
    # we start from (3). Specifically, from the orthogonality of characters, we have
    #   ∑ᵢ χʲ(gᵢ)*χᵏ(gᵢ) = δⱼₖ|g|   ⇔   χ†χ = |g|I                                      (4)
    # Thus, if we multiply (3) by χ†, we immediately obtain (2b) by using (4).
    # To go the other way, i.e. from (2) to (3), is harder since it is *not* true in general
    # that χχ† ∝ I. Instead, the situation is that χ/|g| is the *pseudo*inverse of χ†, in
    # the sense χ†(χ/|g|)χ† = χ†. This could almost certainly be exploited to go from (2b)
    # to (3), but it is not really important.
    #
    # ¹⁾ We can choose to work with a sum of symmetry eigenvalues of multiple bands - to
    #    obtain the irrep multiplicities of a band multiplet - because (1) is a linear 
    #    one-to-one relation.

    inv_g_order = inv(order(first(lgirs)))
    ms = inv_g_order .* (χs' * symvals) # irrep multiplicities obtained from (2b)

    #ms = χs\symvals # Adrian's approach

    # check that imaginary part is numerically zero and that all entries are integer-valued
    msℝ = real.(ms)
    msℕ = round.(Int, msℝ)
    @assert isapprox(msℝ, ms,  atol=DEFAULT_ATOL)
    @assert isapprox(msℕ, msℝ, atol=DEFAULT_ATOL)

    # TODO: Return something more meaningful than an assertion error if ms is not real or
    #       integer; could just mean that the symmetry data fed to the method is incomplete
    #       and needs to be "expanded" (i.e. add bands). E.g., return nothing as a sentinel?

    return msℕ
end
function find_representation(multiplet_symvals::AbstractVector{<:AbstractVector{<:Number}}, 
                             lgirs::AbstractVector{<:LGIrrep}, 
                             αβγ::Union{AbstractVector{<:Real},Nothing}=nothing)
    # If multiple symmetry values are given, as a vector of vectors, we take their sum over
    # the band-multiplet to determine irrep of the overall multiplet (see discussion in main
    # method)
    return find_representation(sum.(multiplet_symvals), lgirs, αβγ) 
end