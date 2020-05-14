"""
    find_representation(symvals::AbstractVector{Number}, 
                        lgirs::AbstractVector{<:AbstractIrrep},
                        αβγ::Union{AbstractVector{<:Real},Nothing}=nothing,
                        assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int))

                        --> Vector{assert_return_T}

From a vector (or vector of vectors) of symmetry eigenvalues `symvals` sampled along all the
operations of a group gᵢ, whose irreps are contained in `irs` (evaluated with optional free 
parameters `αβγ`), return the multiplicities of each irrep.

Optionally, the multiciplities' element type can be specified via the `assert_return_T`
argument (performing checked conversion; returns `nothing` if representation in 
`assert_return_T` is impossible). This can be useful if one suspects a particular band to 
transform like a fraction of an irrep (i.e., the specified symmetry data is incomplete).

# Extended help
Effectively, this applies the projection operator P⁽ʲ⁾ of each irrep's character set
χ⁽ʲ⁾(gᵢ) (j = 1, ... , Nⁱʳʳ) to the symmetry data sᵢ ≡ `symvals`:

    P⁽ʲ⁾  ≡ (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*gᵢ         [characters χ⁽ʲ⁾(gᵢ), irrep dimension dⱼ]
    P⁽ʲ⁾s = (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*sᵢ = nⱼ,   [number of bands that transform like jth irrep]

returning the irrep multiplicities mⱼ ≡ nⱼ/dⱼ.
"""
function find_representation(symvals::AbstractVector{<:Number}, 
                             lgirs::AbstractVector{<:AbstractIrrep},
                             αβγ::Union{AbstractVector{<:Real},Nothing}=nothing,
                             assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int)
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
    # Unfortunately, (1) doesn't seem to hold if we work with "realified" co-reps (i.e. 
    # "doubled" complex or pseudoreal irreps). In that case, we do _not_² have χ†χ = |g|I. 
    # TODO: I'm not sure why this ruins (2) and not (3), but that seems to be the case. 
    #       As such, we just use (3), even if it is slightly slower.
    #
    # ¹⁾ We can choose to work with a sum of symmetry eigenvalues of multiple bands - to
    #    obtain the irrep multiplicities of a band multiplet - because (1) is a linear 
    #    one-to-one relation.
    # ²⁾ Instead, χ†χ is still a diagonal matrix, but the diagonal entries corresponding to
    #    doubled irreps are 2|g| instead of |g| (from testing).

    #inv_g_order = inv(order(first(lgirs)))                         # [not used; see above]
    #ms = inv_g_order .* (χs' * symvals) # irrep multiplicities obtained from (2b)
    
    ms = χs\symvals # Adrian's approach

    # check that imaginary part is numerically zero and that all entries are representible
    # in type assert_return_T
    msℝ = real.(ms)
    @assert isapprox(msℝ, ms,  atol=DEFAULT_ATOL)
    if assert_return_T <: Integer
        msT = round.(assert_return_T, msℝ)
        if isapprox(msT, msℝ, atol=DEFAULT_ATOL)
            return msT
        else
            # if ms cannot be represented by the the requested integer type, we return a
            # sentinel value of `nothing` to indicate this (could e.g. mean that the input
            # symmetry data is "incomplete" and needs to be "expanded" (i.e. add # bands))
            return nothing
        end

    elseif assert_return_T <: AbstractFloat
        return convert.(assert_return_T, msℝ)
    end

end
function find_representation(multiplet_symvals::AbstractVector{<:AbstractVector{<:Number}}, 
                             lgirs::AbstractVector{<:LGIrrep}, 
                             αβγ::Union{AbstractVector{<:Real},Nothing}=nothing,
                             assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int)
    # If multiple symmetry values are given, as a vector of vectors, we take their sum over
    # the band-multiplet to determine irrep of the overall multiplet (see discussion in main
    # method)
    return find_representation(sum.(multiplet_symvals), lgirs, αβγ, assert_return_T) 
end