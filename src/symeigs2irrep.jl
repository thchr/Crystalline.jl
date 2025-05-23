"""
    find_representation(
        symvals::AbstractVector{<:Number}, 
        lgirs::AbstractVector{<:AbstractIrrep},
        assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int);
        αβγ::Union{AbstractVector{<:Real},Nothing}=nothing,
        atol::Real=DEFAULT_ATOL,
        maxresnorm::Real=1e-3,
        verbose::Bool=false
    )                                          --> Union{Nothing, Vector{assert_return_T}}

From a vector (or vector of vectors) of symmetry eigenvalues `symvals` sampled along all the
operations of a group gᵢ, whose irreps are contained in `irs` (evaluated with optional free 
parameters `αβγ`), return the multiplicities of each irrep.

Optionally, the multiciplities' element type can be specified via the `assert_return_T`
argument (performing checked conversion; returns `nothing` if representation in 
`assert_return_T` is impossible). This can be useful if one suspects a particular band to 
transform like a fraction of an irrep (i.e., the specified symmetry data is incomplete).

If no valid set of multiplicities exist (i.e., is solvable, and has real-valued and
`assert_return_T` representible type), the sentinel value `nothing` is returned. Optional
debugging information can in this case be shown by setting `verbose=true`.

# Extended help
Effectively, this applies the projection operator P⁽ʲ⁾ of each irrep's character set
χ⁽ʲ⁾(gᵢ) (j = 1, ... , Nⁱʳʳ) to the symmetry data sᵢ ≡ `symvals`:

    P⁽ʲ⁾  ≡ (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*gᵢ         [characters χ⁽ʲ⁾(gᵢ), irrep dimension dⱼ]
    P⁽ʲ⁾s = (dⱼ/|g|) ∑ᵢ χ⁽ʲ⁾(gᵢ)*sᵢ = nⱼ,   [number of bands that transform like jth irrep]

returning the irrep multiplicities mⱼ ≡ nⱼ/dⱼ.
"""
function find_representation(
    symvals::AbstractVector{<:Number}, 
    lgirs::AbstractVector{<:AbstractIrrep},
    assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int;
    αβγ::Union{AbstractVector{<:Real},Nothing}=nothing,
    atol::Real=DEFAULT_ATOL,
    maxresnorm::Real=1e-3,
    verbose::Bool=false
)
    ct = characters(lgirs, αβγ)
    χs = matrix(ct) # character table as matrix (irreps-as-columns & operations-as-rows)

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
    residual = χs*ms - symvals
    if norm(residual) > maxresnorm
        if verbose
            @info """computed multiplicities have a residual exceeding `maxresnorm`, 
                     suggesting that there is no solution to ``χm = s``; typically, this
                     means that `symdata` is "incomplete" (in the sense that it e.g., 
                     only contains "half" a degeneracy); `nothing` returned as sentinel"""
        end
        return nothing
    end

    # check that imaginary part is numerically zero and that all entries are representible
    # in type assert_return_T
    msℝ = real.(ms)
    if !isapprox(msℝ, ms, atol=atol)
        if verbose 
            @info """non-negligible imaginary components found in irrep multicity; returning
                     `nothing` as sentinel""" maximum(imag, ms)
        end
        return nothing
    end
    if assert_return_T <: Integer
        msT = round.(assert_return_T, msℝ)
        if isapprox(msT, msℝ, atol=atol)
            return msT
        else
            if verbose
                @info """multiplicities cannot be represented by the requested integer type,
                         within the specified absolute tolerance (typically, this means that
                         `symdata` is "incomplete" and needs to be "expanded" (i.e. add
                         bands)"""
            end
            return nothing
        end

    elseif assert_return_T <: AbstractFloat
        # TODO: not sure why we wouldn't want to check for `atol` distance from nearest
        #       integer here (returning `nothing` in the negative case) - should check if
        #       we ever actually rely on this
        return convert.(assert_return_T, msℝ)
    end
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Given a set of _band_-resolved symmetry eigenvalues `symeigs`, a set of irreps `irs`, and a
selection of band indices `bands` into `symeigs`, return the associated multiplicities of
corresponding to `irs`.

The `symeigs` argument is a vector of symmetry eigenvalues, or symmetry characters, for
individual bands. The `symeigs[n][i]`th entry gives the symmetry eigenvalue ``x_n(g_i)`` of
the `n`th band and the ``g_i =`` `group(irs)[i]` symmetry operation, with:

    ```math
    x_n(g_i) = \\langle \\psi_n | g_i | \\psi_n \\rangle
    ```
If unset, `bands` contains all bands in `symeigs`, i.e., is equal to `eachindex(symeigs)`.

## Keyword arguments
- `kws`: Additional keyword arguments passed to
  [`find_representation(::AbstractVector{<:Number})`](@ref) (e.g., `αβγ`, `atol`,
  `maxresnorm`).
"""
function find_representation(
    symeigs::AbstractVector{<:AbstractVector{<:Number}},
    irs::AbstractVector{<:AbstractIrrep},
    bands::AbstractVector{Int} = eachindex(symeigs),
    assert_return_T::Type{<:Union{Integer, AbstractFloat}}=Int;
    kws...
)                                                                           
    # return multiplicities over provided irreps (and across `bands`)
    symeigs_bands = sum(@view symeigs[bands])
    return find_representation(symeigs_bands, irs, assert_return_T; kws...)
end