# Equivalent to (unchecked) `rand(Uniform(low, high))` from Distributions.jl, but copied
# in here to avoid taking on that dependency (long load time; ~7 s)
uniform_rand(low::Real, high::Real) = low + (high - low) * rand()

"""
    relrand(lims::NTuple{2,Real}) --> Float64

Computes a random number in the range specified by the two-element 
tuple `lims`. The random numbers are sampled from two uniform 
distributions, namely [`lims[1]`, 1] and [1, `lims[2]`], in such a
way as to ensure that the sampling is uniform over the joint  
interval [-1/`lims[1]`, -1] ∪ [1, `lims[2]`].

This is useful for ensuring an even sampling of numbers that are
either smaller or larger than unity. E.g. for `x = relrand((0.2,5.0))`,
`x` is equally probable to fall in inv(`x`)∈[1,5] or `x`∈[1,5].
"""
function relrand(lims::NTuple{2,<:Real})
    low, high = lims; invlow = inv(low)
    lowthres = (invlow - 1.0)/(invlow + high - 2.0)
    if rand() < lowthres && low < 1.0   # smaller than 1.0
        r = uniform_rand(low, 1.0)
    elseif high > 1.0                   # bigger than 1.0
        r = uniform_rand(1.0, high)
    else                                # default
        return uniform_rand(low, high)
    end
end


# === frequent error messages ===
# Note: there are some copied/overlapping definitions w/ Crystalline here
@noinline function _throw_invalidcntr(cntr::AbstractChar, D::Integer)
    if D == 3
        throw(DomainError(cntr,
            "centering abbreviation must be either P, I, F, R, A, or C in dimension 3"))
    elseif D == 2
        throw(DomainError(cntr,
            "centering abbreviation must be either p or c in dimension 2"))
    elseif D == 1
        throw(DomainError(cntr,
            "centering abbreviation must be p in dimension 1"))
    else
        _throw_invalid_dim(D)
    end
end
@noinline function _throw_invaliddim(D::Integer)
    throw(DomainError(D, "dimension must be 1, 2, or 3"))
end
@noinline function _throw_invalid_sgnum(sgnum::Integer, D::Integer)
    throw(DomainError(sgnum,
        "sgnum must be between 1 and $((2, 17, 230)[D]) in dimension $D"))
end

@inline function boundscheck_sgnum(sgnum::Integer, D::Integer)
    if D == 3 
        sgnum ∈ 1:230 || _throw_invalid_sgnum(sgnum, 3)
    elseif D == 2
        sgnum ∈ 1:17  || _throw_invalid_sgnum(sgnum, 2)
    elseif D == 1 && 
        sgnum ∈ 1:2   || _throw_invalid_sgnum(sgnum, 1)
    else
        _throw_invaliddim(D)
    end
end