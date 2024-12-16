struct MultiSetPermutationIndices
    f::Vector{Int}
end
Base.eltype(::Type{MultiSetPermutationIndices}) = Vector{Int}
Base.length(c::MultiSetPermutationIndices) = _multinomial(c.f)

# copy of https://github.com/JuliaMath/Combinatorics.jl/pull/170, until merged & public
function _multinomial(k)
    s = zero(eltype(k))
    result = one(eltype(k))
    @inbounds for i in k
        s += i
        result *= binomial(s, i)
    end
    result
end

"""
    multiset_permutation_indices(f :: Vector{Int})

Return a generator over all index permutations corresponding to the multiset permutations of
a multiset `f = {f[1]*1, f[2]*2, ...}`, i.e. of a multiset with `f[1]` elements of `1`,
`f[2]` elements of `2`, and so on.

More explicitly, `f[i]` indicates the frequency of the `i`th element of the underlying
multiset `a`, which can consequently be considered as the sequence or sorted multiset:
`a` = [`1` repeated `f[1]` times, `2` repeated `f[2]` times, ...].

The indices `inds` of the unique permutations of `a` are returned. Compared to the iterator
`multiset_permutations(a, length(a))` - which returns the actual permutations, not the
indices - from Combinatorics.jl, this means that
`a[ind] == collect(multiset_permutations(a), length(a))[ind]` for all `ind in inds`.

## Implementation
The implementation here is adapted from the implementation in `multiset_permutations` in
Combinatorics.jl, which in turn appears to follow the Algorithm L (p. 319) from Knuth's
TAOCP Vol. 4A.
"""
multiset_permutation_indices(f::Vector{Int}) = MultiSetPermutationIndices(f)

function canonical_multiset(f::Vector{Int})
    # build the sequence `a` = [`1` repeated `f[1]` times, `2` repeated `f[2]` times, ...].
    a = Vector{Int}(undef, sum(f))
    N = 1
    for (i, n) in enumerate(f)
        @inbounds a[N:N+n-1] .= i
        N += n
    end
    return a
end

function Base.iterate(
        p::MultiSetPermutationIndices, 
        composite_state = (a=canonical_multiset(p.f); (a, collect(eachindex(a))))
        )
    state, indices = composite_state
    if !isempty(state) && state[1] > length(state)
        return nothing
    end
    return nextpermutation!(state, indices)
end

function nextpermutation!(state::Vector{Int}, indices::Vector{Int})
    i = length(state) - 1
    indices_next = copy(indices)
    while i ≥ 1 && state[i] ≥ state[i+1]
        i -= 1
    end
    if i > 0
        j = length(state)
        while j > i && state[i] ≥ state[j]
            j -= 1
        end
        state[i], state[j] = state[j], state[i]
        indices_next[i], indices_next[j] = indices_next[j], indices_next[i]
        reverse!(state, i+1)
        reverse!(indices_next, i+1)
    else
        state[1] = length(state) + 1 # finished; no more permutations
    end

    return indices, (state, indices_next)
end