module Jagged # TODO: Register as a separate package some day?

# ---------------------------------------------------------------------------------------- #

export JaggedVector
using Base: @propagate_inbounds

# ---------------------------------------------------------------------------------------- #

# ::: struct definition :::
const VectorView{T} = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}
"""
    JaggedVector{T}

A memory-efficient representation of a `Vector{Vector{T}}`, avoiding the overhead of
multiple inner allocations.

The iterants of `JaggedVector` are returned as views into an underlying shared `data`
vector, which is a flattened version of the nested vectors.

## Changing lengths of iterants

Due to `JaggedVector`'s memory-layout, the operation `jv[i] = v` that changes the iterant
`jv[i]` to a new vector `v` may have surprising performance characteristics: in particular,
if `length(v) ≠ length(jv[i])`, this is an O(`length(jv)`) operation, since the underlying
storage vector must be resized to allowed shuffling of adjacent iterants.
If this is a frequent occurence, consider using `Vector{Vector{T}}` instead.

Conversely, if `length(v) == length(jv[i])`, this remains an O(1) operation with
performance nearly identical to that of `Vector{Vector{T}}`.
"""
struct JaggedVector{T} <: AbstractVector{VectorView{T}}
    data :: Vector{T}
    offsets :: Vector{Int}
    # unsafe/unvalidated inner constructor
    function JaggedVector{T}(data::Vector{T}, offsets::Vector{Int}, ::Val{:unsafe}) where T
        return new{T}(data, offsets)
    end
end

# ::: safe/validating constructor :::
"""
    JaggedVector{T}(data::Vector{T}, offsets::Vector{Int})
    JaggedVector(data::Vector{T}, offsets::Vector{Int})

Given a set of `data` and associated `offsets` into `data`, return a `JaggedVector{T}`,
whose iterants return views into the `data` vector, as defined by `offsets`.
Construction is safe and validated.

## Structure of inputs
- `data`: a "flat" vector, representing the concatenations of all iterants of the jagged
  vector.
- `offsets`: a vector of offsets into `data`, delineating the start and end of each iterant.
  The `i`th element of the jagged vector is a view into `data` from `offsets[i]` to
  `offsets[i+1]-1`. The last element of `offsets` must be equal to `length(data) + 1`, and
  the first element must be `1`.
"""
function JaggedVector{T}(data::Vector{T}, offsets::Vector{Int}) where T
    issorted(offsets) || throw(ArgumentError("`offsets` must be sorted"))
    first(offsets) == 1 || throw(ArgumentError("`offsets[1]` ≠ 1"))
    last(offsets) == length(data) + 1 || throw(ArgumentError("`offsets[end] ≠ length(data) + 1`"))
    return JaggedVector{T}(data, offsets, Val(:unsafe))
end
JaggedVector(data::Vector{T}, offsets::Vector{Int}) where T = JaggedVector{T}(data, offsets)

# ::: constructor for "empty" jagged vector with specified element lengths :::
"""
    JaggedVector{T}(lengths::AbstractVector{<:Integer})

Return a `JaggedVector{T}` with iterant lengths given by `lengths`, but uninitialized data.
This is useful for creating empty jagged vectors with known iterant lengths but
uninitialized iterant data. The elements of `lengths` must be non-negative valued.

## Example
The following creates a 3-element jagged vector, whose iterants have lengths 3, 2, and 4,
respectively, and with uninitialized data for each iterant:
```julia
jv = JaggedVector{Int}([3, 2, 4])
```
"""
function JaggedVector{T}(lengths::AbstractVector{<:Integer}) where T
    # create a jagged vector with specified vector lengths, but uninitialized data (i.e., to
    # create empty jagged vectors, where we know the length of the vectors it iterates, but
    # don't specify the data of the iterants initially)
    offsets = Vector{Int}(undef, length(lengths) + 1)
    last_offset = offsets[1] = 1
    for (i, l) in enumerate(lengths)
        l_int = convert(Int, l) :: Int
        l_int < 0 && error(lazy"lengths[$i] = $l is negative: must be non-negative")
        last_offset += l_int
        @inbounds offsets[i + 1] = last_offset
    end
    data = Vector{T}(undef, last_offset - 1) # uninitialized data
    return JaggedVector{T}(data, offsets, Val(:unsafe))
end

JaggedVector{T}() where T = JaggedVector{T}(Vector{T}(), [1]) # empty jagged vector
JaggedVector() = JaggedVector{Any}()

"""
    zeros(::Type{JV}, lengths::AbstractVector{<:Integer}) where JV <: JaggedVector

Return a zero-initialized `JaggedVector{T}` with iterant lengths given by `lengths`.
"""
function Base.zeros(::Type{JV}, lengths::AbstractVector{<:Integer}) where JV <: JaggedVector
    jv = JV(lengths)
    fill!(jv.data, zero(eltype(eltype(JV))))
    return jv
end

Base.zero(jv::JaggedVector) = JaggedVector(zero(jv.data), copy(jv.offsets))

# ::: conversion & constructor from vector of vectors :::
"""
    JaggedVector{T}(vs::AbstractVector{<:AbstractVector})

Return a `JaggedVector{T}` from a vector of vectors `vs`, with the iterant of the resulting
jagged vector being the elements of `vs`: i.e., convert the vector of vectors `vs` into a
`JaggedVector`.
"""
function JaggedVector{T}(vs::AbstractVector{<:AbstractVector}) where {T}
    return convert(JaggedVector{T}, vs)
end
JaggedVector(vs::AbstractVector{<:AbstractVector{T}}) where {T} = JaggedVector{T}(vs)

# NB: we define things in terms of `convert` in order to also get implicit conversion
function Base.convert(::Type{JaggedVector{T}}, vs::AbstractVector{<:AbstractVector}) where T
    offsets = Vector{Int}(undef, length(vs) + 1)
    last_offset = offsets[1] = 1
    for i in eachindex(vs)
        last_offset = offsets[i+1] = last_offset + length(vs[i])
    end
    data = Vector{T}(undef, last_offset - 1)
    start = 1
    for i in eachindex(vs)
        v = @inbounds vs[i]
        _possibly_unsafe_copyto!(data, start, v, 1, length(v)) # also handles S→T conversion
        start += length(v)
    end
    return JaggedVector{T}(data, offsets, Val(:unsafe)) :: JaggedVector{T}
end
_possibly_unsafe_copyto!(dst, i, src, j, N) = copyto!(dst, i, src, j, N)
_possibly_unsafe_copyto!(dst::Vector, i, src::Vector, j, N) = unsafe_copyto!(dst, i, src, j, N)

# ::: conversion & constructor from iterable of iterables :::
"""
    JaggedVector{T}(iter)

Return a `JaggedVector{T}` from an iterable of iterables `iter`, whose inner iterable have
elements convertible to type `T`.
"""
JaggedVector{T}(iter) where T = convert(JaggedVector{T}, iter)
JaggedVector(iter) = JaggedVector{eltype(first(first(iter)))}(iter)

function Base.convert(::Type{JaggedVector{T}}, iter) where T
    Nᵥ = length(iter)
    offsets = Vector{Int}(undef, Nᵥ + 1)
    last_offset = offsets[1] = 1
    data = sizehint!(Vector{T}(), Nᵥ * 4) # guess-tuned for typical SymmetryVector use
    for (i, v) in enumerate(iter)
        append!(data, v)
        last_offset = offsets[i + 1] = last_offset + length(v)
    end
    resize!(data, last_offset - 1) # trim the data vector to the correct size
    return JaggedVector{T}(data, offsets, Val(:unsafe)) :: JaggedVector{T}
end

# ::: AbstractArray interface :::
Base.size(jv::JaggedVector) = (length(jv.offsets) - 1,)
@propagate_inbounds function Base.getindex(jv::JaggedVector, i::Int)
    @boundscheck checkbounds(jv, i)
    start = @inbounds jv.offsets[i]
    stop  = @inbounds jv.offsets[i + 1] - 1
    return Base.unsafe_view(jv.data, start:stop)
end
@propagate_inbounds function Base.setindex!(jv::JaggedVector, v, i::Int)
    @boundscheck checkbounds(jv, i)
    start = @inbounds jv.offsets[i]
    stop  = @inbounds jv.offsets[i + 1] - 1
    Nᵢ = stop - start + 1 # length of `jv[i]` before `setindex!`
    V = length(v)         # length of `v` & length of `jv[i]` after `setindex!`

    # update `jv.data` and `jv.offsets` as needed
    if Nᵢ == V
        # inserted `v` has same length as `jv[i]` before `setindex!`, so just copy over
        # (nothing to update in `jv.offsets`, since length is unchanged)
        _possibly_unsafe_copyto!(jv.data, start, v, 1, V)
        nothing
    else # NB: unlike `Vector{Vector{T}}`, this operation is now O(length(jv)), not O(1)
        Δ = V - Nᵢ # negative if `setindex!` reduces length, positive if it increases
        if Δ > 0         # new length is bigger, V > Nᵢ
            # copy over the first `Nᵢ` elements of `v` directly into `jv.data` at start,
            # then insert remaining data in `v[Nᵢ+1:V]` afterwards
            _possibly_unsafe_copyto!(jv.data, start, v, 1, Nᵢ) # safe (cf. Nᵢ < V)
            splice!(jv.data, (start + Nᵢ):(start + Nᵢ - 1), (@inbounds @view v[(Nᵢ + 1):V]))
        else #= Δ < 0 =#     # new length is smaller, V < Nᵢ
            # copy over the first `V` elements of `v` directly into `jv.data` at start,
            # then delete the now-excess elements in `jv.data` after that
            _possibly_unsafe_copyto!(jv.data, start, v, 1, V) # safe (cf. V < Nᵢ)
            deleteat!(jv.data, (stop + Δ + 1):stop) # remove the now-excess elements
        end
        # adjust offsets (subsequent offsets change by `Δ`)
        for j in (i+1):length(jv.offsets)
            @inbounds jv.offsets[j] += Δ
            # NB: the fact that this is done _subsequent_ to updating `jv.data` is not
            #     super, since it means that e.g., `ctrl+c` (or an error) could leave `jv`
            #     in a malformed state; for now, we don't care
        end
    end
    return jv
end
@inline function Base.checkbounds(::Type{Bool}, jv::JaggedVector, i::Integer)
    # equiv. to `i ≥ 1 && i ≤ length(jv)` but maybe faster (borrowed from iterate(::Array))
    return (i - 1)%UInt < length(jv)%UInt
end
function Base.similar(jv::JaggedVector{T}, ::Type{S}=T) where {T,S} # same offsets but new data
    return JaggedVector{S}(similar(jv.data, S), copy(jv.offsets))
end
@inline function Base.iterate(jv::JaggedVector, (i, start)::Tuple{Int, Int} = (1, 1))
    if checkbounds(Bool, jv, i)
        next_start = @inbounds jv.offsets[i + 1]
        v = Base.unsafe_view(jv.data, start:next_start-1)
        return v, (i + 1, next_start)
    else
        return nothing
    end
end
Base.IndexStyle(::Type{<:JaggedVector}) = IndexLinear()

# ::: push!, append!, pop! (and unsupported methods, e.g., pushfirst!) :::
function Base.push!(jv::JaggedVector, vs...)
    append!(jv.data, vs...)
    push!(jv.offsets, new_offsets(jv, vs...)...)
    return jv
end
# `new_offsets` returns the new offsets as a tuple, used by `push!` above; equivalent to 
#    closure = let last_offset = last(jv.offsets), vs = vs
#        (i) -> (last_offset += length(vs[i]); last_offset)
#    end
#    new_offsets = ntuple(closure, Val(length(vs)))
# but the closure & ntuple approach leads to boxing, which this tail-call recursion doesn't
new_offsets(jv::AbstractVector, vs...) = _new_offsets(last(jv.offsets), vs)
_new_offsets(_, ::Tuple{}) = ()
_new_offsets(x, vs) = (x += length(vs[1]); (x, _new_offsets(x, Base.tail(vs))...))

function Base.append!(jv::JaggedVector, iters...)
    for iter in iters
        isempty(iter) && continue
        sizehint!(jv.data, length(jv.data) + sum(length, iter))
        sizehint!(jv.offsets, length(jv.offsets) + length(iter))
        last_offset = last(jv.offsets)
        for v in iter
            append!(jv.data, v)
            last_offset += length(v)
            push!(jv.offsets, last_offset)
        end
        # NB: while it's tempting to do the above with `resize!` instead of `sizehint!` and
        #     with `unsafe_copyto!` instead of `append!` etc., this would not be safe
        #     for malformed input that is not convertible to `jv.data`'s element type `T`.
        #     (e.g., an input that is a vector of strings), which could leave `jv` in a
        #     malformed state that might subsequently lead to out-of-memory access.
    end
    return jv
end

function Base.pop!(jv::JaggedVector)
    isempty(jv) && throw(ArgumentError("`pop!` called on empty `JaggedVector`"))
    v = copy(last(jv))
    resize!(jv.data, length(jv.data) - length(v))
    resize!(jv.offsets, length(jv.offsets) - 1)
    return v
end

Base.pushfirst!(::JaggedVector, vs...) = error("`pushfirst!` is not supported for `JaggedVector`")
Base.popfirst!(::JaggedVector) = error("`popfirst!` is not supported for `JaggedVector`")
Base.popat!(::JaggedVector, i::Integer) = error("`popat!` is not supported for `JaggedVector`")
Base.deleteat!(::JaggedVector, inds) = error("`deleteat!` is not supported for `JaggedVector`")
Base.keepat!(::JaggedVector, inds) = error("`keepat!` is not supported for `JaggedVector`")

# ::: misc niceties & performance optimizations :::
Base.parent(jv::JaggedVector) = jv.data # return the underlying "flat" vector with no copy
Base.reduce(::typeof(vcat), jv::JaggedVector) = copy(jv.data) # copied "flat" vector
function Base.convert(::Type{Vector{Vector{T}}}, jv::JaggedVector) where T
    # unlike `Vector(jv)`, this unaliases the underlying data of `jv`, so that the returned
    # vector is a true copy of the data (and can be modified without affecting `jv`); this
    # is just a small performance optimization over the fallback implementation
    vs = Vector{Vector{T}}(undef, length(jv))
    start = 1
    for i in eachindex(jv)
        next_start = (@inbounds jv.offsets[i + 1])
        N = next_start - start
        v = @inbounds vs[i] = Vector{T}(undef, N)
        Base.unsafe_copyto!(v, 1, jv.data, start, N)
        start = next_start
    end
    return vs
end
Base.copy(jv::JaggedVector{T}) where T = JaggedVector{T}(copy(jv.data), copy(jv.offsets))
Base.empty!(jv::JaggedVector) = (empty!(jv.data); resize!(jv.offsets, 1); jv)
Base.empty(::JaggedVector{T}, ::Type{S}=T) where {T,S} = JaggedVector{S}()

# ::: arithmetic operations :::
function Base.:+(jv1::JaggedVector, jv2::JaggedVector)
    jv1.offsets == jv2.offsets || error(ArgumentError("mismatched offsets: cannot add"))
    return JaggedVector(jv1.data + jv2.data, jv1.offsets)
end
function Base.:-(jv1::JaggedVector, jv2::JaggedVector)
    jv1.offsets == jv2.offsets || error(ArgumentError("mismatched offsets: cannot subtract"))
    return JaggedVector(jv1.data - jv2.data, jv1.offsets)
end
Base.:-(jv::JaggedVector) = JaggedVector(-jv.data, jv.offsets)
Base.:*(s::Number, jv::JaggedVector) = JaggedVector(s * jv.data, jv.offsets)
Base.:*(jv::JaggedVector, s::Number) = JaggedVector(jv.data * s, jv.offsets)
Base.:/(jv::JaggedVector, s::Number) = JaggedVector(jv.data / s, jv.offsets)

# ---------------------------------------------------------------------------------------- #

end # module Jagged