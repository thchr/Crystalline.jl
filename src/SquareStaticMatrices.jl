module SquareStaticMatrices
# This module defines `SqSMatrix{D,T}` an analogue of `SMatrix{D,D,T,D*D}` that only needs
# a single dimensional type parameter, avoiding the redundant type parameter `L=D*D`. 
# This is desirable if we want to have a square static matrix in some struct but want to
# avoid dragging the redundant `L` parameter around *everywhere*. 
# The trick here is to store the matrix as a nested `NTuple{D,NTuple{D,T}}` rather than as
# `NTuple{D*D, T}` as in StaticArrays.
# The performance is pretty much on par with an `SMatrix` (e.g. there are no allocations),
# only very slightly slower because we didn't bother to implemented certain things as 
# generated functions. Other than that, the memory layout should be essentially the same.
# Some relevant discussion of related problems is e.g.:
#   https://github.com/JuliaLang/julia/issues/18466 
#   https://discourse.julialang.org/t/addition-to-parameter-of-parametric-type/20059

using LinearAlgebra: checksquare
using Base: @propagate_inbounds

import StaticArrays: SMatrix, MMatrix
import Base: convert, eltype, size, getindex, firstindex, lastindex, eachcol, one, zero,
             IndexStyle

export SqSMatrix

# ---------------------------------------------------------------------------------------- #
# struct definition

# an equivalent of `SMatrix{D,D,T}`, but requiring only a single dimension type
# parameter, rather than two
struct SqSMatrix{D, T} <: AbstractMatrix{T}
    cols::NTuple{D, NTuple{D, T}} # tuple of columns (themselves stored as tuples)
end

# ---------------------------------------------------------------------------------------- #
# AbstractArray interface
size(::SqSMatrix{D}) where D = (D,D)
eltype(::SqSMatrix{D,T}) where {D,T} = T
firstindex(::SqSMatrix) = 1
lastindex(::SqSMatrix{D}) where D = D
lastindex(::SqSMatrix{D}, d::Int64) where D = d == 1 ? D : (d == 2 ? D : 1)
eachcol(A::SqSMatrix) = A.cols
@propagate_inbounds function getindex(A::SqSMatrix{D}, i::Int) where D
    @boundscheck 1 ≤ i ≤ D*D || throw(BoundsError(A, (i,)))
    i, j = (idx+D-1)÷D, mod1(idx, D)
    return @inbounds A.cols[j][i]
end
@propagate_inbounds function getindex(A::SqSMatrix{D}, i::Integer, j::Integer) where D
    @boundscheck (1 ≤ i ≤ D && 1 ≤ j ≤ D) || throw(BoundsError(A, (i,j)))
    return @inbounds A.cols[j][i]
end
IndexStyle(::Type{<:SqSMatrix}) = IndexCartesian()

# ---------------------------------------------------------------------------------------- #
# constructors and converters 
@propagate_inbounds @inline function convert(::Type{SqSMatrix{D, T}}, A::AbstractMatrix) where {D,T}
    # TODO: this could be a little bit faster if we used a generated function as they do in
    #       for the StaticArrays ` unroll_tuple(a::AbstractArray, ::Length{L})` method...
    @boundscheck checksquare(A) == D
    cols = ntuple(Val{D}()) do j
        ntuple(i->convert(T, @inbounds A[i,j]), Val{D}())
    end
    SqSMatrix{D,T}(cols)
end

# allow an NTuple{D,NTuple{D,T}} to be converted automatically to SqSMatrix{D,T} if relevant
function convert(::Type{SqSMatrix{D, T}}, cols::NTuple{D, NTuple{D, T}}) where {D,T}
    SqSMatrix{D,T}(cols)
end

@propagate_inbounds SqSMatrix{D,T}(A::AbstractMatrix) where {D,T}  = convert(SqSMatrix{D,T}, A)
@propagate_inbounds SqSMatrix{D}(A::AbstractMatrix{T}) where {D,T} = convert(SqSMatrix{D,T}, A)
@propagate_inbounds function SqSMatrix(A::AbstractMatrix{T}) where T
    D = checksquare(A)
    @inbounds SqSMatrix{D}(A::AbstractMatrix{T})
end

function flatten_nested(cols::NTuple{D, NTuple{D, T}}) where {D,T}
    ntuple(Val{D*D}()) do idx
       i, j = (idx+D-1)÷D, mod1(idx, D)
       @inbounds cols[i][j]
    end
end
flatten(A::SqSMatrix{D,T}) where {D,T} = flatten_nested(A.cols)
# equivalent recursive implementation (equal performance for small N but much worse for large N)
# flatten_nested(cols::NTuple{N, NTuple{N, T}}) where {N,T} = flatten_nested(cols...)
# flatten_nested(x) = x
# flatten_nested(x,y) = (x...,y...)
# flatten_nested(x,y,z...) = (x..., flatten_nested(y, z...)...)

for M in (:SMatrix, :MMatrix)
    @eval $M(A::SqSMatrix{D, T}) where {D,T} = $M{D, D, T}(flatten(A))
end

# use a generated function to stack a "square" NTuple very efficiently: allows efficient
# conversion from an `SMatrix` to an `SqSmatrix`
@generated function stack_square_tuple(xs::NTuple{N, T}) where {N,T}
    D = isqrt(N)
    if D*D != N 
        return :(throw(DomainError($N, "called with tuple of non-square length $N")))
    else
        exs=[[:(xs[$i+($j-1)*$D]) for i in 1:D] for j in 1:D]
        quote
            Base.@_inline_meta
            @inbounds return $(Expr(:tuple, (map(ex->Expr(:tuple, ex...), exs))...))
        end
    end
end
@inline function convert(::Type{SqSMatrix{D,T}}, A::Union{SMatrix{D,D,T}, MMatrix{D,D,T}}) where {D,T}
    SqSMatrix{D,T}(@inbounds stack_square_tuple(A.data))
end
SqSMatrix(A::SMatrix{D,D,T}) where {D,T} = convert(SqSMatrix{D,T}, A)
SqSMatrix{D,T}(A::SMatrix{D,D,T′}) where {D,T,T′} = convert(SqSMatrix{D,T}, convert(SMatrix{D,D,T}, A))


# ---------------------------------------------------------------------------------------- #
# linear algebra and other methods

# by the magic of julia's compiler, this reduces to the exact same machine code as if we
# were writing out ntuples ourselves (i.e. the initial construction as an SMatrix is
# optimized out completely)
zero(::Type{SqSMatrix{D,T}}) where {D,T} = SqSMatrix{D,T}(zero(SMatrix{D,D,T}))
one( ::Type{SqSMatrix{D,T}}) where {D,T} = SqSMatrix{D,T}(one(SMatrix{D,D,T}))

end # module