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

import StaticArrays: SMatrix
import Base: convert, eltype, size, getindex, firstindex, lastindex

export SqSMatrix

# ---------------------------------------------------------------------------------------- #
# struct definition

# an equivalent of `SMatrix{D,D,T,D*D}`, but requiring only a single dimension type
# parameter, rather than two
struct SqSMatrix{D, T} <: AbstractMatrix{T}
    data::NTuple{D, NTuple{D, T}} # tuple of column-elements (as tuple elements)
end
function getindex(A::SqSMatrix{D}, i, j) where D
    @boundscheck (1 ≤ i ≤ D && 1 ≤ j ≤ D) || throw(BoundsError(A, (i,j)))
    return @inbounds A.data[j][i]
end

# ---------------------------------------------------------------------------------------- #
# AbstractArray interface
size(::SqSMatrix{D}) where D = (D,D)
eltype(::SqSMatrix{D,T}) where D where T = T
firstindex(::SqSMatrix) = 1
lastindex(::SqSMatrix{D}) where D = D
lastindex(::SqSMatrix{D}, d::Int64) where D = d == 1 ? D : (d == 2 ? D : 1)

# ---------------------------------------------------------------------------------------- #
# constructors and converters 
@propagate_inbounds function convert(::Type{SqSMatrix{D, T}}, A::AbstractMatrix{T}) where D where T
    # TODO: this could be a little bit faster if we used a generated function as they do in
    #       for the StaticArrays ` unroll_tuple(a::AbstractArray, ::Length{L})` method...
    @boundscheck checksquare(A) == D
    data = @inbounds ntuple(Val{D}()) do j
        ntuple(i->A[i,j], Val{D}())
    end
    SqSMatrix{D,T}(data)
end

@propagate_inbounds function SqSMatrix{D}(A::AbstractMatrix{T}) where D where T
    convert(SqSMatrix{D,T}, A)
end
@propagate_inbounds function SqSMatrix(A::AbstractMatrix{T}) where T
    D = checksquare(A)
    @inbounds SqSMatrix{D}(A::AbstractMatrix{T})
end

function flatten_nested_ntuples(data::NTuple{D, NTuple{D, T}}) where D where T
    ntuple(Val{D*D}()) do idx
       i, j = (idx+D-1)÷D, mod1(idx, D)
       data[i][j]
    end
end
flatten(A::SqSMatrix{D,T}) where D where T = flatten_nested_ntuples(A.data)

function SMatrix(A::SqSMatrix{D, T}) where D where T
    SMatrix{D, D, T, D*D}(flatten(A))
end

end # module