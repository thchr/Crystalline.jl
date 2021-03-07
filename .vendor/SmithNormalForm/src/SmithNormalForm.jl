module SmithNormalForm

# import Base: show, getindex
using LinearAlgebra
using SparseArrays
using Base.CoreLogging

import Base: show
import LinearAlgebra: diagm

export snf, smith

include("bezout.jl")
include("snf.jl")

struct Smith{P,Q<:AbstractMatrix} <: Factorization{P}
    S::Q
    Sinv::Q
    T::Q
    Tinv::Q
    SNF::AbstractVector{P}
    Smith{P,Q}(S::AbstractMatrix{P}, Sinv::AbstractMatrix{P},
               T::AbstractMatrix{P}, Tinv::AbstractMatrix{P},
               A::AbstractVector{P}) where {P,Q} = new(S, Sinv, T, Tinv, A)
end
Smith(S::AbstractMatrix{P}, T::AbstractMatrix{P}, A::AbstractVector{P}) where {P} =
    Smith{P,typeof(A)}(S, similar(S, 0, 0), T, similar(T, 0, 0), A)

function smith(X::AbstractMatrix{P}; inverse=true) where {P}
    S, T, A, Sinv, Tinv = snf(X, inverse=inverse)
    return Smith{P, typeof(X)}(S, Sinv, T, Tinv, diag(A))
end

"""Retrive SNF matrix from the factorization"""
function diagm(F::Smith{P,Q}) where {P,Q}
    base = issparse(F.SNF) ? spzeros(P, size(F.S,1), size(F.T,1)) : zeros(P, size(F.S,1), size(F.T,1))
    for i in 1:length(F.SNF)
        base[i,i] = F.SNF[i]
    end
    return base
end

function show(io::IO, F::Smith)
    println(io, "Smith normal form:")
    show(io, diagm(F))
end

end # module
