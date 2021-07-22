module SmithNormalForm

using LinearAlgebra
using SparseArrays
using Base.CoreLogging

import Base: show, summary
import LinearAlgebra: diagm, diag

export snf, smith, Smith

include("bezout.jl")
include("snf.jl")

# ---------------------------------------------------------------------------------------- #

struct Smith{P,Q<:AbstractMatrix{P},V<:AbstractVector{P}} <: Factorization{P}
    S::Q
    Sinv::Q
    T::Q
    Tinv::Q
    SNF::V
    Smith{P,Q,V}(S::AbstractMatrix{P}, Sinv::AbstractMatrix{P},
                 T::AbstractMatrix{P}, Tinv::AbstractMatrix{P},
                 D::AbstractVector{P}) where {P,Q,V} = new(S, Sinv, T, Tinv, D)
end
Smith(S::AbstractMatrix{P}, T::AbstractMatrix{P}, SNF::AbstractVector{P}) where {P} =
    Smith{P,typeof(S),typeof(SNF)}(S, similar(S, 0, 0), T, similar(T, 0, 0), SNF)

# ---------------------------------------------------------------------------------------- #

"""
    smith(X::AbstractMatrix{P}; inverse::Bool=true) --> Smith{P,Q,V}

Return a Smith normal form of an integer matrix `X` as a `Smith` structure (of element type
`P`, matrix type `Q`, and invariant factor type `V`).

The Smith normal form is well-defined for any matrix ``m×n`` matrix `X` with elements in a
principal domain (PID; e.g., integers) and provides a decomposition of `X` into ``m×m``,
``m×n``, and `S`, `Λ`, and `T` as `X = SΛT`, where `Λ` is a diagonal matrix with entries 
("invariant factors") Λᵢ ≥ Λᵢ₊₁ ≥ 0 with nonzero entries divisible in the sense Λᵢ | Λᵢ₊₁.
The invariant factors can be obtained from [`diag(::Smith)`](@ref).

`S` and `T` are invertible matrices; if the keyword argument `inverse` is true (default),
the inverse matrices are computed and returned as part of the `Smith` factorization.
"""
function smith(X::AbstractMatrix{P}; inverse::Bool=true) where {P}
    S, T, D, Sinv, Tinv = snf(X, inverse=inverse)
    SNF = diag(D)
    return Smith{P, typeof(X), typeof(SNF)}(S, Sinv, T, Tinv, SNF)
end

# ---------------------------------------------------------------------------------------- #

"""
    diagm(F::Smith) --> AbstractMatrix

Return the Smith normal form diagonal matrix `Λ` from a factorization `F`.
"""
function diagm(F::Smith{P}) where P
    rows = size(F.S, 1)
    cols = size(F.T, 1)
    D    = issparse(F.SNF) ? spzeros(P, rows, cols) : zeros(P, rows, cols)
    for (i,Λᵢ) in enumerate(diag(F))
        D[i,i] = Λᵢ
    end
    return D
end

"""
    diag(F::Smith) --> AbstractVector

Return the invariant factors (or, equivalently, the elementary divisors) of a Smith normal
form `F`.
"""
diag(F::Smith) = F.SNF

summary(io::IO, F::Smith{P,Q,V}) where {P,Q,V} = print(io, "Smith{", P, ",", Q, ",", V, "}")
function show(io::IO, ::MIME"text/plain", F::Smith)
    summary(io, F)
    println(io, " with:")
    Base.print_matrix(io, diagm(F), " Λ = "); println(io)
    Base.print_matrix(io, F.S, " S = "); println(io)
    Base.print_matrix(io, F.T, " T = ")
    return nothing
end

end # module
