@deprecate reciprocalbasis(Rs::Union{AbstractVector, NTuple}) dualbasis(DirectBasis(Rs))
@deprecate reciprocalbasis(Rs::DirectBasis) dualbasis(Rs)
@deprecate directbasis(Rs::ReciprocalBasis) dualbasis(Rs)