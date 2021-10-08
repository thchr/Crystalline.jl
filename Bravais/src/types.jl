# --- DirectBasis and ReciprocalBasis for crystalline lattices ---
"""
    AbstractBasis <: StaticVector{D, SVector{D,Float64}}

Abstract supertype of a `D`-dimensional basis in `D`-dimensional space.
"""
abstract type AbstractBasis{D} <: StaticVector{D, SVector{D,Float64}} end
for (T, space_type) in zip((:DirectBasis, :ReciprocalBasis), ("direct", "reciprocal"))
    @eval begin
        @doc """
            $($T){D} <: AbstractBasis{D}

        A wrapper type over `D` distinct `D`-dimensional vectors (given as a
        `SVector{D, SVector{D,Float64}}`), defining a lattice basis in $($(space_type))
        space.
        """
        struct $T{D} <: AbstractBasis{D}
            vecs::SVector{D, SVector{D,Float64}}
            $T{D}(vecs::SVector{D, SVector{D,Float64}}) where D = new{D}(vecs)
            $T(vecs::SVector{D, SVector{D,Float64}}) where D    = new{D}(vecs)
        end
    end
    @eval function convert(::Type{$T{D}}, Vs::StaticVector{D, <:StaticVector{D, <:Real}}) where D
        $T{D}(convert(SVector{D, SVector{D, Float64}}, Vs))
    end
    @eval $T{D}(Vs::NTuple{D, SVector{D, Float64}}) where D = $T{D}(SVector{D}(Vs))
    @eval $T(Vs::NTuple{D, SVector{D, Float64}}) where D = $T{D}(Vs)
    @eval $T{D}(Vs::NTuple{D, NTuple{D,<:Real}}) where D = $T{D}(SVector{D,Float64}.(Vs))
    @eval $T(Vs::NTuple{D, NTuple{D,<:Real}}) where D = $T{D}(Vs)
    @eval $T{D}(Vs::NTuple{D, <:AbstractVector{<:Real}}) where D = $T{D}(Vs...)
    @eval $T(Vs::NTuple{D, <:AbstractVector{<:Real}}) where D = $T{D}(Vs...)
    @eval $T{D}(Vs::AbstractVector{<:AbstractVector{<:Real}}) where D = $T{D}(Vs...)
    @eval $T(Vs::AbstractVector{<:AbstractVector{<:Real}}) = $T(Vs...)
    @eval $T{D}(Vs::AbstractVector{<:Real}...) where D = $T{D}(convert(SVector{D, SVector{D, Float64}}, Vs))
    @eval $T(Vs::AbstractVector{<:Real}...) = $T{length(Vs)}(Vs...)
    @eval $T{D}(Vs::StaticVector{D,<:Real}...) where D = $T{D}(Vs) # resolve ambiguities w/
    @eval $T(Vs::StaticVector{D,<:Real}...) where D = $T{D}(Vs)    # `::StaticArray` methods
end

parent(Vs::AbstractBasis) = Vs.vecs
# define the AbstractArray interface for DirectBasis{D}
@propagate_inbounds getindex(Vs::AbstractBasis, i::Int) = parent(Vs)[i]
size(::AbstractBasis{D}) where D = (D,)
IndexStyle(::Type{<:AbstractBasis}) = IndexLinear()

_angle(rA,rB) = acos(dot(rA,rB)/(norm(rA)*norm(rB)))
function angles(Rs::AbstractBasis{D}) where D
    D == 1 && return nothing
    γ = _angle(Rs[1], Rs[2])
    if D == 3
        α = _angle(Rs[2], Rs[3])
        β = _angle(Rs[3], Rs[1])
        return α,β,γ
    end
    return γ
end

"""
    stack(Vs::AbstractBasis{D}) where D

Return a matrix `[Vs[1] Vs[2] .. Vs[D]]` from `Vs::AbstractBasis{D}`, i.e. the matrix whose
columns are the basis vectors in `Vs`. 
"""
stack(Vs::AbstractBasis{D}) where D = reduce(hcat, parent(Vs))
# TODO: At some point, this should hopefully no longer be necessary to do manually (and
# `stack` may end up exported by Base): https://github.com/JuliaLang/julia/issues/21672