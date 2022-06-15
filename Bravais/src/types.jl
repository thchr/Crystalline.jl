# --- DirectBasis and ReciprocalBasis for crystalline lattices ---
"""
    AbstractBasis <: StaticVector{D, SVector{D,Float64}}

Abstract supertype of a `D`-dimensional basis in `D`-dimensional space.
"""
abstract type AbstractBasis{D, T} <: StaticVector{D, SVector{D, T}} end
for (T, space_type) in zip((:DirectBasis, :ReciprocalBasis), ("direct", "reciprocal"))
    @eval begin
        @doc """
            $($T){D} <: AbstractBasis{D}

        A wrapper type over `D` distinct `D`-dimensional vectors (given as a
        `SVector{D, SVector{D,Float64}}`), defining a lattice basis in $($space_type)
        space.
        """
        struct $T{D} <: AbstractBasis{D, Float64}
            vs::SVector{D, SVector{D, Float64}}
            # ambiguity-resolving methods relative to StaticArrays's methods
            $T{D}(vs::StaticVector{D}) where D = new{D}(convert(SVector{D, SVector{D, Float64}}, vs))
            $T{D}(vs::NTuple{D})       where D = new{D}(convert(SVector{D, SVector{D, Float64}}, vs))
            $T{D}(vs::AbstractVector)  where D = new{D}(convert(SVector{D, SVector{D, Float64}}, vs))
            # special-casing for D=1 (e.g., to make `$T([1.0])` work)
            $T{1}(vs::StaticVector{1,<:Real}) = new{1}(SVector((convert(SVector{1, Float64}, vs),)))
            $T{1}(vs::NTuple{1,Real})         = new{1}(SVector((convert(SVector{1, Float64}, vs),)))
            $T{1}(vs::AbstractVector{<:Real}) = new{D}(SVector((convert(SVector{1, Float64}, vs),)))
        end
    end
    @eval function convert(::Type{$T{D}}, vs::StaticVector{D, <:StaticVector{D, <:Real}}) where D
        $T{D}(convert(SVector{D, SVector{D, Float64}}, vs))
    end
    @eval $T(vs::StaticVector{D})    where D = $T{D}(vs) # resolve more ambiguities, both
    @eval $T(vs::StaticVector{D}...) where D = $T{D}(vs) # internally and with StaticArrays,
    @eval $T(vs::NTuple{D})          where D = $T{D}(vs) # and make most reasonable accessor
    @eval $T(vs::NTuple{D}...)       where D = $T{D}(vs) # patterns functional
    @eval $T(vs::AbstractVector)             = $T{length(vs)}(vs) # [type-unstable]
    @eval $T(vs::AbstractVector...)          = $T{length(vs)}(promote(vs...))
end

parent(Vs::AbstractBasis) = Vs.vs
# define the AbstractArray interface for DirectBasis{D}
@propagate_inbounds getindex(Vs::AbstractBasis, i::Int) = parent(Vs)[i]
size(::AbstractBasis{D}) where D = (D,)
IndexStyle(::Type{<:AbstractBasis}) = IndexLinear()

_angle(rA, rB) = acos(dot(rA, rB) / (norm(rA) * norm(rB)))
angles(Rs::AbstractBasis{2}) = _angle(Rs[1], Rs[2])
function angles(Rs::AbstractBasis{3})
    α = _angle(Rs[2], Rs[3])
    β = _angle(Rs[3], Rs[1])
    γ = _angle(Rs[1], Rs[2])
    return α, β, γ
end
angles(::AbstractBasis{D}) where D = _throw_invalid_dim(D)

"""
    stack(Vs::AbstractBasis)

Return a matrix `[Vs[1] Vs[2] .. Vs[D]]` from `Vs::AbstractBasis{D}`, i.e. the matrix whose
columns are the basis vectors of `Vs`.
"""
stack(Vs::AbstractBasis) = reduce(hcat, parent(Vs))
# TODO: At some point, this should hopefully no longer be necessary to do manually (and
# `stack` may end up exported by Base): https://github.com/JuliaLang/julia/issues/21672

"""
    volume(Vs::AbstractBasis)

Return the volume ``V`` of the unit cell associated with the basis `Vs::AbstractBasis{D}`.

The volume is computed as ``V = \\sqrt{\\mathrm{det}G}`` with with ``G`` denoting the metric
matrix of `Vs` (cf. the International Tables of Crystallography, Volume A, Section 5.2.2.3).

See also [`Bravais.metricmatrix`](@ref).
"""
volume(Vs::AbstractBasis) = sqrt(det(metricmatrix(Vs))) # TODO: wrong. TODO: export.

"""
    metricmatrix(Vs::AbstractBasis)

Return the (real, symmetric) metric matrix of a basis `Vs`, i.e. the matrix with elements
``G_{ij} =`` `dot(Vs[i], Vs[j])`, as defined in the International Tables of Crystallography,
Volume A, Section 5.2.2.3.

Equivalently, this is the Gram matrix of `Vs`, and so can also be expressed as `Vm' * Vm`
with `Vm` denoting the columnwise concatenation of the basis vectors in `Vs`.

See also [`volume`](@ref).
"""
function metricmatrix(Vs::AbstractBasis{D}) where D
    Vm = stack(Vs)
    return Vm' * Vm # equivalent to [dot(v, w) for v in Vs, w in Vs]
end

# ---------------------------------------------------------------------------------------- #

abstract type AbstractPoint{D, T} <: StaticVector{D, T} end

parent(p::AbstractPoint) = p.v

@propagate_inbounds getindex(v::AbstractPoint, i::Int) = parent(v)[i]
size(::AbstractPoint{D}) where D = (D,)
IndexStyle(::Type{<:AbstractPoint}) = IndexLinear()

for (PT, BT, space_type) in zip((:DirectPoint, :ReciprocalPoint),
                                (:DirectBasis, :ReciprocalBasis),
                                ("direct", "reciprocal"))
    @eval begin
        @doc """
            $($PT){D} <: AbstractPoint{D}

        A wrapper type over an `SVector{D, Float64}`, defining a single point in
        `D`-dimensional $($space_type) space. 
        
        The coordinates of a $($PT) are generally assumed specified relative to an
        associated $($BT).
        """
        struct $PT{D} <: AbstractPoint{D, Float64}
            v::SVector{D, Float64}
            # ambiguity-resolving methods relative to StaticArray's
            $PT{D}(v::StaticVector{D, <:Real}) where D = new{D}(convert(SVector{D, Float64}, v))
            $PT{D}(v::NTuple{D, Real}) where D         = new{D}(convert(SVector{D, Float64}, v))
            $PT{D}(v::AbstractVector{<:Real}) where D  = new{D}(convert(SVector{D, Float64}, v))
        end
        @eval function convert(::Type{$PT{D}}, v::AbstractVector{<:Real}) where D
            $PT{D}(convert(SVector{D, Float64}, v))
        end
        @eval $PT(v::StaticVector{D}) where D = $PT{D}(v) # resolve internal/StaticArrays
        @eval $PT(v::NTuple{D, Real}) where D = $PT{D}(v) # ambiguities & add accessors
        @eval $PT(v::AbstractVector) = $PT{length(v)}(v)
        @eval $PT(v::Real...) = $PT{length(v)}(v)
    end
end