# --- DirectBasis and ReciprocalBasis for crystalline lattices ---
"""
    AbstractBasis <: StaticVector{D, SVector{D, E}}

Abstract supertype of a `D`-dimensional basis in `D`-dimensional space with coordinate
values of type `E`.
"""
abstract type AbstractBasis{D, E} <: StaticVector{D, SVector{D, E}} end

for (T, space_type) in zip((:DirectBasis, :ReciprocalBasis), ("direct", "reciprocal"))
    @eval Bravais begin
        @doc """
            $($T){D, E} <: AbstractBasis{D, E}

        A wrapper type over `D` distinct `D`-dimensional vectors (given as a
        `SVector{D, SVector{D, E}}`), defining a lattice basis in $($space_type)
        space. By default (i.e., if omitted), `E` is `Float64`.
        """
        struct $T{D, E} <: AbstractBasis{D, E}
            vs::SVector{D, SVector{D, E}}
            # ambiguity-resolving methods relative to StaticArrays's methods
            $T{D}(vs::StaticVector{D}) where D = (E=eltype(first(vs)); new{D,E}(convert(SVector{D, SVector{D, E}}, vs)))
            $T{D}(vs::NTuple{D})       where D = (E=eltype(first(vs)); new{D,E}(convert(SVector{D, SVector{D, E}}, vs)))
            $T{D}(vs::AbstractVector)  where D = (E=eltype(first(vs)); new{D,E}(convert(SVector{D, SVector{D, E}}, vs)))
            $T{D,E}(vs::StaticVector{D}) where {D,E} = new{D,E}(convert(SVector{D, SVector{D, E}}, vs))
            $T{D,E}(vs::NTuple{D})       where {D,E} = new{D,E}(convert(SVector{D, SVector{D, E}}, vs))
            $T{D,E}(vs::AbstractVector)  where {D,E} = new{D,E}(convert(SVector{D, SVector{D, E}}, vs))
            # special-casing for D=1 (e.g., to make `$T([1.0])` work)
            $T{1}(vs::StaticVector{1,E}) where E<:Number = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
            $T{1}(vs::NTuple{1,E})       where E<:Number = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
            $T{1}(vs::AbstractVector{E}) where E<:Number = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
            $T{1,E}(vs::StaticVector{1,<:Number}) where E = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
            $T{1,E}(vs::NTuple{1,<:Number})       where E = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
            $T{1,E}(vs::AbstractVector{<:Number}) where E = new{1,E}(SVector((convert(SVector{1, E}, vs),)))
        end
    end
    @eval function convert(::Type{$T{D}}, vs::StaticVector{D, <:StaticVector{D}}) where D
        E = eltype(first(vs))
        $T{D,E}(convert(SVector{D, SVector{D, E}}, vs))
    end
    @eval function convert(::Type{$T{D,E}}, vs::StaticVector{D, <:StaticVector{D, E}}) where {D,E}
        $T{D,E}(convert(SVector{D, SVector{D, E}}, vs))
    end
    @eval $T(vs::StaticVector{D})      where D     = $T{D}(vs)   # resolve more ambiguities, both
    @eval $T(vs::StaticVector{D,E}...) where {D,E} = $T{D,E}(vs) # internally and with StaticArrays,
    @eval $T(vs::NTuple{D})            where D     = $T{D}(vs)   # and make most reasonable accessor
    @eval $T(vs::NTuple{D,E}...)       where {D,E} = $T{D,E}(vs) # patterns functional
    @eval $T(vs::AbstractVector)    = $T{length(vs)}(vs) # [type-unstable]
    @eval $T(vs::AbstractVector...) = $T{length(vs), eltype(eltype(promote(vs...)))}(promote(vs...)) # [type-unstable]
end

parent(Vs::AbstractBasis) = Vs.vs
# define the AbstractArray interface for AbstractBasis{D}
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

if VERSION > v"1.9.0-DEV.1163" 
    # since https://github.com/JuliaLang/julia/pull/43334, Julia defines its own `stack`;
    # however, it is still much slower than a naive implementation based on `reduce` cf.
    # https://github.com/JuliaLang/julia/issues/52590. As such, we extend `Base.stack` even
    # on more recent versions; when the issue is fixed, it would be enough to only define
    # `stack` on earlier versions of Julia, falling back to `Base.stack` on later versions.
    import Base: stack
end
"""
    stack(Vs::AbstractBasis)

Return a matrix `[Vs[1] Vs[2] .. Vs[D]]` from `Vs::AbstractBasis{D}`, i.e., the matrix whose
columns are the basis vectors of `Vs`.
"""
stack(Vs::AbstractBasis) = reduce(hcat, parent(Vs))
stack(Vs::AbstractBasis{1, E}) where E = SMatrix{1, 1, E, 1}(@inbounds Vs[1][1])

"""
    volume(Vs::AbstractBasis)

Return the volume ``V`` of the unit cell associated with the basis `Vs::AbstractBasis{D}`.

The volume is computed as ``V = \\sqrt{\\mathrm{det}\\mathbf{G}}`` with with ``\\mathbf{G}``
denoting the metric matrix of `Vs` (cf. the International Tables of Crystallography, 
Volume A, Section 5.2.2.3).

See also [`metricmatrix`](@ref).
"""
volume(Vs::AbstractBasis) = sqrt(det(metricmatrix(Vs)))

"""
    metricmatrix(Vs::AbstractBasis)

Return the (real, symmetric) metric matrix of a basis `Vs`, i.e., the matrix with elements
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

"""
    dualtype(Vs::AbstractBasis)
    dualtype(::Type{<:AbstractBasis})

Return the dual type of a basis `Vs`, i.e., the type of the dual basis of `Vs`.

The default basis types `DirectBasis` and `ReciprocalBasis` are dual to each other. If no
dual type is defined, `nothing` is returned.
"""
dualtype(::T) where T = dualtype(T)
dualtype(::Type{<:DirectBasis{D,E}}) where {D,E} = ReciprocalBasis{D,_invtype(E)}
dualtype(::Type{<:ReciprocalBasis{D,E}}) where {D,E} = DirectBasis{D,_invtype(E)}
dualtype(::Type) = nothing
_invtype(::Type{T}) where T<:Number = typeof(inv(oneunit(T)))

# ---------------------------------------------------------------------------------------- #

"""
    AbstractPoint{D, T} <: StaticVector{D, T}

Abstract supertype of a `D`-dimensional point with elements of type `T`.
"""
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
        associated $($BT). To convert to Cartesian coordinates, see [`cartesianize`](@ref).
        """
        struct $PT{D} <: AbstractPoint{D, Float64}
            v::SVector{D, Float64}
            # ambiguity-resolving methods relative to StaticArray's
            $PT{D}(v::StaticVector{D, <:Real}) where D = new{D}(convert(SVector{D, Float64}, v))
            $PT{D}(v::NTuple{D, Real}) where D         = new{D}(convert(SVector{D, Float64}, v))
            $PT{D}(v::AbstractVector{<:Real}) where D  = new{D}(convert(SVector{D, Float64}, v))
        end
        @eval function convert(::Type{$PT{D}}, v::StaticVector{D, <:Real}) where D
            $PT{D}(convert(SVector{D, Float64}, v))
        end
        @eval $PT(v::StaticVector{D}) where D = $PT{D}(v) # resolve internal/StaticArrays
        @eval $PT(v::NTuple{D, Real}) where D = $PT{D}(v) # ambiguities & add accessors
        @eval $PT(v::AbstractVector) = $PT{length(v)}(v)
        @eval $PT(v::Real...) = $PT{length(v)}(v)
    end
end

# arithmetic
Base.:+(v::T, w::T) where T<:AbstractPoint = T(v.v + w.v)
Base.:-(v::T, w::T) where T<:AbstractPoint = T(v.v - w.v)
Base.:-(v::T) where T<:AbstractPoint = T(-v.v)
Base.:*(v::T, c::Real) where T<:AbstractPoint = T(v.v * c)
Base.:*(c::Real, v::AbstractPoint) = v*c
Base.:/(v::T, c::Real) where T<:AbstractPoint = T(v.v / c)
Base.zero(::Type{<:T}) where T<:AbstractPoint{D} where D = T(zero(SVector{D, Float64}))