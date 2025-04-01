# ---------------------------------------------------------------------------------------- #
# DirectBasis

function Base.show(io::IO, ::MIME"text/plain", Vs::AbstractBasis)
    _print_basis_typename(io, typeof(Vs))
    print(io, " ($(crystalsystem(Vs))):")
    for V in Vs
        print(io, "\n ", V)
    end
end

@inline _print_basis_typename(io, T::Type{<:AbstractBasis}) = print(io, T)
@inline function _print_basis_typename(io, ::Type{T}) where T <: AbstractBasis{D, Float64} where D
    print(io, T.name.name, "{", D, "}")
end