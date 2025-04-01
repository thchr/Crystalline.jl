# ---------------------------------------------------------------------------------------- #
# AbstractBasis

function Base.show(io::IO, ::MIME"text/plain", Vs::AbstractBasis)
    _print_basis_typename(io, Vs)
    print(io, " ($(crystalsystem(Vs))):")
    for V in Vs
        print(io, "\n ")
        # same as `show(io, Vs)`, but suppress explicit printing of typeinfo before "[...]"
        # (in case it isn't a primitive type; e.g., for unitful types)
        Base.show_delim_array(io, V, '[', ',', ']', false)
    end
end

@inline _print_basis_typename(io, Vs::AbstractBasis) = print(io, typeof(Vs))
@inline function _print_basis_typename(io, Vs::AbstractBasis{D, Float64}) where D
    # omit `Float64` from type name: `DirectBasis{3, Float64}` → `DirectBasis{3}`
    print(io, _basis_typename(Vs), "{", D, "}")
end
@inline _basis_typename(::DirectBasis) = "DirectBasis"
@inline _basis_typename(::ReciprocalBasis) = "ReciprocalBasis"