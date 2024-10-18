# ---------------------------------------------------------------------------------------- #
# DirectBasis

function Base.show(io::IO, ::MIME"text/plain", Vs::AbstractBasis)
    print(io, typeof(Vs))
    print(io, " ($(crystalsystem(Vs))):")
    for V in Vs
        print(io, "\n ", V)
    end
end