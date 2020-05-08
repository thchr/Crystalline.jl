# === frequent error messages ===

@noinline _throw_1d2d_not_yet_implemented(D::Integer) = 
    throw(DomainError(D, "dimensions D≠3 not yet supported"))
@noinline _throw_invaliddim(D::Integer) = 
    throw(DomainError(D, "dimension must be 1, 2, or 3"))
@noinline _throw_invalidcntr(cntr::Char) = 
    throw(DomainError(cntr, "centering must be {P,I,F,R,A,C} in 3D, {p,c} in 2D, or p in 1D"))

# === string manipulation ===
""" 
    parsefraction(str::AbstractString)

Parse a string `str`, allowing fraction inputs (e.g. `"1/2"`), return as `Float64`.
"""
function parsefraction(str::AbstractString)
    slashidx = findfirst(==('/'), str)
    if slashidx === nothing
        return parse(Float64, str)
    else
        num = SubString(str, firstindex(str), prevind(str, slashidx))
        den = SubString(str, nextind(str, slashidx), lastindex(str))
        return parse(Float64, num)/parse(Float64, den)
    end
end

"""
    fractionify!(io::IO, x::Real, forcesign::Bool=true, tol::Real=1e-4)

Write a string representation of the nearest fraction (within a tolerance `tol`) of `x` to 
`io`. If `forcesign` is true, the sign character of `x` is printed whether `+` or `-` 
(otherwise, only printed if `-`).
"""
function fractionify!(io::IO, x::Number, forcesign::Bool=true, tol::Real=1e-4)
    if forcesign || signbit(x)
        print(io, signaschar(x))
    end
    t = rationalize(float(x), tol=tol) # convert to "minimal" Rational fraction (within nearest 1e-4 neighborhood)
    if !isinteger(t)
        print(io, abs(numerator(t)), '/', denominator(t))
    else
        print(io, abs(numerator(t)))
    end
    return nothing
end
function fractionify(x::Number, forcesign::Bool=true, tol::Real=1e-4)
    buf = IOBuffer()
    fractionify!(buf, x, forcesign, tol)
    return String(take!(buf))
end

# --- UNICODE FUNCTIONALITY ---
const SUBSCRIPT_MAP = Dict('1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', '5'=>'₅',  # digits
                           '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉', '0'=>'₀',
                           'a'=>'ₐ', 'e'=>'ₑ', 'h'=>'ₕ', 'i'=>'ᵢ', 'j'=>'ⱼ',  # letters (missing several)
                           'k'=>'ₖ', 'l'=>'ₗ',  'm'=>'ₘ', 'n'=>'ₙ', 'o'=>'ₒ', 
                           'p'=>'ₚ', 'r'=>'ᵣ', 's'=> 'ₛ', 't'=>'ₜ', 'u'=>'ᵤ', 
                           'v'=>'ᵥ', 'x'=>'ₓ', 
                           '+'=>'₊', '-'=>'₋', '='=>'₌', '('=>'₍', ')'=>'₎',  # special characters
                           'β'=>'ᵦ', 'γ'=>'ᵧ', 'ρ'=>'ᵨ', 'ψ'=>'ᵩ', 'χ'=>'ᵪ',  # greek
                           # missing letter subscripts: b, c, d, f, g, q, w, y, z
                           )                                          
const SUPSCRIPT_MAP = Dict('1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', '5'=>'⁵',  # digits
                           '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '0'=>'⁰',
                           'a'=>'ᵃ', 'b'=>'ᵇ', 'c'=>'ᶜ', 'd'=>'ᵈ', 'e'=>'ᵉ', 
                           'f'=>'ᶠ', 'g'=>'ᵍ', 'h'=>'ʰ', 'i'=>'ⁱ', 'j'=>'ʲ',  # letters (only 'q' missing)
                           'k'=>'ᵏ', 'l'=>'ˡ', 'm'=>'ᵐ', 'n'=>'ⁿ', 'o'=>'ᵒ', 
                           'p'=>'ᵖ', 'r'=>'ʳ', 's'=>'ˢ', 't'=>'ᵗ', 'u'=>'ᵘ', 
                           'v'=>'ᵛ', 'w'=>'ʷ', 'x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ',
                           '+'=>'⁺', '-'=>'⁻', '='=>'⁼', '('=>'⁽', ')'=>'⁾',  # special characters
                           'α'=>'ᵅ', 'β'=>'ᵝ', 'γ'=>'ᵞ', 'δ'=>'ᵟ', 'ε'=>'ᵋ',  # greek
                           'θ'=>'ᶿ', 'ι'=>'ᶥ', 'φ'=>'ᶲ', 'ψ'=>'ᵠ', 'χ'=>'ᵡ',
                           # missing letter superscripts: q
                           )                                          
const SUBSCRIPT_MAP_REVERSE = Dict(v=>k for (k,v) in SUBSCRIPT_MAP)
const SUPSCRIPT_MAP_REVERSE = Dict(v=>k for (k,v) in SUPSCRIPT_MAP)

subscriptify(str::AbstractString) = map(subscriptify, str)
function subscriptify(c::Char)
    if c ∈ keys(SUBSCRIPT_MAP)
        return SUBSCRIPT_MAP[c]
    else
        return c
    end
end

supscriptify(str::AbstractString) = map(supscriptify, str)
function supscriptify(c::Char) 
    if c ∈ keys(SUPSCRIPT_MAP)
        return SUPSCRIPT_MAP[c]
    else
        return c
    end
end

function formatirreplabel(str::AbstractString)
    buf = IOBuffer()
    for c in str
        if c ∈ ['+','-']
            write(buf, supscriptify(c))
        elseif isdigit(c)
            write(buf, subscriptify(c))
        else
            write(buf, c)
        end
    end
    return String(take!(buf))
end


normalizesubsup(str::AbstractString) = map(normalizesubsup, str)
function normalizesubsup(c::Char)
    if c ∈ keys(SUBSCRIPT_MAP_REVERSE)
        return SUBSCRIPT_MAP_REVERSE[c]
    elseif c ∈ keys(SUPSCRIPT_MAP_REVERSE)
        return SUPSCRIPT_MAP_REVERSE[c]
    else 
        return c
    end
end

issubdigit(c::AbstractChar) = (c >= '₀') & (c <= '₉')

function unicode_frac(x::Number)
    xabs=abs(x)
    if     xabs == 0;   return "0" # early termination for common case & avoids undesirable sign for -0.0
    elseif xabs ≈ 1/2;  xstr = "½"
    elseif xabs ≈ 1/3;  xstr = "⅓"
    elseif xabs ≈ 2/3;  xstr = "⅔"
    elseif xabs ≈ 1/4;  xstr = "¼"
    elseif xabs ≈ 3/4;  xstr = "¾"
    elseif xabs ≈ 1/5;  xstr = "⅕"
    elseif xabs ≈ 2/5;  xstr = "⅖"
    elseif xabs ≈ 3/5;  xstr = "⅗"
    elseif xabs ≈ 4/5;  xstr = "⅘"
    elseif xabs ≈ 1/6;  xstr = "⅙"
    elseif xabs ≈ 5/6;  xstr = "⅚"
    elseif xabs ≈ 1/7;  xstr = "⅐"
    elseif xabs ≈ 1/8;  xstr = "⅛"
    elseif xabs ≈ 3/8;  xstr = "⅜"
    elseif xabs ≈ 5/8;  xstr = "⅝"
    elseif xabs ≈ 7/8;  xstr = "⅞"
    elseif xabs ≈ 1/9;  xstr = "⅑"
    elseif xabs ≈ 1/10; xstr = "⅒"
    else xstr = string(xabs) # return a conventional string representation
    end
    return signbit(x) ? "-"*xstr : xstr
end

const roman2greek_dict = Dict("LD"=>"Λ", "DT"=>"Δ", "SM"=>"Σ", "GM"=>"Γ", "GP"=>"Ω")
                              #"LE"=>"Λ′", "DU"=>"Δ′", "SN"=>"Σ′",  # These are the awkwardly annoted analogues of the pairs (Z,ZA), (W,WA) etc. 
                              #"ZA"=>"Z′", "WA"=>"W′")              # They "match" a simpler k-vector, by reducing their second character by one,
                                                                    # alphabetically (e.g. LE => LD = Λ). The primed notation is our own (actually,
                                                                    # it is also used in B&C, e.g. p. 412).
function roman2greek(label::String)
    idx = findfirst(!isletter, label)
    if idx !== nothing
        front=label[firstindex(label):prevind(label,idx)]
        if front ∈ keys(roman2greek_dict)
            return roman2greek_dict[front]*label[idx:lastindex(label)]
        end
    end
    return label
end


function printboxchar(io, i, N)
    if i == 1
        print(io, "╭") #┌
    elseif i == N
        print(io, "╰") #└
    else
        print(io, "│")
    end
end



function readuntil(io::IO, delim::F; keep::Bool=false) where F<:Function
    buf = IOBuffer()
    while !eof(io)
        c = read(io, Char)
        if delim(c)
            keep && write(buf, c)
            break
        end
        write(buf, c)
    end
    return String(take!(buf))
end


"""
    compact_print_matrix(io)

Canibalized and adapted from Base.print_matrix, specifically to allow a `prerow` input.

Should never be used for printing very large matrices, as it will not wrap or abbreviate
rows/columns.
"""
function compact_print_matrix(io, X::Matrix, prerow, elformat=(x->round(x,digits=2)), sep="")
    X_formatted = elformat.(X) # allocates; can't be bothered... (could be fixed using MappedArrays)
    screenheight = screenwidth = typemax(Int)
    rowsA, colsA = UnitRange(axes(X,1)), UnitRange(axes(X,2))

    A = Base.alignment(io, X_formatted, rowsA, colsA, screenwidth, screenwidth, length(sep))
    for i in rowsA
        i != first(rowsA) && print(io, prerow)
        print(io, i == first(rowsA) ? '┌' : (i == last(rowsA) ? '└' : '│'), ' ')
        Base.print_matrix_row(IOContext(io, :compact=>true), X_formatted, A, i, colsA, sep)
        print(io, ' ', i == first(rowsA) ? '┐' : (i == last(rowsA) ? '┘' : '│'))
        if i != last(rowsA); println(io); end
    end
end

# === misc functionality ===

"""
    isapproxin(x, itr) --> Bool

Determine whether `x` ∈ `itr` with approximate equality.
"""
isapproxin(x, itr, optargs...; kwargs...) = any(y -> isapprox(y, x, optargs...; kwargs...), itr)


"""
    uniquetol(a; kwargs)

Computes approximate-equality unique with tolerance specifiable
via keyword arguments `kwargs` in O(n²) runtime.

Copied from https://github.com/JuliaLang/julia/issues/19147#issuecomment-256981994
"""
function uniquetol(A::AbstractArray{T}; kwargs...) where T
    S = Vector{T}()
    for a in A
         if !any(s -> isapprox(s, a; kwargs...), S)
             push!(S, a)
         end
    end
    return S
end


"""
    get_kvpath(kvs::T, Ninterp::Integer) 
        where T<:AbstractVector{<:AbstractVector{<:Real}} --> Vector{Vector{Float64}}

Compute an interpolated k-path between discrete k-points given in `kvs` (a vector of
vectors of `Real`s), so that the interpolated path has `Ninterp` points in total.

Note that, in general, it is not possible to do this so that all points are equidistant; 
but points are equidistant in-between the initial discrete points provided in `kvs`.
"""
function get_kvpath(kvs::AbstractVector{<:AbstractVector{<:Real}}, Ninterp::Integer)
    Nkpairs = length(kvs)-1
    dists = Vector{Float64}(undef, Nkpairs)
    @inbounds for i in Base.OneTo(Nkpairs)
        dists[i] = norm(kvs[i] .- kvs[i+1])
    end
    mindist = mean(dists)

    kvpath = [float.(kvs[1])]
    @inbounds for i in  Base.OneTo(Nkpairs)
        # try to maintain an even distribution of k-points along path
        Ninterp_i = round(Int64, dists[i]./mindist*Ninterp)
        # new k-points
        newkvs = range(kvs[i],kvs[i+1],length=Ninterp_i)
        # append new kvecs to kpath
        append!(kvpath, (@view newkvs[2:end]))
    end
    return kvpath
end


"""
    ImmutableDict(ps::Pair...)

Construct an `ImmutableDict` from any number of `Pair`s; a convenience function that extends
`Base.ImmutableDict` which otherwise only allows construction by iteration.

# TODO: Will be included in Base at some point? Maybe v1.4?
"""
function ImmutableDict(ps::Pair{K,V}...) where {K,V}
    d = ImmutableDict{K,V}()
    for p in ps # construct iteratively (linked list)
        d = ImmutableDict(d, p)
    end
    return d
end
ImmutableDict(ps::Pair...) = ImmutableDict(ps)


# Distributions.jl provides a nice Uniform distribution type - but it is not worth adding
# the high compilation-time of Distributions (~7 s) just for that functionality, so we just
# copy a subset of the methods and the struct here, distinguishing the struct-pirating by an 
# underscore. We only copy the scalar rand(..) methods (i.e. no array generators).
struct _Uniform{T<:Real}
    a::T # low
    b::T # high (unchecked...)
end
_Uniform(a::Real, b::Real) = _Uniform(promote(a, b)...)
_Uniform(a::Integer, b::Integer) = _Uniform(float(a), float(b))
rand(u::_Uniform) = rand(Random.GLOBAL_RNG, u)
rand(rng::Random.AbstractRNG, u::_Uniform) = u.a + (u.b - u.a) * rand(rng)