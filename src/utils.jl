""" 
    parsefraction(str::AbstractString)

Parse a string `str`, allowing fraction inputs (e.g. `"1/2"`), return as `Float64`.
"""
function parsefraction(str::AbstractString)
    slashidx = findfirst(x->x=='/',str)
    if isnothing(slashidx)
        return parse(Float64, str)
    else
        num=str[1:prevind(str, slashidx)]
        den=str[nextind(str, slashidx):end]
        return parse(Float64, num)/parse(Float64, den)
    end
end


"""
    isapproxin(x, itr) --> Bool

Determine whether `x` ∈ `itr` with approximate equality.
"""
isapproxin(x, itr, optargs...; kwargs...) = any(y -> isapprox(y, x, optargs...; kwargs...), itr)



# --- UNICODE FUNCTIONALITY ---
const SUBSCRIPT_MAP = Dict('1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', '5'=>'₅',  # digits
                           '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉', '0'=>'₀',
                           'a'=>'ₐ', 'e'=>'ₑ', 'h'=>'ₕ', 'i'=>'ᵢ', 'j'=>'ⱼ',  # letters (missing several)
                           'k'=>'ₖ', 'l'=>'ₗ', 'm'=>'ₘ', 'n'=>'ₙ', 'o'=>'ₒ', 
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
                           'p'=>'ᵖ', 'r'=>'ʳ', 's'=> 'ˢ', 't'=>'ᵗ', 'u'=>'ᵘ', 
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
    else xstr = @sprintf("%g",xabs) # return a conventional string representation
    end
    return signbit(x) ? "-"*xstr : xstr
end

const roman2greek_dict = Dict("LD"=>"Λ", "DT"=>"Δ", "SM"=>"Σ", "GM"=>"Γ", "GP"=>"Ω",
                              "LE"=>"Λ′", "DU"=>"Λ′", "SN"=>"Σ′")  # These are the awkwardly annoted analogues of the pairs (Z,ZA), (W,WA) etc. 
                                                                   # They "match" a simpler k-vector, by reducing their second character by one,
                                                                   # alphabetically (e.g. LE => LD = Λ). The primed notation is our own.
function roman2greek(label::String)
    idx = findfirst(!isletter, label)
    if idx != nothing
        front=label[firstindex(label):prevind(label,idx)]
        if length(front) == 2 # have to check length rather than just idx, in case of non-ascii input
            return roman2greek_dict[front]*label[idx:lastindex(label)]
        end
    end
    return label
end


function printboxchar(io, i, N)
    if i == 1
        print(io, "╭") #┌
    elseif i == N
        print(io, "╰") #┕
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