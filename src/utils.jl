""" 
    parsefraction(str)

    Parse a string, allowing fraction inputs (e.g. "1/2"), return as Float64
"""
function parsefraction(str)
    slashidx = findfirst(x->x=='/',str)
    if isnothing(slashidx)
        return parse(Float64, str)
    else
        num=str[1:prevind(str, slashidx)]
        den=str[nextind(str, slashidx):end]
        return parse(Float64, num)/parse(Float64, den)
    end
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