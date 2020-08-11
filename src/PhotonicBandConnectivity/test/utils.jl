function convert_irreplabel2latex(str::AbstractString)
    buf = IOBuffer()
    previous_was_digit = false
    for c in str
        if issubdigit(c)
            previous_was_digit || write(buf, "_{")
            write(buf, normalizesubsup(c)) 
            previous_was_digit = true
        else
            previous_was_digit && (write(buf, '}'); previous_was_digit = false)
            if c ∈ ('⁺', '⁻')
                write(buf, "^{", c=='⁺' ? '+' : '-', "}")
            else
                write(buf, latexifygreek(c))
            end
        end
    end
    previous_was_digit && write(buf, '}')
    return String(take!(buf))
end

const GREEK_SYMBOL_TO_LATEX_MAP = Dict{Char, String}(
    'Γ'=>"\\Gamma", 'Ω' => "\\Omega", 'Σ' => "\\Sigma", 'Λ' => "\\Lambda",
    'Δ'=>"\\Delta"
)

function latexifygreek(c::Char)
    if c ∈ keys(GREEK_SYMBOL_TO_LATEX_MAP)
        return GREEK_SYMBOL_TO_LATEX_MAP[c]
    else
        return string(c)
    end
end
