# ---------------------------------------------------------------------------------------- #
# SymOperation
function show(io::IO, ::MIME"text/plain", op::AbstractOperation{D}) where D
    opseitz, opxyzt = seitz(op), xyzt(op)
    print(io, opseitz)
    
    # don't print triplet & matrix format if the IOContext is :compact=>true
    if get(io, :compact, false)
        return nothing
    end

    # --- print triplet expression ---
    printstyled(io, " ", repeat('─',max(38-length(opseitz)-length(opxyzt), 1)),
                    " (", opxyzt, ")"; color=:light_black)
    println(io)

    # --- print matrix ---
    # info that is needed before we start writing by column
    τstrs = fractionify.(translation(op), false)
    Nsepτ = maximum(length, τstrs)
    firstcol_hasnegative = any(_has_negative_sign_and_isnonzero, @view matrix(op)[:,1])
    for i in 1:D
        printstyled(io, " ", i == 1 ? '┌' : (i == D ? '└' : '│'), color=:light_black) # open brace char
        for j in 1:D
            c = matrix(op)[i,j]
            # assume and exploit that a valid symop (in the lattice basis only!) never has an 
            # entry that is more than two characters long (namely, -1) in its rotation parts
            sep = repeat(' ', 1 + (j ≠ 1 || firstcol_hasnegative) - _has_negative_sign_and_isnonzero(c))
            if isinteger(c)
                cᴵ = convert(Int, matrix(op)[i,j])
                printstyled(io, sep, cᴵ, color=:light_black)
            else
                # just use the same sep even if the symop is specified in a nonstandard basis (e.g.
                # cartesian); probably isn't a good general solution, but good enough for now
                printstyled(io, sep, round(c; digits=4), color=:light_black)
            end
        end
        printstyled(io, " ", i == 1 ? "╷" : (i == D ? "╵" : "┆"), " ", repeat(' ', Nsepτ-length(τstrs[i])), τstrs[i], " ", color=:light_black)
        printstyled(io, i == 1 ? '┐' : (i == D ? '┘' : '│'), color=:light_black) # close brace char
        op isa MSymOperation && i == 1 && timereversal(op) && print(io, '′')
        i ≠ D && println(io)
    end
    return nothing
end
_has_negative_sign_and_isnonzero(x) = !iszero(x) && signbit(x)
# print vectors of `SymOperation`s compactly
show(io::IO, op::AbstractOperation) = print(io, seitz(op))

# ---------------------------------------------------------------------------------------- #
# MultTable

function show(io::IO, ::MIME"text/plain", mt::MultTable)
    summary(io, mt)
    println(io, ":")
    seitz_ops = seitz.(mt.operations)
    pretty_table(io,
        getindex.(Ref(seitz_ops), mt.table);
        row_labels = seitz_ops,
        header = seitz_ops,
        vlines = [1,],
        hlines = [:begin, 1, :end]
        )
    return nothing
end

# ---------------------------------------------------------------------------------------- #
# AbstractVec

function show(io::IO, ::MIME"text/plain", v::AbstractVec)
    cnst, free = parts(v)
    print(io, '[')
    if isspecial(v)
        for i in eachindex(cnst) 
            coord = cnst[i] == -0.0 ? 0.0 : cnst[i] # normalize -0.0 to 0.0
            prettyprint_scalar(io, coord)
            # prepare for next coordinate/termination
            i == length(cnst) ? print(io, ']') : print(io, ", ")
        end
    else
        for i in eachindex(cnst)
            # constant/fixed parts
            if !iszero(cnst[i]) || iszero(@view free[i,:]) # don't print zero, if it adds unto anything nonzero
                coord = cnst[i] == -0.0 ? 0.0 : cnst[i] # normalize -0.0 to 0.0
                prettyprint_scalar(io, coord)
            end
            # free-parameter parts
            for j in eachindex(cnst) 
                if !iszero(free[i,j])
                    sgn = signaschar(free[i,j])
                    if !(iszero(cnst[i]) && sgn=='+' && iszero(free[i,1:j-1])) # don't print '+' if nothing precedes it
                        print(io, sgn)
                    end
                    if abs(free[i,j]) != oneunit(eltype(free)) # don't print prefactors of 1
                        prettyprint_scalar(io, abs(free[i,j]))
                    end
                    print(io, j==1 ? 'α' : (j == 2 ? 'β' : 'γ'))
                end
            end
            # prepare for next coordinate/termination
            i == length(cnst) ? print(io, ']') : print(io, ", ")
        end
    end
    return
end
# print arrays of `AbstractVec`s compactly
show(io::IO, v::AbstractVec) = show(io, MIME"text/plain"(), v)

function prettyprint_scalar(io, v::Real)
    if isinteger(v)
        print(io, Int(v))
    else
        # print all fractions divisible by 2, ..., 10 as fractions, and everything else
        # as decimal
        rv = rationalize(Int, v; tol=1e-2)
        if isapprox(v, rv; atol=DEFAULT_ATOL)
            print(io, rv.num, "/", rv.den)
        else
            print(io, v)
        end
    end
end

function show(io::IO, ::MIME"text/plain", wp::WyckoffPosition)
    print(io, wp.mult, wp.letter, ": ") # TODO: This is `=` elsewhere; change?
    show(io, MIME"text/plain"(), parent(wp))
end

# ---------------------------------------------------------------------------------------- #
# AbstractGroup

function summary(io::IO, g::AbstractGroup)
    print(io, typeof(g))
    _print_group_descriptor(io, g; prefix=" ")
    print(io, " with ", order(g), " operations")
end
function show(io::IO, ::MIME"text/plain", g::AbstractGroup)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    summary(io, g)
    println(io, ':')
    for (i, op) in enumerate(g)
        print(io, ' ')
        show(io, MIME"text/plain"(), op)
        if i < order(g); println(io); end
    end
end
function show(io::IO, g::AbstractGroup)
    print(io, '[')
    join(io, g, ", ")
    print(io, ']')
end
function show(io::IO, g::Union{LittleGroup, SiteGroup})
    print(io, '[')
    join(io, g, ", ")
    print(io, ']')
    printstyled(io, " (", position(g), ")", color=:light_black)
end

function _print_group_descriptor(io::IO, g::AbstractGroup; prefix::AbstractString="")
    print(io, prefix)
    g isa GenericGroup && return nothing
    print(io, "⋕")
    join(io, num(g), '.') # this slightly odd approach to treat magnetic groups also
    print(io, " (", label(g), ")")
    if position(g) !== nothing
        print(io, " at ")
        print(io, fullpositionlabel(g))
    end
    return nothing
end
function _group_descriptor(g; prefix::AbstractString="")
    return sprint( (io, _g) -> _print_group_descriptor(io, _g; prefix), g)
end

# ---------------------------------------------------------------------------------------- #
# AbstractIrrep

function show(io::IO, ::MIME"text/plain", ir::AbstractIrrep)
    irlab = label(ir)
    lablen = length(irlab)
    nindent = lablen+1
    prettyprint_header(io, irlab)
    prettyprint_irrep_matrices(io, ir, nindent)
end
function show(io::IO, ir::AbstractIrrep)
    print(io, label(ir))
end

# ... utilities to print PGIrreps and LGIrreps
function prettyprint_group_header(io::IO, g::AbstractGroup)
    print(io, "⋕", num(g), " (", iuc(g), ")")
    if g isa LittleGroup
        print(io, " at " , klabel(g), " = ")
        show(io, MIME"text/plain"(), position(g))
    end
    println(io)
end

function prettyprint_scalar_or_matrix(io::IO, printP::AbstractMatrix, prefix::AbstractString,
                                      ϕabc_contrib::Bool=false, digits::Int=4)
    if size(printP) == (1,1) # scalar case
        v = @inbounds printP[1]
        prettyprint_irrep_scalars(io, v, ϕabc_contrib; digits)

    else # matrix case
        formatter(x) = _stringify_characters(x; digits)
        # FIXME: not very optimal; e.g. makes a whole copy and doesn't handle displaysize
        compact_print_matrix(io, printP, prefix, formatter)
    end
end

function prettyprint_irrep_scalars(
        io::IO, v::Number, ϕabc_contrib::Bool=false;
        atol::Real=DEFAULT_ATOL, digits::Int=4
    )

    if norm(v) < atol
        print(io, 0)
    elseif isapprox(v, real(v), atol=atol)     # real scalar
        if ϕabc_contrib && isapprox(abs(real(v)), 1.0, atol=atol)
            signbit(real(v)) && print(io, '-')
        else
            print(io, _stringify_characters(real(v); digits))
        end
    elseif isapprox(v, imag(v)*im, atol=atol)   # imaginary scalar
        if ϕabc_contrib && isapprox(abs(imag(v)), 1.0, atol=atol)
            signbit(imag(v)) && print(io, '-')
        else
            print(io, _stringify_characters(imag(v); digits))
        end
        print(io, "i")
    else                                        # complex scalar (print as polar)
        vρ, vθ = abs(v), angle(v)
        vθ /= π
        isapprox(vρ, 1.0, atol=atol) || print(io, _stringify_characters(vρ; digits))
        print(io, "exp(") 
        if isapprox(abs(vθ), 1.0, atol=atol)
            signbit(vθ) && print(io, '-')
        else
            print(io, _stringify_characters(vθ; digits))
        end
        print(io, "iπ)")
        #print(io, ϕabc_contrib ? "(" : "", v, ϕabc_contrib ? ")" : "")
    end
end

function prettyprint_irrep_matrix(
        io::IO, lgir::LGIrrep, i::Integer, prefix::AbstractString;
        digits::Int=4
    )
    # unpack
    k₀, kabc = parts(position(group(lgir)))
    P = lgir.matrices[i]
    τ = lgir.translations[i]

    # phase contributions
    ϕ₀ = dot(k₀, τ)                                   # constant phase
    ϕabc = [dot(kabcⱼ, τ) for kabcⱼ in eachcol(kabc)] # variable phase
    ϕabc_contrib = norm(ϕabc) > sqrt(dim(lgir))*DEFAULT_ATOL

    # print the constant part of the irrep that is independent of α,β,γ
    printP = abs(ϕ₀) < DEFAULT_ATOL ? P : cis(2π*ϕ₀)*P # avoids copy if ϕ₀≈0; copies otherwise
    prettyprint_scalar_or_matrix(io, printP, prefix, ϕabc_contrib)

    # print the variable phase part that depends on the free parameters α,β,γ 
    if ϕabc_contrib
        nnzabc = count(c->abs(c)>DEFAULT_ATOL, ϕabc)
        print(io, "exp")
        if nnzabc == 1
            print(io, "(")
            i = findfirst(c->abs(c)>DEFAULT_ATOL, ϕabc)
            c = ϕabc[i]
            signbit(c) && print(io, "-")
            if !(abs(c) ≈ 0.5) # do not print if multiplicative factor is 1
                print(io, _stringify_characters(abs(2c); digits))
            end

            print(io, "iπ", 'ΰ'+i, ")") # prints 'α', 'β', and 'γ' for i = 1, 2, and 3, respectively ('ΰ'='α'-1)

        else
            print(io, "[iπ(")
            first_nzidx = true
            for (i,c) in enumerate(ϕabc)
                if abs(c) > DEFAULT_ATOL
                    if first_nzidx 
                        signbit(c) && print(io, '-')
                        first_nzidx = false
                    else
                        print(io, signaschar(c))
                    end
                    if !(abs(c) ≈ 0.5) # do not print if multiplicative factor is 1
                        print(io, _stringify_characters(abs(2c); digits))
                    end
                    print(io, 'ΰ'+i) # prints 'α', 'β', and 'γ' for i = 1, 2, and 3, respectively ('ΰ'='α'-1)
                end
            end
            print(io, ")]")
        end
    end
end

function prettyprint_irrep_matrix(
        io::IO, ir::Union{<:PGIrrep, <:SiteIrrep}, i::Integer, prefix::AbstractString
    )
    P = ir.matrices[i]
    prettyprint_scalar_or_matrix(io, P, prefix, false)
end

function prettyprint_irrep_matrices(
        io::IO, ir::Union{<:LGIrrep, <:PGIrrep, <:SiteIrrep}, nindent::Integer,
        nboxdelims::Integer=45
    )
    indent = repeat(" ", nindent)
    boxdelims = repeat("─", nboxdelims)
    linelen = nboxdelims + 4 + nindent
    Nₒₚ = order(ir)
    for (i,op) in enumerate(operations(ir))
        print(io, indent, " ├─ ")
        opseitz, opxyzt  = seitz(op), xyzt(op)
        printstyled(io, opseitz, ": ", 
                        repeat("─", linelen-11-nindent-length(opseitz)-length(opxyzt)),
                        " (", opxyzt, ")\n"; color=:light_black)
        print(io, indent, " │     ")
        prettyprint_irrep_matrix(io, ir, i, indent*" │     ")
        if i < Nₒₚ; println(io, '\n', indent, " │"); end
    end
    print(io, "\n", indent, " └", boxdelims)
end

function prettyprint_header(io::IO, irlab::AbstractString, nboxdelims::Integer=45)
    println(io, irlab, " ─┬", repeat("─", nboxdelims))
end

# ---------------------------------------------------------------------------------------- #
# Collection{<:AbstractIrrep}

function summary(io::IO, c::Collection{T}) where T <: AbstractIrrep
    print(io, length(c), "-element Collection{", T, "}")
end
function show(io::IO, ::MIME"text/plain", c::Collection{T}) where T <: AbstractIrrep
    summary(io, c)
    isassigned(c, firstindex(c)) && _print_group_descriptor(io, group(first(c)); prefix=" for ")
    println(io, ":")
    for i in eachindex(c)
        if isassigned(c, i)
            show(io, MIME"text/plain"(), c[i])
        else
            print(io, " #undef")
        end
        i ≠ length(c) && println(io)
    end
end
function show(io::IO, c::Collection{T}) where T <: AbstractIrrep
    show(io, c.vs)
    g = group(first(c))
    if position(g) !== nothing
        printstyled(io, " (", fullpositionlabel(g), ")"; color=:light_black)
    end
end

# ---------------------------------------------------------------------------------------- #
# CharacterTable

function show(io::IO, ::MIME"text/plain", ct::AbstractCharacterTable)
    chars = matrix(ct)
    chars_formatted = _stringify_characters.(chars; digits=4)

    ops = operations(ct)
    println(io, typeof(ct), " for ", tag(ct), ":") # type name and space group/k-point tags
    pretty_table(io,
        chars_formatted;
        # row/column names
        row_labels = seitz.(ops), # seitz labels
        header = labels(ct),     # irrep labels
        tf = tf_unicode,
        vlines = [1,], hlines = [:begin, 1, :end]
        )

    if ct isa ClassCharacterTable
        _print_class_representatives(io, ct)
    end
end

function _print_class_representatives(io::IO, ct::ClassCharacterTable)
    maxlen = maximum(length∘seitz∘first, classes(ct))
    print(io, "Class representatives:")
    for class in classes(ct)
        print(io, "\n ")
        for (i, op) in enumerate(class)#
            op_str = seitz(op)
            if i == 1
                printstyled(io, op_str; bold=true)
                length(class) ≠ 1 && print(io, " "^(maxlen-length(op_str)), " : ")
            else
                printstyled(io, op_str; color=:light_black)
                i ≠ length(class) && printstyled(io, ", "; color=:light_black)
            end
        end
    end
end

function _stringify_characters(c::Number; digits::Int=4)
    c′ = round(c; digits)
    cr, ci = reim(c′)
    if iszero(ci)     # real
        isinteger(cr) && return string(Int(cr))
        return string(cr)

    elseif iszero(cr) # imaginary
        isinteger(ci) && return string(Int(ci))*"im"
        return string(ci)*"im"

    else              # complex
        (isinteger(cr) && isinteger(ci)) && return _complex_as_compact_string(Complex{Int}(cr,ci))
        return _complex_as_compact_string(c′)
    end
end
function _complex_as_compact_string(c::Complex) # usual string(::Complex) has spaces; avoid that
    io = IOBuffer()
    print(io, real(c), signaschar(imag(c)), abs(imag(c)), "im")
    return String(take!(io))
end

# ---------------------------------------------------------------------------------------- #
# BandRep

function prettyprint_symmetryvector(
            io::IO, 
            irvec::AbstractVector{<:Real},
            irlabs::Vector{String};
            braces::Bool=true)

    Nⁱʳʳ  = length(irlabs)
    Nⁱʳʳ′ = length(irvec) 
    if !(Nⁱʳʳ′ == Nⁱʳʳ || Nⁱʳʳ′ == Nⁱʳʳ+1)
        # we allow irvec to exceed the dimension of irlabs by 1, in case it includes dim(BR)
        throw(DimensionMismatch("irvec and irlabs must have matching dimensions"))
    end
    braces && print(io, '[')

    first_nz   = true
    group_klab = klabel(first(irlabs))
    for idx in 1:Nⁱʳʳ
        # shared prepwork
        irlab = irlabs[idx]
        klab  = klabel(irlab)

        if klab ≠ group_klab    # check if this irrep belongs to a new k-label group
            first_nz = true
            group_klab = klab
            print(io, ", ") 
        end

        c  = irvec[idx]         # coefficient of current term
        c == 0 && continue      # early stop if the coefficient is zero
        
        absc = abs(c)
        if first_nz             # first nonzero entry in a k-label group
            if abs(c) == 1
                c == -1 && print(io, '-')
            else
                print(io, c)
            end
            first_nz = false
        else                    # entry that is added/subtracted from another irrep
            print(io, signaschar(c))
            if absc ≠ 1
                print(io, absc)
            end
        end
        print(io, irlab)
    end
    braces && print(io, ']')
end
function symvec2string(irvec::AbstractVector{<:Real}, irlabs::Vector{String}; 
            braces::Bool=true)
    io = IOBuffer()
    prettyprint_symmetryvector(io, irvec, irlabs; braces=braces)
    return String(take!(io))
end

summary(io::IO, BR::BandRep) = print(io, dim(BR), "-band BandRep (", label(BR), " at ", position(BR), ")")
function show(io::IO, ::MIME"text/plain", BR::BandRep)
    summary(io, BR)
    print(io, ":\n ")
    prettyprint_symmetryvector(io, BR, irreplabels(BR))
end
function show(io::IO, BR::BandRep)
    prettyprint_symmetryvector(io, BR, irreplabels(BR))
end

# ---------------------------------------------------------------------------------------- #
# BandRepSet

function show(io::IO, ::MIME"text/plain", brs::BandRepSet)
    Nⁱʳʳ = length(irreplabels(brs))
    Nᵉᵇʳ = length(brs)

    # print a "title" line and the irrep labels
    println(io, "BandRepSet (⋕", num(brs), "): ",
                length(brs), " BandReps, ",
                "sampling ", Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(brs) ? "½" : "1", " ",
                brs.timereversal ? "w/" : "w/o", " TR)")

    # print band representations as table
    k_idx = (i) -> findfirst(==(klabel(irreplabels(brs)[i])), klabels(brs)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_μ   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")
    pretty_table(io, 
        # table contents
        stack(brs);
        # row/column names
        row_labels = vcat(irreplabels(brs), "μ"),
        header = (position.(brs), chop.(label.(brs), tail=2)), # remove repetitive "↑G" postfix
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,], hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_label_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_μ), 
        header_crayon = crayon"bold"
        # TODO: Would be nice to highlight the `row_labels` in a style matching the contents,
        #       but not possible atm (https://github.com/ronisbr/PrettyTables.jl/issues/122)
        )

    # print k-vec labels
    print(io, "  KVecs: ")
    join(io, klabels(brs), ", ")
end

# ---------------------------------------------------------------------------------------- #
# SymmetryVector

function Base.show(io :: IO, ::MIME"text/plain", n :: SymmetryVector)
    print(io, length(n)-1, "-irrep ", typeof(n), ":\n ")
    show(io, n)
end
function Base.show(io :: IO, n :: SymmetryVector)
    print(io, "[")
    for (i, (mults_k, lgirs_k)) in enumerate(zip(multiplicities(n), irreps(n)))
        str = if !iszero(mults_k)
            Crystalline.symvec2string(mults_k, label.(lgirs_k); braces=false)
        else # if there are no occupied irreps at the considered k-point print "0kᵢ"
            "0" * klabel(first(lgirs_k)) * "ᵢ"
        end
        printstyled(io, str; color=iseven(i) ? :normal : :light_blue)
        i ≠ length(multiplicities(n)) && print(io, ", ")
    end
    print(io, "]")
    printstyled(io, " (", occupation(n), " band", abs(occupation(n)) ≠ 1 ? "s" : "", ")"; 
                    color=:light_black)
end

# ---------------------------------------------------------------------------------------- #
# NewBandRep

function Base.show(io :: IO, ::MIME"text/plain", br :: NewBandRep)
    print(io, length(br.n)-1, "-irrep ", typeof(br), ":\n ")
    print(io, "(", )
    printstyled(io, label(position(br.siteir)); bold=true)
    print(io, "|")
    printstyled(io, label(br.siteir); bold=true)
    print(io, "): ")
    show(io, br.n)
end

function Base.show(io :: IO, br :: NewBandRep)
    print(io, "(", label(position(br.siteir)), "|", label(br.siteir), ")")
end

# ---------------------------------------------------------------------------------------- #
# Collection{<:NewBandRep}

function Base.show(io :: IO, ::MIME"text/plain", brs :: Collection{<:NewBandRep})
    irlabs = irreplabels(brs)
    Nⁱʳʳ = length(irlabs)

    # print a "summary" line
    print(io, length(brs), "-element ", typeof(brs), " for ⋕", num(brs))
    print(io, " (", iuc(num(brs), dim(brs)), ") ")
    print(io, "over ", Nⁱʳʳ, " irreps")
    print(io, " (spin-", first(brs).spinful ? "½" : "1", 
              " w/", first(brs).timereversal ? "" : "o", "TR):")
    println(io)

    # print band representations as table
    k_idx = (i) -> findfirst(==(klabel(irreplabels(brs)[i])), klabels(brs)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_μ   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")
    pretty_table(io, 
        # table contents
        stack(brs);
        # row/column names
        row_labels = vcat(irlabs, "μ"),
        header = (label.(position.(brs)), label.(getfield.(brs, :siteir))),
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,], hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_label_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_μ), 
        header_crayon = crayon"bold"
        # TODO: Would be nice to highlight the `row_labels` in a style matching the contents,
        #       but not possible atm (https://github.com/ronisbr/PrettyTables.jl/issues/122)
        )
end

# ---------------------------------------------------------------------------------------- #
# CompositeBandRep

function Base.show(io::IO, cbr::CompositeBandRep{D}) where D
    first = true
    for (j, c) in enumerate(cbr.coefs)
        iszero(c) && continue
        absc = abs(c)
        if first
            first = false
            c < 0 && print(io, "-")
        else
            print(io, " ", Crystalline.signaschar(c), " ")
        end
        if !isone(absc)
            if isinteger(absc)
                print(io, Int(absc))
            else
                print(io, "(", numerator(absc), "/", denominator(absc), ")×")
            end
        end
        print(io, cbr.brs[j])
    end
    first && print(io, "0")
end

function Base.show(io::IO, ::MIME"text/plain", cbr::CompositeBandRep{D}) where D
    println(io, length(irreplabels(cbr)), "-irrep ", typeof(cbr), ":")
    print(io, " ")
    show(io, cbr)
    μ = occupation(cbr)
    printstyled(io, " (", μ, " band", abs(μ) ≠ 1 ? "s" : "", ")", color=:light_black)
end