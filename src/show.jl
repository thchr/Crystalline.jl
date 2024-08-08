# --- DirectBasis ---
function show(io::IO, ::MIME"text/plain", Vs::DirectBasis)
    # cannot use for ReciprocalBasis at the moment (see TODO in `crystalsystem`)
    print(io, typeof(Vs))
    print(io, " ($(crystalsystem(Vs))):")
    for (i,V) in enumerate(Vs)
        print(io, "\n   ", V)
    end
end

# --- SymOperation ---
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

# --- MultTable ---
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

# --- AbstractVec ---
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

# --- AbstractGroup ---
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

# --- LGIrrep & PGIrrep ---
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
        #Base.print_matrix(IOContext(io, :compact=>true), ir, indent*(i == Nₒₚ ? " ╰" : " │")*"    ")
        print(io, indent, " │     ")
        prettyprint_irrep_matrix(io, ir, i, indent*" │     ")
        if i < Nₒₚ; println(io, '\n', indent, " │"); end
    end
    print(io, "\n", indent, " └", boxdelims)
end

function prettyprint_header(io::IO, irlab::AbstractString, nboxdelims::Integer=45)
    println(io, irlab, " ─┬", repeat("─", nboxdelims))
end

# --- IrrepCollection ---
function summary(io::IO, c::IrrepCollection{T}) where T
    print(io, length(c), "-element IrrepCollection{", T, "}")
end
function show(io::IO, ::MIME"text/plain", c::IrrepCollection)
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
function show(io::IO, c::IrrepCollection)
    show(io, c.irs)
    g = group(first(c))
    if position(g) !== nothing
        printstyled(io, " (", fullpositionlabel(g), ")"; color=:light_black)
    end
end


# --- CharacterTable ---
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


# --- BandRep ---
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


# --- BandRepSet ---
function show(io::IO, ::MIME"text/plain", BRS::BandRepSet)
    Nⁱʳʳ = length(irreplabels(BRS))
    Nᵉᵇʳ = length(BRS)

    # print a "title" line and the irrep labels
    println(io, "BandRepSet (⋕", num(BRS), "): ",
                length(BRS), " BandReps, ",
                "sampling ", Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(BRS) ? "½" : "1", " ",
                BRS.timereversal ? "w/" : "w/o", " TR)")

    # print band representations as table
    k_idx = (i) -> findfirst(==(klabel(irreplabels(BRS)[i])), klabels(BRS)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_μ   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")
    pretty_table(io, 
        # table contents
        matrix(BRS; includedim=true);
        # row/column names
        row_labels = vcat(irreplabels(BRS), "μ"),
        header = (position.(BRS), chop.(label.(BRS), tail=2)), # remove repetitive "↑G" postfix
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
    print(io, "  KVecs (", hasnonmax(BRS) ? "incl. non-maximal" : "maximal only", "): ")
    join(io, klabels(BRS), ", ")
    
    # EARLIER MANUAL LAYOUTS: DISCONTINUED

    #=
    # === EBRS-by-column layout ===
    # prep-work to figure out how many bandreps we can write to the io
    cols_avail = displaysize(io)[2]   # available cols in io (cannot write to all of it; subtract 2)
    indent     = 3
    μ_maxdigs  = maximum(ndigits∘dim, BRS) # implicitly assuming this to be ≤2...
    maxcols_irr   = maximum(length, irreplabels(BRS))
    cols_brlabs = length.(label.(BRS))
    cumcols_brlabs = cumsum(cols_brlabs .+ 1)
    cols_requi  = cumcols_brlabs .+ (indent + maxcols_irr + 2)
    @show cols_requi
    room_for = findlast(≤(cols_avail), cols_requi)
    abbreviate = room_for < length(BRS) ? true : false
    rowend = abbreviate ? '…' : '║'

    # print EBR header
    print(io, ' '^(indent+maxcols_irr+1), '║')
    for idxᵇʳ in 1:room_for
        br = BRS[idxᵇʳ]
        print(io, ' ', chop(label(br), tail=2),' ', idxᵇʳ ≠ room_for ? '│' : rowend)
    end
    println(io)

    # print irrep content, line by line
    for idxⁱʳʳ in 1:Nⁱʳʳ
        irlab = irreplabels(BRS)[idxⁱʳʳ]
        print(io, ' '^indent, irlab, ' '^(maxcols_irr + 1 - length(irlab)), '║')
        for idxᵇʳ in 1:room_for
            x = BRS[idxᵇʳ][idxⁱʳʳ]
            addspace = div(cols_brlabs[idxᵇʳ], 2)
            print(io, ' '^addspace)
            iszero(x) ? print(io, '·') : print(io, x)
            print(io, ' '^(cols_brlabs[idxᵇʳ] - addspace-1))
            print(io, idxᵇʳ ≠ room_for ? '│' : rowend)
        end
        println(io)
    end

    # print band filling
    print(io, ' '^indent, 'μ', ' '^maxcols_irr, '║')
    for idxᵇʳ in 1:room_for
        μ = dim(BRS[idxᵇʳ])
        addspace = div(cols_brlabs[idxᵇʳ], 2)+1
        print(io, ' '^(addspace-ndigits(μ)), μ)
        print(io, ' '^(cols_brlabs[idxᵇʳ] - addspace))
        print(io, idxᵇʳ ≠ room_for ? '│' : rowend)
    end
    =#

    #=
    # === EBRs-by-rows layout ===
    μ_maxdigs = maximum(ndigits∘dim, BRS)
    cols_brlab = maximum(x->length(label(x)), BRS)+1
    cols_irstart = cols_brlab+4
    cols_avail = displaysize(io)[2]-2                                 # available cols in io (cannot write to all of it; subtract 2)
    cols_requi = sum(x->length(x)+3, irreplabels(BRS))+cols_irstart+μ_maxdigs+3 # required cols for irrep labels & band reps
    if cols_requi > cols_avail
        cols_toomany    = ceil(Int, (cols_requi-cols_avail)/2) + 2  # +2 is to make room for '  …  ' extender
        cols_midpoint   = div(cols_requi-cols_irstart,2)+cols_irstart
        cols_skipmin    = cols_midpoint - cols_toomany
        cols_skipmax    = cols_midpoint + cols_toomany
        cols_eachstart  = [0; cumsum(length.(irreplabels(BRS)).+3)].+cols_irstart
        iridx_skiprange = [idx for (idx, col_pos) in enumerate(cols_eachstart) if cols_skipmin ≤ col_pos ≤ cols_skipmax]
        abbreviate = true
    else
        abbreviate = false
    end

    print(io, " "^(cols_irstart-1),'║'); # align with spaces
    for (iridx,lab) in enumerate(irreplabels(BRS)) # irrep labels
        if abbreviate && iridx ∈ iridx_skiprange
            if iridx == first(iridx_skiprange)
                print(io, "\b  …  ")
            end
        else
            print(io, ' ', lab, " │")
        end
    end
    println(io, ' '^μ_maxdigs, "μ", " ║") # band-filling column header
    # print each bandrep
    for (bridx,BR) in enumerate(BRS)
        μ = dim(BR)
        print(io, "   ", label(BR),                      # bandrep label
                  " "^(cols_brlab-length(label(BR))), '║')
        for (iridx,x) in enumerate(BR) # iterate over vector representation of band rep
            if abbreviate && iridx ∈ iridx_skiprange
                if iridx == first(iridx_skiprange)
                    print(io, mod(bridx,4) == 0 ? "\b  …  " : "\b     ")
                end
            else
                print(io, "  ")
                !iszero(x) ? print(io, x) : print(io, '·')
                print(io, " "^(length(irreplabels(BRS)[iridx])-1), '│') # assumes we will never have ndigit(x) != 1
            end
        end
        
        print(io, ' '^(1+μ_maxdigs-ndigits(μ)), μ, " ║") # band-filling
        if bridx != length(BRS); println(io); end
    end
    =#
end