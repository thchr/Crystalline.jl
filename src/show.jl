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
function show(io::IO, ::MIME"text/plain", op::SymOperation{D}) where D
    opseitz, opxyzt = seitz(op), xyzt(op)
    print(io, opseitz, " ")
    printstyled(io, repeat('─',max(38-length(opseitz)-length(opxyzt), 1)),
                    " (", opxyzt, ")"; color=:light_black)
    
    # don't print matrix format if we IOContext is :compact=>true or dimension is 1
    if get(io, :compact, false) || D == 1
        return nothing
    end 
    
    println(io)
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
                cᴵ = convert(Int64, matrix(op)[i,j])
                printstyled(io, sep, cᴵ, color=:light_black)
            else
                # just use the same sep even if the symop is specified in a nonstandard basis (e.g.
                # cartesian); probably isn't a good general solution, but good enough for now
                printstyled(io, sep, round(c, digits=4), color=:light_black)
            end
        end
        printstyled(io, " ", i == 1 ? "╷" : (i == D ? "╵" : "┆"), " ", repeat(' ', Nsepτ-length(τstrs[i])), τstrs[i], " ", color=:light_black)
        printstyled(io, i == 1 ? '┐' : (i == D ? '┘' : '│'), color=:light_black) # close brace char
        i ≠ D && println(io)
    end
    return nothing
end
_has_negative_sign_and_isnonzero(x) = !iszero(x) && signbit(x)


# --- MultTable ---
function show(io::IO, ::MIME"text/plain", mt::MultTable)
    summary(io, mt)
    println(io, ":")
    mt.isgroup || (println(io, "  invalid multiplication table"); return nothing)
    seitz_ops = seitz.(mt.operations)
    pretty_table(io,
        [seitz_ops getindex.(Ref(seitz_ops), mt.table)], # 1st column and table itself
        vcat("", seitz_ops);                             # 1st row (header)
        highlighters = Highlighter((data,i,j) -> i==1 || j==1; bold=true),
        vlines = [1,], hlines = [:begin, 1, :end]
        )
    return nothing
end

# --- KVec ---
function show(io::IO, ::MIME"text/plain", v::AbstractVec)
    cnst, free = parts(v)
    print(io, '[')
    if isspecial(v)
        for i in eachindex(cnst) 
            coord = cnst[i] == -0.0 ? 0.0 : cnst[i] # normalize -0.0 to 0.0
            print(io, coord)
            # prepare for next coordinate/termination
            i == length(cnst) ? print(io, ']') : print(io, ", ")
        end
    else
        for i in eachindex(cnst)
            # constant/fixed parts
            if !iszero(cnst[i]) || iszero(@view free[i,:]) # don't print zero, if it adds unto anything nonzero
                coord = cnst[i] == -0.0 ? 0.0 : cnst[i] # normalize -0.0 to 0.0
                print(io, coord)
            end
            # free-parameter parts
            for j in eachindex(cnst) 
                if !iszero(free[i,j])
                    sgn = signaschar(free[i,j])
                    if !(iszero(cnst[i]) && sgn=='+' && iszero(free[i,1:j-1])) # don't print '+' if nothing precedes it
                        print(io, sgn)
                    end
                    if abs(free[i,j]) != oneunit(eltype(free)) # don't print prefactors of 1
                        print(io, abs(free[i,j]))
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
string(kv::AbstractVec) = (io=IOBuffer(); show(io, MIME"text/plain"(), kv); String(take!(io)))


# --- AbstractGroup ---
function summary(io::IO, g::T) where T<:AbstractGroup 
    print(io, T)
    if !(T <: GenericGroup)
        print(io, " #", num(g), " (", label(g), ")")
    end
    print(io, " with ", order(g), " operations")
end
function show(io::IO, ::MIME"text/plain", g::T) where T<:AbstractGroup
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    summary(io, g)
    println(io, ':')
    for (i,op) in enumerate(g)
        print(io, ' ')
        show(io, MIME"text/plain"(), op)
        if i < order(g); println(io); end
    end
end
function show(io::IO, ::MIME"text/plain", gs::AbstractVector{<:AbstractGroup})
    # TODO: This kind of show extension is bad style, afaik...
    Ngs = length(gs)
    for (i,g) in enumerate(gs)
        show(io, MIME"text/plain"(), g); 
        if i < Ngs; println(io); end
    end
end


# --- LGIrrep & PGIrrep ---
function show(io::IO, ::MIME"text/plain", plgir::Union{<:LGIrrep, <:PGIrrep})
    lgirlab = formatirreplabel(label(plgir))
    lablen = length(lgirlab)
    nindent = lablen+1
    prettyprint_header(io, lgirlab)
    prettyprint_irrep_matrices(io, plgir, nindent)
end
function show(io::IO, ::MIME"text/plain", plgirs::AbstractVector{T}) where T<:Union{<:LGIrrep, <:PGIrrep}
    # TODO: This kind of show extension is bad style, afaik...
    # Header line
    plg = group(first(plgirs))
    print(io, "$T: ")
    prettyprint_group_header(io, plg)

    Nᵢᵣ = length(plgirs)
    for (i,plgir) in enumerate(plgirs)
        show(io, MIME"text/plain"(), plgir)
        if i != Nᵢᵣ; println(io); end
    end
end

# ... utilities to print PGIrreps and LGIrreps
function prettyprint_group_header(io::IO, plg::AbstractGroup)
    print(io, "#", num(plg), " (", iuc(plg), ")")
    if plg isa LittleGroup
        print(io, " at " , klabel(plg), " = ")
        show(io, MIME"text/plain"(), kvec(plg))
    end
    println(io)
end
function prettyprint_scalar_or_matrix(io::IO, printP::AbstractMatrix, prefix::AbstractString,
                                      ϕabc_contrib::Bool=false)
    if size(printP) == (1,1) # scalar case
        v = printP[1]
        if isapprox(v, real(v), atol=DEFAULT_ATOL)          # real scalar
            if ϕabc_contrib && abs(real(v)) ≈ 1.0
                signbit(real(v)) && print(io, '-')
            else
                print(io, real(v))
            end
        elseif isapprox(v, imag(v)*im, atol=DEFAULT_ATOL)   # imaginary scalar
            if ϕabc_contrib && abs(imag(v)) ≈ 1.0
                signbit(imag(v)) && print(io, '-')
            else
                print(io, imag(v))
            end
            print(io, "i")
        else                                                # complex scalar (print as polar)
            vρ, vθ = abs(v), angle(v)
            vθ /= π
            print(io, vρ  ≈ 1.0 ? "" : vρ, "exp(") 
            if abs(vθ) ≈ 1.0
                signbit(vθ) && print(io, '-')
            else
                print(io, vθ)
            end
            print(io, "iπ)")
            #print(io, ϕabc_contrib ? "(" : "", v, ϕabc_contrib ? ")" : "")
        end

    else # matrix case
        formatter = x->(xr = real(x); xi = imag(x);
                        ComplexF64(abs(xr) > DEFAULT_ATOL ? xr : 0.0,
                                   abs(xi) > DEFAULT_ATOL ? xi : 0.0)) # round small complex components to zero
        # FIXME: not very optimal; e.g. makes a whole copy and doesn't handle displaysize
        compact_print_matrix(io, printP, prefix, formatter)
    end
end
function prettyprint_irrep_matrix(io::IO, lgir::LGIrrep, i::Integer, prefix::AbstractString)
    # unpack
    k₀, kabc = parts(lgir.lg.kv)
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
            abs(c) ≈ 0.5 || print(io, abs(2c)) # do not print if multiplicative factor is 1

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
                    abs(c) ≈ 0.5 || print(io, abs(2c)) # do not print if multiplicative factor is 1
                    print(io, 'ΰ'+i) # prints 'α', 'β', and 'γ' for i = 1, 2, and 3, respectively ('ΰ'='α'-1)
                end
            end
            print(io, ")]")
        end
    end
    
    # Least-effort way to indicate nontrivial (pseudo-real/complex) co-representations
    # TODO: Improve printing of pseudo-real and complex LGIrrep co-representations?
    if iscorep(lgir) 
        if type(lgir) == 2     # pseudo-real
            print(io, " + block-repetition")
        elseif type(lgir) == 3 # complex
            print(io, " + conjugate-block-repetition")
        else
            throw(DomainError(type, "Unexpected combination of iscorep=true and type≠{2,3}"))
        end
    end
end
function prettyprint_irrep_matrices(io::IO, plgir::Union{<:LGIrrep, <:PGIrrep}, 
                                  nindent::Integer, nboxdelims::Integer=45)  
    indent = repeat(" ", nindent)
    boxdelims = repeat("─", nboxdelims)
    linelen = nboxdelims + 4 + nindent
    Nₒₚ = order(plgir)
    for (i,op) in enumerate(operations(plgir))
        print(io, indent, " ├─ ")
        opseitz, opxyzt  = seitz(op), xyzt(op)
        printstyled(io, opseitz, ": ", 
                        repeat("─", linelen-11-nindent-length(opseitz)-length(opxyzt)),
                        " (", opxyzt, ")\n"; color=:light_black)
        #Base.print_matrix(IOContext(io, :compact=>true), ir, indent*(i == Nₒₚ ? " ╰" : " │")*"    ")
        print(io, indent, " │     ")
        prettyprint_irrep_matrix(io, plgir, i, indent*" │     ")
        if i < Nₒₚ; println(io, '\n', indent, " │     "); end
    end
    print(io, "\n", indent, " └", boxdelims)
end
function prettyprint_header(io::IO, plgirlab::AbstractString, nboxdelims::Integer=45)
    println(io, plgirlab, " ─┬", repeat("─", nboxdelims))
end


# --- CharacterTable ---
function show(io::IO, ::MIME"text/plain", ct::CharacterTable)
    chars = characters(ct)
    chars_formatted = Array{Union{Float64, ComplexF64, Int64, Complex{Int64}}}(undef, size(chars))
    for (idx, c) in enumerate(chars)
        chars_formatted[idx] = if isreal(c)
            isinteger(real(c)) ? convert(Int64, real(c)) : real(c)
        else
            ((isinteger(real(c)) && isinteger(imag(c))) 
                      ? convert(Int64, real(c)) + convert(Int64, imag(c))
                      : c)
        end
    end
    println(io, typeof(ct), ": ", tag(ct)) # type name and space group/k-point tags
    pretty_table(io,
        [seitz.(operations(ct)) chars_formatted], # 1st column: seitz operations; then formatted character table
        ["" formatirreplabel.(labels(ct))...];    # 1st row (header): irrep labels
        tf = tf_unicode,
        highlighters = Highlighter((data,i,j) -> i==1 || j==1; bold=true),
        vlines = [1,], hlines = [:begin, 1, :end]
        )
end


# --- BandRep ---
function prettyprint_symmetryvector(io::IO, irvec::Vector{<:Real}, irlabs::Vector{String})
    Nⁱʳʳ  = length(irlabs)
    Nⁱʳʳ′ = length(irvec) 
    if !(Nⁱʳʳ′ == Nⁱʳʳ || Nⁱʳʳ′ == Nⁱʳʳ+1)
        # we allow irvec to exceed the dimension of irlabs by 1, in case it includes dim(BR)
        throw(DimensionMismatch("irvec and irlabs must have matching dimensions"))
    end
    print(io, '[')

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
    print(io, ']')
end

summary(io::IO, BR::BandRep) = print(io, dim(BR), "-band BandRep (", label(BR), " at ", wyck(BR), ")")
function show(io::IO, ::MIME"text/plain", BR::BandRep)
    summary(io, BR)
    print(io, ":\n ")
    prettyprint_symmetryvector(io, vec(BR), irreplabels(BR))
end


# --- BandRepSet ---
function show(io::IO, ::MIME"text/plain", BRS::BandRepSet)
    Nⁱʳʳ = length(irreplabels(BRS))
    Nᵉᵇʳ = length(BRS)

    # print a "title" line and the irrep labels
    println(io, "BandRepSet (#", num(BRS), "): ",
                length(BRS), " BandReps, ", 
                "sampling ", Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(BRS) ? "½" : "1", " ",
                istimeinvar(BRS) ? "w/" : "w/o", " TR)")

    # print band representations as table
    k_idx = (i) -> findfirst(==(klabel(irreplabels(BRS)[i])), klabels(BRS)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_ν   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")
    pretty_table(io, 
        # table contents
        matrix(BRS, true),
        # header
        permutedims([wyck.(BRS) chop.(label.(BRS), tail=2)], (2,1)); # get rid of the repeating "↑G" part
        # row names
        row_names = vcat(irreplabels(BRS), "ν"),
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,], hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_name_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_ν), 
        header_crayon = crayon"bold"
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
    ν_maxdigs  = maximum(ndigits∘dim, BRS) # implicitly assuming this to be ≤2...
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
    print(io, ' '^indent, 'ν', ' '^maxcols_irr, '║')
    for idxᵇʳ in 1:room_for
        ν = dim(BRS[idxᵇʳ])
        addspace = div(cols_brlabs[idxᵇʳ], 2)+1
        print(io, ' '^(addspace-ndigits(ν)), ν)
        print(io, ' '^(cols_brlabs[idxᵇʳ] - addspace))
        print(io, idxᵇʳ ≠ room_for ? '│' : rowend)
    end
    =#

    #=
    # === EBRs-by-rows layout ===
    ν_maxdigs = maximum(ndigits∘dim, BRS)
    cols_brlab = maximum(x->length(label(x)), BRS)+1
    cols_irstart = cols_brlab+4
    cols_avail = displaysize(io)[2]-2                                 # available cols in io (cannot write to all of it; subtract 2)
    cols_requi = sum(x->length(x)+3, irreplabels(BRS))+cols_irstart+ν_maxdigs+3 # required cols for irrep labels & band reps
    if cols_requi > cols_avail
        cols_toomany    = ceil(Int64, (cols_requi-cols_avail)/2) + 2  # +2 is to make room for '  …  ' extender
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
    println(io, ' '^ν_maxdigs, "ν", " ║") # band-filling column header
    # print each bandrep
    for (bridx,BR) in enumerate(BRS)
        ν = dim(BR)
        print(io, "   ", label(BR),                      # bandrep label
                  " "^(cols_brlab-length(label(BR))), '║')
        for (iridx,x) in enumerate(vec(BR)) # vector representation of band rep
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
        
        print(io, ' '^(1+ν_maxdigs-ndigits(ν)), ν, " ║") # band-filling
        if bridx != length(BRS); println(io); end
    end
    =#
end