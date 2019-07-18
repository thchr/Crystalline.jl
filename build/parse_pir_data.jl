# magic numbers
const BYTES_PER_KCHAR = 3
const BYTES_PER_KVEC  = 16*BYTES_PER_KCHAR;

function parsepir(nirreps=10294; verbose=false)
    i = 1 # how many irreps we parsed; should agree with irnum
    io = open((@__DIR__)*"/../data/ISOTROPY_PIR/PIR_data.txt","r")
    irlabels = Vector{String}(undef, nirreps); highsympoints = Vector{Char}()
    while i ≤ nirreps && !eof(io)
        # --- READ BASIC INFO, LIKE SG & IR #, NAMES, DIMENSIONALITY, ORDER ---
        # read IR# (next 5 characters)
        irnum = parse(Int64, String(read(io, 5)))
        # read SG# (next 4 characters)
        sgnum = parse(Int64, String(read(io, 4)))

        # read SGlabel, IRlabel, 
        skip(io, 2)
        sglabel = filter(!isspace, readuntil(io, "\""))
        skip(io, 2)
        irlabel = filter(!isspace, readuntil(io, "\""))
        irlabels[i] = irlabel
        # read irdim, irtype, knum, pmknum, opnum (rest of line; split at spaces)
        #       irdim  : dimension of IR
        #       irtype : type of IR (see Sec. 5 of Acta Cryst. (2013). A69, 388)
        #                   - type 1 : intrinsically real IR
        #                   - type 2 : intrinsically complex IR, but equiv. to its 
        #                              own complex conjugate; pseudoreal
        #                   - type 3 : intrinsically complex IR, inequivalent to 
        #                              its own complex conjugate
        #       knum   : # of k-points in k-star
        #       pmknum : # of k-points in ±k-star
        #       opnum  : number of symmetry elements in little group of k-star 
        #                (i.e. order of little group of k-star)
        irdim, irtype, knum, pmknum, opnum = parsespaced(String(readline(io)))

        # --- READ VECTORS IN THE ±K-STAR (those not related by ±symmetry) ---
        k = [Vector{Float64}(undef, 3) for n=1:pmknum]
        kαβγ = [Matrix{Float64}(undef, 3,3) for n=1:pmknum]
        kspecial = true
        if kspecial == true
            highsympoints = append!(highsympoints, irlabel[1])
            unique!(highsympoints)
        end
        for n = 1:pmknum # loop over distinct kvecs in ±k-star
            if n == 1 # we don't have to worry about \n then, and can use faster read()
                kmat = reshape(parsespaced(String(read(io, BYTES_PER_KVEC))), (4,4))
            else
                kmat = reshape(parsespaced(readexcept(io, BYTES_PER_KVEC)), (4,4))
            end

            # a general k-vector is specified as a pair (k, kαβγ), denoting a k-vector
            #       \sum_i=1^3 (k[i] + a[i]α+b[i]β+g[i]γ)*G[i]     (w/ recip. basis vecs. G[i])
            # here the matrix kαβγ is decomposed into vectors (a,b,g) while α,β,γ are 
            # free parameters ranging over (-½,½) [TODO: I think!, but cannot confirm].
            k[n] = (@view kmat[1:3,1])./kmat[4,1]
            kαβγ[n] =  (@view kmat[1:3,2:4])./(@view kmat[4,2:4])' # free parameter coefficients
            if n == 1
                kspecial = iszero(kαβγ[n]) # if no free parameters, this is a high-symmetry point
            end
        end
        checked_read2eol(io) # read to end of line & check we didn't miss anything

        # --- READ OPERATORS AND IRREPS IN LITTLE GROUP OF ±K-STAR ---
        opmatrix = [Array{Float64}(undef, 3,4) for i=1:opnum]
        irtranslation = []
        irmatrix = [Array{Float64}(undef, irdim,irdim) for i=1:opnum]
        for i = 1:opnum
            # --- OPERATOR ---
            # matrix form of the symmetry operation (originally in a 4x4 form; the [4,4] idx is the common denominator)
            optempvec = parsespaced(readline(io))
            opmatrix[i] = reshape(optempvec, (4,4))[1:3,:]./optempvec[16]
            # TODO opshorthand[i] = 

            # --- ASSOCIATED IRREP ---
            if !kspecial # if this is a general position, we have to incorporate a translational contribution to the irreps
                irtranslation = append!(irtranslation, parsespaced(readline(io))) 
                # TODO: Interpret this according to the rules in https://stokes.byu.edu/iso/irtableshelp.php
            else
                irtranslation = append!(irtranslation, 1)
            end
            # irrep matrix "base"
            tempvec = Vector{Float64}()
            while length(tempvec) != irdim^2
                tempvec = append!(tempvec, parsespaced(Float64, readline(io)))
            end
            irmatrix[i] = reshape(tempvec, (irdim,irdim))
            # TODO: Might want to "up" the precision of truncated floating point numbers according to this table
            #     0,1,-1,0.5,-0.5,0.25,-0.25, (don't require any adjustments)
            #     ±0.866025403784439 => ±sqrt(3)/2
            #     ±0.707106781186548 => ±sqrt(2)/2
            #     ±0.433012701892219 => ±sqrt(3)/4
            #     ±0.683012701892219 => ±cos(15)/sqrt(2)
            #     ±0.183012701892219 => ±sin(15)/sqrt(2)

            # TODO: ir_character[i] = trace(irmatrix[i]*irtranslation[i])
        end
    
        if verbose
            @show irnum, sgnum
            @show sglabel, irlabel
            @show irdim, irtype, knum, pmknum, opnum  
            @show k, kspecial
            @show opmatrix
            @show irtranslation
            @show irmatrix
            println()
            println(irnum)
        end
        i += 1
    end
    close(io)
    return unique!(irlabels)
end

""" 
    parsespaced(T::Type, s::AbstractString)

    Parses a string `s` with spaces interpreted as delimiters, split-
    ting at every contiguious block of spaces and returning a vector
    of the split elements, with elements parsed as type `T`.
    E.g. `parsespaced(Int64, "  1  2  5") = [1, 2, 5]`
"""
parsespaced(T::Type, s::AbstractString) = parse.(T, split(s, r"\s+",keepempty=false))
parsespaced(s::AbstractString) = parsespaced(Int64, s)

""" readexcept(s::IO, nb::Integer, except::Char; all=true)

Same as `read(s::IO, nb::Integer; all=true)` but allows us ignore byte matches to 
those in `except`.
"""
function readexcept(io::IO,  nb::Integer, except::Union{Char,Nothing}='\n'; all=true)
    out = IOBuffer(); n = 0
    while n < nb && !eof(io)
        c = read(io, Char)
        if c == except
            continue
        end
        write(out, c)
        n += 1
    end
    return String(take!(out))
end


function checked_read2eol(io) # read to end of line & check that this is indeed the last character
    s = readuntil(io, "\n")
    if !isempty(s); # move to next line & check we didn't miss anything
        error(s,"Parsing error; unexpected additional characters after expected end of line"); 
    end
end


@time parsepir(10294)