using SGOps, LinearAlgebra

# magic numbers
const BYTES_PER_KCHAR = 3
const BYTES_PER_KVEC  = 16*BYTES_PER_KCHAR;


parseisoir(T::Type{Real}) = parseisoir(Float64)         # just for being able to call it without thinking about Float64 or not...
parseisoir(T::Type{Complex}) = parseisoir(ComplexF64)

function parseisoir(::Type{T}) where T<:Union{Float64,ComplexF64} #, verbose=false
    datatag = if T <: Real; "PIR"; elseif T <: Complex; "CIR"; end   
    io = open((@__DIR__)*"/../data/ISOTROPY/"*datatag*"_data.txt","r") # open file for reading

    irreps = Vector{Vector{Irrep{T}}}()
    while !eof(io)
        # --- READ BASIC INFO, LIKE SG & IR #, NAMES, DIMENSIONALITY, ORDER ---
        irnum = parse(Int64, String(read(io, 5))) # read IR# (next 5 characters)
        sgnum = parse(Int64, String(read(io, 4))) # read SG# (next 4 characters)
        # read SGlabel, IRlabel, 
        skip(io, 2)
        sglabel = filter(!isequal(' '), readuntil(io, "\""))
        skip(io, 2)
        irlabel = filter(!isequal(' '), readuntil(io, "\""))
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

        # --- READ VECTORS IN THE (±)𝐤-STAR (those not related by ±symmetry) ---
        # this is a weird subtlelty: Stokes et al store their 4×4 𝐤-matrices
        # in column major format; but they store their operators and irreps
        # in row-major format - we take care to follow their conventions on 
        # a case-by-case basis
        Nstoredk = T <: Real ? pmknum : knum # number of stored k-points depend on whether we load real or complex representations
        k = [Vector{Float64}(undef, 3) for n=Base.OneTo(Nstoredk)]
        kabc = [Matrix{Float64}(undef, 3,3) for n=Base.OneTo(Nstoredk)]
        for n = Base.OneTo(Nstoredk) # for Complex, loop over distinct vecs in 𝐤-star; for Real loop over distinct vecs in ±k-star
            if n == 1 # we don't have to worry about '\n' then, and can use faster read()
                kmat = reshape(parsespaced(String(read(io, BYTES_PER_KVEC))), (4,4))     # load as column major matrix
            else
                kmat = reshape(parsespaced(readexcept(io, BYTES_PER_KVEC, '\n')), (4,4)) # load as column major matrix
            end
            # Stokes' conventions implicitly assign NaN columns to unused free parameters when any 
            # other free parameters are in play; here, we prefer to just have zero-columns. To do 
            # that, we change common denominators to 1 if they were 0 in Stokes' data.
            kdenom = [(iszero(origdenom) ? 1 : origdenom) for origdenom in kmat[4,:]] 

            # 𝐤-vectors are specified as a pair (k, kabc), denoting a 𝐤-vector
            #       \sum_i=1^3 (k[i] + a[i]α+b[i]β+c[i]γ)*𝐆[i]     (w/ recip. basis vecs. 𝐆[i])
            # here the matrix kabc is decomposed into vectors (a,b,c) while α,β,γ are free
            # parameters ranging over all non-special values (i.e. not coinciding with high-sym 𝐤)
            k[n] = (@view kmat[1:3,1])./kdenom[1]                  # coefs of "fixed" parts 
            kabc[n] =  (@view kmat[1:3,2:4])./(@view kdenom[2:4])' # coefs of free parameters (α,β,γ)
        end
        kspecial = iszero(kabc[1]) # if no free parameters, i.e. a=b=c=0 ==> high-symmetry 𝐤-point (i.e. "special") 
        checked_read2eol(io) # read to end of line & check we didn't miss anything

        # --- READ OPERATORS AND IRREPS IN LITTLE GROUP OF ±𝐤-STAR ---
        opmatrix = [Matrix{Float64}(undef, 3,4) for i=1:opnum]
        opxyzt   = Vector{String}(undef, opnum)
        irtranslation = Vector{Union{Nothing, Vector{Float64}}}(nothing, opnum)
        irmatrix = [Matrix{T}(undef, irdim,irdim) for i=1:opnum]
        for i = 1:opnum
            # --- OPERATOR ---
            # matrix form of the symmetry operation (originally in a 4x4 form; the [4,4] idx is a common denominator)
            optempvec   = parsespaced(readline(io))
            opmatrix[i] = rowmajorreshape(optempvec, (4,4))[1:3,:]./optempvec[16] # surprisingly, this is in row-major form..!
            opxyzt[i]   = matrix2xyzt(opmatrix[i]) 

            # --- ASSOCIATED IRREP ---
            if !kspecial # if this is a general position, we have to incorporate a translational modulation in the point-part of the irreps
                transtemp = parsespaced(readline(io))
                irtranslation[i] = transtemp[1:3]./transtemp[4]
                # TODO: Use this to create the appropriate "translation-modulation" matrix for nonspecial kvecs (see rules in https://stokes.byu.edu/iso/irtableshelp.php)
            end
            
            # irrep matrix "base" (read next irdim^2 elements into matrix)
            elcount1 = elcount2 = 0
            while elcount1 != irdim^2
                tempvec = parsespaced(T, readline(io))
                elcount2 += length(tempvec)
                irmatrix[i][elcount1+1:elcount2] = tempvec
                elcount1 = elcount2
            end
            irmatrix[i] = permutedims(irmatrix[i], (2,1)) # we loaded as column-major, but Stokes et al use row-major (unlike conventional Fortran)
            irmatrix[i] .= reprecision_data.(irmatrix[i])
            
            # TODO: ir_character[i] = tr(irmatrix[i]*irtranslation[i])

        end

        # --- WRITE DATA TO VECTOR OF IRREPS ---
        irrep = Irrep{T}(irnum,    irlabel,    irdim,
                         sgnum,    sglabel,
                         irtype,   opnum,      
                         knum,     pmknum,   kspecial,
                         collect(zip(k, kabc)),  
                         SymOperation.(opxyzt, opmatrix),
                         irtranslation, irmatrix)
        if length(irreps) < sgnum; push!(irreps, Vector{Irrep{T}}()); end # new space group idx
        push!(irreps[sgnum], irrep)
      

        # --- TESTING OF IRREP ---
        #=
        irchar = tr.(irmatrix) 
        sumchar2 = sum(abs2.(irchar))
        #symmorph = all([iszero(op[:,end]) for op in opmatrix])
        
        if irlabel[1:2] == "GM"
            if sumchar2 != opnum
                println(sumchar2 == opnum, ": ", sum(abs.(irchar).^2), " | ", opnum, " | ", "(", pmknum," ", knum,")")
            end
        end
        #if sgnum == 5 && irlabel[1:2] == "LD"
        #    for (i,irm) in enumerate(irmatrix)
        #        println(irlabel)
        #        println(opxyzt[i])
        #        display(irm)
        #        display(irtranslation)
        #        println()
        #    end
        #end
        =#
        
        #--- DEBUG PRINTING ---
        #=
        if verbose
            @show irnum, sgnum
            @show sglabel, irlabel
            @show irdim, irtype, knum, pmknum, opnum  
            @show k, kspecial
            @show opmatrix
            @show irtranslation
            @show irmatrix
            println()
        end
        =#

        # --- FINISHED READING CURRENT IRREP; MOVE TO NEXT ---
    end
    close(io)
    return irreps
end

""" 
    parsespaced(T::Type, s::AbstractString)

    Parses a string `s` with spaces interpreted as delimiters, split-
    ting at every contiguious block of spaces and returning a vector
    of the split elements, with elements parsed as type `T`.
    E.g. `parsespaced(Int64, "  1  2  5") = [1, 2, 5]`
"""
#@inline parsespaced(T::Type, s::AbstractString) = parse.(T, split(s, r"\s+", keepempty=false))
@inline function parsespaced(T::Type{<:Number}, s::AbstractString) 
    spacesplit=split(s, r"\s+", keepempty=false)
    if T <: Complex
        for (i,el) in enumerate(spacesplit)
            if el[1]=='('
                @inbounds spacesplit[i] = replace(replace(el[2:end-1],",-"=>"-"),','=>'+')*"i"
            end
        end
    end 
    return parse.(T, spacesplit)
end
@inline parsespaced(s::AbstractString) = parsespaced(Int64, s)



""" 
    readexcept(s::IO, nb::Integer, except::Char; all=true)

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

function rowmajorreshape(v::AbstractVector, dims::Tuple)
    return PermutedDimsArray(reshape(v, dims), reverse(ntuple(i->i, length(dims))))
end

const datafloats = (sqrt(3)/2, sqrt(2)/2, sqrt(3)/4, cos(15)/sqrt(2), sin(15)/sqrt(2))
""" 
    reprecision_data(x::Float64) --> Float64

    Stokes et al. used a table to convert integers to floats; in addition, 
    the floats were truncated on writing. We can restore their precision 
    by checking if any of the relevant entries occur, and then returning 
    their untruncated floating point value. See also `datafloats::Tuple`.
      | The possible floats that can occur in the irrep tables are:
      |     0,1,-1,0.5,-0.5,0.25,-0.25, (parsed with full precision)
      |     ±0.866025403784439 => ±sqrt(3)/2
      |     ±0.707106781186548 => ±sqrt(2)/2
      |     ±0.433012701892219 => ±sqrt(3)/4
      |     ±0.683012701892219 => ±cos(15)/sqrt(2)
      |     ±0.183012701892219 => ±sin(15)/sqrt(2)
"""
function reprecision_data(x::T) where T<:Real
    absx = abs(x)
    for preciseabsx in datafloats
        if isapprox(absx, preciseabsx, atol=1e-4) 
            return copysign(preciseabsx, x)
        end
    end
    return x
end
reprecision_data(z::T) where T<:Complex = complex(reprecision_data(real(z)), reprecision_data(imag(z)))



function littlegroupirrep(ir::Irrep{<:Complex})

    if !isspecial(ir)
        warning("didn't deal with non-special k-points for little groups yet")
        return nothing
    end

    k₀ = kstar(ir)[1][1]
    littleidx, littleops = littlegroup(operations(ir), k₀)
    knum = ir.knum; 
    lirdim = ir.dim/knum
    if !(lirdim ≈ round(lirdim)); 
        error("..."); 
    else
        lirdim = round(Int64,lirdim)
    end

    lir = getindex.((@view irreps(ir)[littleidx]), Ref(Base.OneTo(lirdim)), Ref(Base.OneTo(lirdim)))

    return littleops, lir
end




IR = parseisoir(Complex);


#CIR = parseisoir(Complex)
#PIR = parseisoir(Real)

## 
#[translation.(operations(ir)) for ir in IR[136]]
#parsepir(200)
#@btime parsepir(10294)
#parsepir(100)