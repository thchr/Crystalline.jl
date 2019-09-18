# magic numbers
const BYTES_PER_KCHAR = 3
const BYTES_PER_KVEC  = 16*BYTES_PER_KCHAR;


parseisoir(T::Type{Real}) = parseisoir(Float64)         # just for being able to call it with Real or Complex
parseisoir(T::Type{Complex}) = parseisoir(ComplexF64)   # as input rather than Float64 and ComplexF64

function parseisoir(::Type{T}) where T<:Union{Float64,ComplexF64}
    datatag = if T <: Real; "PIR"; elseif T <: Complex; "CIR"; end   
    io = open((@__DIR__)*"/../data/ISOTROPY/"*datatag*"_data.txt","r") # open file for reading

    irreps = Vector{Vector{SGIrrep{T}}}()
    while !eof(io)
        # --- READ BASIC INFO, LIKE SG & IR #, NAMES, DIMENSIONALITY, ORDER ---
        irnum = parse(Int64, String(read(io, 5))) # read IR# (next 5 characters)
        sgnum = parse(Int64, String(read(io, 4))) # read SG# (next 4 characters)
        # read SGlabel, IRlabel, 
        skip(io, 2)
        sglabel = filter(!isequal(' '), readuntil(io, "\""))
        skip(io, 2)
        irlabel = roman2greek(filter(!isequal(' '), readuntil(io, "\"")))
        # read irdim, irtype, knum, pmknum, opnum (rest of line; split at spaces)
        #       irdim  : dimension of IR
        #       irtype : type of IR (see Sec. 5 of Acta Cryst. (2013). A69, 388)
        #                   - type 1 : intrinsically real IR
        #                   - type 2 : intrinsically complex IR, but equiv. to its 
        #                              own complex conjugate; pseudoreal
        #                   - type 3 : intrinsically complex IR, inequivalent to 
        #                              its own complex conjugate
        #       knum   : # of 𝐤-points in 𝐤-star
        #       pmknum : # of 𝐤-points in ±𝐤-star
        #       opnum  : number of symmetry elements in little group of 𝐤-star
        #                (i.e. order of little group of 𝐤-star)
        irdim, irtype, knum, pmknum, opnum = parsespaced(String(readline(io)))

        # --- READ VECTORS IN THE (±)𝐤-STAR (those not related by ±symmetry) ---
        # this is a weird subtlelty: Stokes et al store their 4×4 𝐤-matrices
        # in column major format; but they store their operators and irreps
        # in row-major format - we take care to follow their conventions on 
        # a case-by-case basis
        Nstoredk = T <: Real ? pmknum : knum # number of stored 𝐤-points depend on whether we load real or complex representations
        k = [Vector{Float64}(undef, 3) for n=Base.OneTo(Nstoredk)]
        kabc = [Matrix{Float64}(undef, 3,3) for n=Base.OneTo(Nstoredk)]
        for n = Base.OneTo(Nstoredk) # for Complex, loop over distinct vecs in 𝐤-star; for Real loop over distinct vecs in ±𝐤-star
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
            #       ∑³ᵢ₌₁ (kᵢ + aᵢα+bᵢβ+cᵢγ)*𝐆ᵢ     (w/ recip. basis vecs. 𝐆ᵢ)
            # here the matrix kabc is decomposed into vectors (𝐚,𝐛,𝐜) while α,β,γ are free
            # parameters ranging over all non-special values (i.e. not coinciding with high-sym 𝐤)
            k[n] = (@view kmat[1:3,1])./kdenom[1]                  # coefs of "fixed" parts 
            kabc[n] =  (@view kmat[1:3,2:4])./(@view kdenom[2:4])' # coefs of free parameters (α,β,γ)
        end
        kspecial = iszero(kabc[1]) # if no free parameters, i.e. 𝐚=𝐛=𝐜=𝟎 ⇒ high-symmetry 𝐤-point (i.e. "special") 
        checked_read2eol(io) # read to end of line & check we didn't miss anything

        # --- READ OPERATORS AND IRREPS IN LITTLE GROUP OF ±𝐤-STAR ---
        opmatrix = [Matrix{Float64}(undef, 3,4) for _=1:opnum]
        opxyzt   = Vector{String}(undef, opnum)
        irtranslation = [zeros(Float64, 3) for _=1:opnum]
        irmatrix = [Matrix{T}(undef, irdim,irdim) for _=1:opnum]
        for i = 1:opnum
            # --- OPERATOR ---
            # matrix form of the symmetry operation (originally in a 4×4 form; the [4,4] idx is a common denominator)
            optempvec   = parsespaced(readline(io))
            opmatrix[i] = rowmajorreshape(optempvec, (4,4))[1:3,:]./optempvec[16] # surprisingly, this is in row-major form..!
            # note the useful convention that the nonsymmorphic translation always ∈[0,1[; in parts of Bilbao, components are 
            # occasionally negative; this makes construction of multtables unnecessarily cumbersome
            opxyzt[i]   = matrix2xyzt(opmatrix[i]) 

            # --- ASSOCIATED IRREP ---
            if !kspecial # if this is a general position, we have to incorporate a translational modulation in the point-part of the irreps
                transtemp = parsespaced(readline(io))
                irtranslation[i] = transtemp[1:3]./transtemp[4]
            else
                irtranslation[i] = zeros(Float64, 3)
                # TODO: Use this to create the appropriate "translation-modulation" matrix for 
                #       nonspecial kvecs (see rules in https://stokes.byu.edu/iso/irtableshelp.php)
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
        irrep = SGIrrep{T}(irnum,    irlabel,    irdim,
                           sgnum,    sglabel,
                           irtype,   opnum,      
                           knum,     pmknum,   kspecial,
                           KVec.(k, kabc),  
                           SymOperation.(opxyzt, opmatrix),
                           irtranslation, irmatrix)
        if length(irreps) < sgnum; push!(irreps, Vector{SGIrrep{T}}()); end # new space group idx
        push!(irreps[sgnum], irrep)
      
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

const tabfloats = (sqrt(3)/2, sqrt(2)/2, sqrt(3)/4, cos(π/12)/sqrt(2), sin(π/12)/sqrt(2))
""" 
    reprecision_data(x::Float64) --> Float64

Stokes et al. used a table to convert integers to floats; in addition, 
the floats were truncated on writing. We can restore their precision 
by checking if any of the relevant entries occur, and then returning 
their untruncated floating point value. See also `tabfloats::Tuple`.

The possible floats that can occur in the irrep tables are:

        ┌ 0,1,-1,0.5,-0.5,0.25,-0.25 (parsed with full precision)
        │ ±0.866025403784439 => ±sqrt(3)/2
        │ ±0.707106781186548 => ±sqrt(2)/2
        │ ±0.433012701892219 => ±sqrt(3)/4
        │ ±0.683012701892219 => ±cos(π/12)/sqrt(2)
        └ ±0.183012701892219 => ±sin(π/12)/sqrt(2)
"""
function reprecision_data(x::T) where T<:Real
    absx = abs(x)
    for preciseabsx in tabfloats
        if isapprox(absx, preciseabsx, atol=1e-4) 
            return copysign(preciseabsx, x)
        end
    end
    return x
end
reprecision_data(z::T) where T<:Complex = complex(reprecision_data(real(z)), reprecision_data(imag(z)))

function littlegroupirrep(ir::SGIrrep{<:Complex})
    lgidx, lgops = littlegroup(operations(ir), kstar(ir)[1], centering(ir.sgnum,3))
    lgirdim′ = ir.dim/ir.knum; lgirdim = div(ir.dim, ir.knum)
    @assert lgirdim′ == lgirdim "The dimension of the little group irrep must be an integer, equaling "*
                                "the dimension of the space group irrep divided by the number of vectors "*
                                "in star{𝐤}"
    # broadcasting to get all the [1:lgirdim, 1:lgirdim] blocks of every irrep assoc. w/ the lgidx list
    lgir = getindex.((@view irreps(ir)[lgidx]), Ref(Base.OneTo(lgirdim)), Ref(Base.OneTo(lgirdim))) 
    lgtrans = ir.translations[lgidx]

    return LGIrrep(num(ir), label(ir), kstar(ir)[1], lgops, lgir, lgtrans, type(ir))
end

parselittlegroupirreps() = parselittlegroupirreps.(parseisoir(Complex))
function parselittlegroupirreps(irvec::Vector{SGIrrep{ComplexF64}})
    lgirvec = Vector{Tuple{LGIrrep,Vararg{LGIrrep}}}()
    curlab = nothing; accidx = Int64[]
    for (idx, ir) in enumerate(irvec) # loop over distinct irreps (e.g., Γ1, Γ2, Γ3, Z1, Z2, ..., GP1)
        if curlab == klabel(ir)
            push!(accidx, idx)
        else
            if curlab != nothing
                lgirs = Vector{LGIrrep}(undef, length(accidx))
                for (pos, kidx) in enumerate(accidx) # write all irreps of a specific k-point to a vector (e.g., Z1, Z2, ...)
                    lgirs[pos] = littlegroupirrep(irvec[kidx])
                end
                push!(lgirvec, (lgirs...,))
            end

            curlab = klabel(ir)
            accidx = [idx,]
        end
    end
    # after the loop finishes, one batch of k-point irreps still needs 
    # incorporation (because we're always _writing_ a new batch, when 
    # we've moved into the next one); for ISOTROPY's default sorting, 
    # this is the GP=Ω=[α,β,γ]ᵀ point)
    lgirs = Vector{LGIrrep}(undef, length(accidx))
    for (pos, kidx) in enumerate(accidx)
        lgirs[pos] = littlegroupirrep(irvec[kidx])
    end
    push!(lgirvec, (lgirs...,))

    return lgirvec
end

"""
    write_littlegroupirreps(lgirsvec::Vector{Tuple{LGIrrep}})
                                                    --> Nothing

Write all little group small irreps associated with a specific space 
group to disk, as JSON files, to ease subsequent loading of little group 
small irreps. Takes a vector of little group small irreps of the sort
    `lgirsvec::Vector{Tuple{LGIrrep}}`
i.e., vector-indexed across distinct k-points and tuple-indexed across
distinct irreps; in practice, calling 
    `write_littlegroupirreps.(parselittlegroupirreps())`
will write **all** the little group irreps to disk.
"""
function write_littlegroupirreps(lgirsvec)
    sgnum = num(first(first(lgirsvec)))
    Nk = length(lgirsvec)

    # build up lists of KVec and SymOperation info
    klab_list = Vector{String}(undef, Nk)
    kv_list   = Vector{KVec}(undef, Nk)
    ops_list  = Vector{T where T<:Vector{String}}(undef, Nk) 
    for (kidx, lgirs) in enumerate(lgirsvec) # lgirs is a tuple of LGIrreps, all at the same 𝐤-point
        lgir = first(lgirs) # 𝐤-info is the same for each LGIrrep in tuple lgirs
        klab_list[kidx] = klabel(lgir)
        kv_list[kidx] = kvec(lgir)
        ops_list[kidx] = xyzt.(operations(lgir))
    end

    filename_kvecs = (@__DIR__)*"/../data/lgirreps/3d/kinfo"*string(sgnum)*".jld"

    bson(filename_kvecs, klab_list = klab_list, 
                         kv_list = kv_list,
                         ops_list = ops_list)

    # write irreps
    matrices_list = [Vector{Matrix{ComplexF64}}() for _=1:Nk]
    translations_list = [Vector{Vector{Float64}}() for _=1:Nk]
    type_list = [Vector{Int}() for _=1:Nk]
    matrices_list = [[lgir.matrices for lgir in lgirs] for lgirs in lgirsvec]
    translations_list = [[lgir.translations for lgir in lgirs] for lgirs in lgirsvec]
    type_list = [[lgir.type for lgir in lgirs] for lgirs in lgirsvec]

    filename_irreps = (@__DIR__)*"/../data/lgirreps/3d/irreps"*string(sgnum)*".jld"
    bson(filename_irreps, matrices_list = matrices_list, 
                          translations_list = translations_list,
                          type_list = type_list)

    return nothing
end
write_littlegroupirreps() = write_littlegroupirreps.(parselittlegroupirreps())



const TEST_αβγ = [0.123,0.456,0.789] # arbitrary test numbers for KVecs
# TODO: This implementation should follow the discussion on p. 650-652 in Bradley 
#       & Cracknell's book (there's some discussion in 622-626 as well, but that's 
#       for point groups). Their discussion is for magnetic groups but is generally 
#       applicable, and is by far the most clear and thorough discussion that I've 
#       found so far.
#       Cornwell also does a good job of explicating this.
#       Inui on p. 296-299 also discuss it, but is less clear overall.
function realify(irs::NTuple{Nirr, LGIrrep}, verbose::Bool=false) where Nirr
    kv = kvec(first(irs)) # must be the same for all irreps in list
    kv_αβγ = kv(TEST_αβγ)
    sgnum = num(first(irs))
    lgops = operations(first(irs))
    Nops = order(first(irs)) # order of little group (= # of operations)

    d = dim(kv)
    cntr = centering(sgnum, d)
    sgops = operations(get_symops(sgnum, d))
    star = starofk(sgops, kv, cntr)

    verbose && print(klabel(first(irs)), " │ ")

    # Check if -𝐤 is in the star of 𝐤, or if 𝐤 is equivalent to -𝐤: 
    # if so, TR is an element of the little group; if not, it isn't 
    # ║ 𝐑𝐞𝐚𝐬𝐨𝐧: if there is an element g of the (unitary) 𝑠𝑝𝑎𝑐𝑒 group G   
    # ║   that takes 𝐤 to -𝐤 mod 𝐆, then (denoting the TR element by Θ, 
    # ║   acting as θ𝐤 = -𝐤) the antiunitary element θg will take 𝐤 to  
    # ║   𝐤 mod 𝐆, i.e. θg will be an element of the little group of 𝐤
    # ║   M(k) associated with the 𝑔𝑟𝑎𝑦 space group M ≡ G + θG.
    # ║   Conversely, if no such element g exists, there can be no anti-
    # ║   unitary elements in the little group derived from M; as a result, 
    # ║   TR is not part of the little group and so does not modify its 
    # ║   small irreps (called "co-reps" for magnetic groups).
    # ║   There can then only be type 'x' degeneracy (between 𝐤 and -𝐤)
    # ║   but TR will not change the degeneracy at 𝐤 itself.
    if !isapproxin(-kv, star, cntr; atol=DEFAULT_ATOL)
        corep_idxs = [[i] for i in Base.OneTo(Nirr)] # TR ∉ M(k) ⇒ smalls irrep (... small co-reps) not modified by TR
        verbose && println(klabel(first(irs)), "ᵢ ∀i (type x) ⇒  no additional degeneracy (star{k} ∌ -k)")

    else
        # Test if 𝐤 is equivalent to -𝐤, i.e. if 𝐤 = -𝐤 + 𝐆
        k_equiv_kv₋ = isapprox(-kv, kv, cntr; atol=DEFAULT_ATOL)

        # Find an element in G that takes 𝐤 → -𝐤 (if 𝐤 is equivalent to -𝐤, 
        # then this is just the unit-element I (if `sgops` is sorted conven-
        # tionally, with I first, this is indeed what the `findfirst(...)`  
        # bits below will find)
        if !k_equiv_kv₋
            g₋ = sgops[findfirst(g-> isapprox(g∘kv, -kv, cntr; atol=DEFAULT_ATOL), sgops)]
        else
            # This is a bit silly: if k_equiv_kv₋ = true, we will never use g₋; but I'm not sure if 
            # the compiler will figure that out, or if it will needlessly guard against missing g₋?
            g₋ = SymOperation(hcat(I, zeros(d))) # ... the unit element I
        end

        # -𝐤 is part of star{𝐤}; we infer reality of irrep from ISOTROPY's data (could also 
        # be done using `herring(...)`). ⇒ deduce new small irreps (... small co-reps).
        corep_idxs = Vector{Vector{Int64}}()
        skiplist = Vector{Int64}()
        for (i, ir) in enumerate(irs)
            if i ∈ skiplist; continue; end # already matched to this irrep previously; i.e. already included now
            verbose && i ≠ 1 && print("  │ ")

            if type(ir) == 1     # real
                push!(corep_idxs, [i])
                if verbose
                    println(formatirreplabel(label(ir)), " (real) ⇒  no additional degeneracy")
                end

            elseif type(ir) == 2 # pseudo-real
                # doubles irrep on its own
                push!(corep_idxs, [i, i])
                if verbose
                    println(formatirreplabel(label(ir)^2), " (pseudo-real) ⇒  doubles degeneracy"); 
                end

            elseif type(ir) == 3 # complex
                # In this case, there must exist a "partner" irrep (say, Dⱼ) which is 
                # equal to the complex conjugate of the current irrep (say, Dᵢ); we 
                # next search for this equivalence.
                # When we check for equivalence between irreps Dᵢ* and Dⱼ we must
                # account for the possibility of a 𝐤-dependence in the matrix-form
                # of the irreps; specifically, for an element g, its small irrep is
                #     Dᵢ[g] = exp(2πik⋅τᵢ[g])Pᵢ[g],
                # where, crucially, for symmetry lines, planes, and general points
                # 𝐤 depends on (one, two, and three) free parameters (α,β,γ).
                # Thus, for equivalence of irreps Dᵢ* and Dⱼ we require that
                #     Dᵢ*[g] ~ Dⱼ[g]       ∀g ∈ G(k)
                #  ⇔ exp(-2πik⋅τᵢ[g])Pᵢ*[g] ~ exp(2πik⋅τⱼ[g])Pⱼ[g]
                # It seems rather tedious to prove that this is the case for all 𝐤s
                # along a line/plane (α,β,γ). Rather than attempt this, we simply test
                # against an arbitrary value of (α,β,γ) [superfluous entires are ignored]
                # that is non-special (i.e. not ={0,0.5,1}); this is `TEST_αβγ`.

                # Characters of the conjugate of Dᵢ, i.e. tr(Dᵢ*) = tr(Dᵢ)*
                θχᵢ = conj.(tr.(irreps(ir, TEST_αβγ))) 
                
                # Find matching complex partner
                partner = 0
                for j = i+1:Nirr
                    if j ∉ skiplist && type(irs[j]) == 3 # only check if j has not previously matched; 
                                                         # similarly, only check if the jth irrep is complex.

                        # Note that we require only equivalence of Dᵢ* and Dⱼ; not equality. 
                        # Cornwell describes (p. 152-153 & 188) a neat trick for checking this 
                        # efficiently: specifically, Dᵢ* and Dⱼ are equivalent irreps if 
                        #     χⁱ(g)* = χʲ(g₋⁻¹gg₋) ∀g ∈ G(k)
                        # with g₋ an element of G that takes 𝐤 to -𝐤, and where χⁱ (χʲ) denotes
                        # the characters the respective irreps.
                        χⱼ = tr.(irreps(irs[j], TEST_αβγ))
                        match = true
                        for n in Base.OneTo(Nops)
                            if k_equiv_kv₋ # 𝐤 = -𝐤 + 𝐆 ⇒ g₋ = I (the unit element), s.t. g₋⁻¹gg₋ = I⁻¹gI = g
                                χⱼ_g₋⁻¹gg₋ = χⱼ[n]
                            else           # 𝐤 not equivalent to -𝐤, i.e. 𝐤 ≠ -𝐤 + 𝐆
                                g₋⁻¹gg₋ = compose(compose(inv(g₋), lgops[n], false), g₋, false)
                                n′, Δw = findequiv(g₋⁻¹gg₋, lgops, cntr)
                                χⱼ_g₋⁻¹gg₋ = exp(2π*1im*dot(kv_αβγ, Δw)) .* χⱼ[n′]
                            end
                            
                            match = isapprox(θχᵢ[n], χⱼ_g₋⁻¹gg₋; atol=DEFAULT_ATOL)
                            if !match # ⇒ not a match
                                break
                            end
                        end

                        if match # ⇒ a match
                            partner = j
                            if verbose; 
                                println(formatirreplabel(label(ir)*label(irs[j])), " (complex) ⇒  doubles degeneracy")
                            end
                        end
                    end
                end
                partner === 0 && throw(ErrorException("Didn't find a matching complex partner for $(label(ir))"))
                push!(skiplist, partner)

                push!(corep_idxs, [i, partner])
                
            else
                throw(ArgumentError("Invalid real/pseudo-real/complex type = $(type(ir))"))
            end
        end
    end

    Ncoreps = length(corep_idxs)

    # New small co-rep labels (composite)
    newlabs = Tuple(join(label(irs[i]) for i in corep_idxs[i′]) for i′ in Base.OneTo(Ncoreps))

    # TODO: New small irreps (small co-reps)
    #=
    for i′ in Base.OneTo(Ncoreps)
        idxs = coreps_idxs[i′]
        if length(idxs) == 1      # real or type x
            # same as before
        elseif idxs[1] == idxs[2] # pseudoreal 
            # doubles self
        else                      # complex
            # doubles with complex conjugate
            # what to do about exp(ikτ) dependence? Need new type, different from LGIrrep?
        end
    end
    =#
    return corep_idxs
end


"""
    herring(ir::LGIrrep, sgops::AbstractVector{SymOperation},
            αβγ::Union{Vector{<:Real},Nothing}=nothing)        --> Tuple{Int, Int}

Computes the Herring criterion for a little group irrep `ir`, from 

        ∑ χ({β|b}²) 
over symmetry operations {β,b} that take k → -k.

The provided space group operations `sgops` **must** be the set reduced by 
primitive translation vectors; i.e. using `get_symops(...)` directly is **not** 
allowable in general. Using the operations from the Γ point of ISOTROPY's 
dataset is, however, fine.

As a sanity check, a value of `αβγ` can be provided to check for invariance
along a symmetry line/plane/general point in k-space. Obviously, the reality 
type should invariant to this choice.

**Implementation:** 
See e.g. Inui's Eq. (13.48), Dresselhaus, p. 618, and 
and Herring's original paper at https://doi.org/10.1103/PhysRev.52.361.
We mainly followed Cornwell, p. 150-152 & 187-188.
"""
function herring(ir::LGIrrep, sgops::AbstractVector{SymOperation}, αβγ::Union{Vector{<:Real},Nothing}=nothing)

    lgops = operations(ir)
    kv = kvec(ir)
    kv₋ = -kv
    dim = length(kv.k₀)
    cntr = centering(num(ir), dim)
    Ds = irreps(ir, αβγ) # irrep matrices
    kv_αβγ = kv(αβγ)

    s = zero(ComplexF64)
    for op in sgops
        if isapprox(op∘kv, kv₋, cntr, atol=DEFAULT_ATOL) # check if op∘k == -k; if so, include in sum
            op² = compose(op, op, false) # this is op∘op, _including_ trivial lattice translation parts
            # find the equivalent of `op²` in `lgops`; this may differ by a number of 
            # primitive lattice vectors `w_op²`; the difference must be included when 
            # we calculate the trace of the irrep 𝐃: the irrep matrix 𝐃 is ∝exp(2πi𝐤⋅𝐭)
            idx_of_op²_in_lgops, Δw_op² = findequiv(op², lgops, cntr)
            ϕ_op² = exp(2π*1im*dot(kv_αβγ, Δw_op²)) # phase accumulated by "trivial" lattice translation parts
            χ_op² = ϕ_op²*tr(Ds[idx_of_op²_in_lgops]) # χ(op²)

            s += χ_op²
        end
    end

    pgops = pointgroup(sgops) # point group assoc. w/ space group
    g₀ = length(pgops) # order of pgops (denoted h, or macroscopic order, in Bradley & Cracknell)
    Mk = length(starofk(pgops, kv, cntr)) # order of star of k (denoted qₖ in Bradley & Cracknell)
    normalization = round(Int, g₀/Mk) # order of G₀ᵏ; the point group derived from the little group Gᵏ (denoted b in Bradley & Cracknell; [𝐤] in Inui)
    if !isapprox(normalization, g₀/Mk)
        throw(ErrorException("The little group is not factored by its point group and star{k}: this should never happen"))
    end

    # check that output is a real integer and then convert to that for output...
    if norm(imag(s)) < DEFAULT_ATOL 
        sInt = round(Int,real(s)); 
    else 
        throw(error("Herring criterion should yield a real value; obtained complex s=$(s)")) 
    end
    if norm(sInt-real(s)) > DEFAULT_ATOL 
        throw(error("Herring criterion should yield an integer; obtained s=$(s)"))
    end
    return sInt, normalization # this is ∑ χ({β|b}²) and g₀/M(k) in Cornwell's Eq. (7.18)
end

"""
    findequiv(op::SymOperation, ops::AbstractVector{SymOperation}, cntr::Char) 
                                                --> Tuple{Int, Vector{Float64}}

Search for an operator `op′` in `ops` which is equivalent, modulo differences
by **primitive** lattice translations `Δw`, to `op`. Return the index of `op′` in 
`ops`, as well as the primitive translation difference `Δw`.

The small irreps of `op` at wavevector k, Dⱼᵏ[`op`], can be computed from 
the small irreps of `op′`, Dⱼᵏ[`op′`], via Dⱼᵏ[`op`] = exp(2πik⋅`Δw`)Dⱼᵏ[`op′`]
"""
function findequiv(op::SymOperation, ops::AbstractVector{SymOperation}, cntr::Char)
    W = rotation(op)
    w = translation(op)

    P = primitivebasismatrix(cntr, dim(op))
    w′ = P\w    # `w` in its primitive basis

    for (j, opⱼ) in enumerate(ops)
        Wⱼ = rotation(opⱼ)
        wⱼ = translation(opⱼ)
        wⱼ′ = P\w

        if W == Wⱼ # rotation-part of op and opⱼ is identical
            # check if translation-part of op and opⱼ is equivalent, modulo a primitive lattice translation
            if all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), w′.-wⱼ′)
                return j, w.-wⱼ
            end
        end
    end
    return nothing, nothing # didn't find any match
end


"""
    isapprox(kv1::KVec, kv2::KVec, cntr::Char; kwargs...) 
                                                            --> Bool
                                            
Compute approximate equality of two KVec's `k1` and `k2` modulo any 
primitive G-vectors. To ensure that primitive G-vectors are used, 
the centering type `cntr` (see `centering(cntr, dim)`) must be given
(the dimensionality is inferred from `kv1` and `kv2`).
Optionally, keyword arguments (e.g., `atol` and `rtol`) can be 
provided, to include in calls to `Base.isapprox`.
"""
function isapprox(kv1::KVec, kv2::KVec, cntr::Char; kwargs...)
    k₀1, kabc1 = parts(kv1) # ... unpacking
    k₀2, kabc2 = parts(kv2)

    dim1, dim2 = length(k₀1), length(k₀2)
    if dim1 ≠ dim2
        throw(ArgumentError("dim(kv1)=$(dim1) and dim(kv2)=$(dim2) must be equal"))
    end

    # check if k₀ ≈ k₀′ differ by a _primitive_ 𝐆 vector
    diff = primitivebasismatrix(cntr, dim1)' * (k₀1 .- k₀2)
    kbool = all(el -> isapprox(el, round(el); kwargs...), diff) 
    # check if kabc1 ≈ kabc2; no need to check for difference by a 
    # 𝐆 vector, since kabc is in interior of BZ
    abcbool = isapprox(kabc1, kabc2;  kwargs...)

    return kbool && abcbool
end