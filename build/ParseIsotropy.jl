module ParseIsotropy

using Crystalline
using Crystalline: AbstractIrrep
using Pkg.Artifacts: ensure_artifact_installed

export parselittlegroupirreps, isotropy_datadir

"""
    isotropy_datadir() --> String

The directory holding ISOTROPY's `CIR_data.txt` and `PIR_data.txt`, downloading them first
if necessary: they are not kept in the repository, but published as a release asset and
declared as the `isotropy` artifact (see `Artifacts.toml` and `build/DATA-RELEASE.md`).

`@artifact_str` is not usable here — it searches upwards from this file for an
`Artifacts.toml` and stops at `build/Project.toml` — so the declaration is named explicitly.
"""
function isotropy_datadir()
    toml = joinpath(dirname(@__DIR__), "Artifacts.toml")
    return ensure_artifact_installed("isotropy", toml)
end

# --- Magic numbers ---
const BYTES_PER_KCHAR = 3
const BYTES_PER_KVEC  = 16*BYTES_PER_KCHAR;

# --- 3D Space group irrep struct ---
struct SGIrrep3D{T} <: AbstractIrrep{3} where T
    iridx::Int      # sequential index assigned to ir by Stokes et al
    cdml::String    # CDML label of irrep (including 𝐤-point label)
    irdim::Int      # dimensionality of irrep (i.e. size)
    sgnum::Int      # space group number
    sglabel::String # Hermann-Mauguin label of space group
    reality::Reality # real, pseudo-real, or complex
    order::Int      # number of operations
    knum::Int       # number of 𝐤-vecs in star
    pmknum::Int     # number of ±𝐤-vecs in star
    special::Bool   # whether star{𝐤} describes high-symmetry points
    pmkstar::Vector{KVec{3}}     # star{𝐤} for Complex, star{±𝐤} for Real
    ops::Vector{SymOperation{3}} # every symmetry operation in space group
    translations::Vector{Vector{Float64}} # translations assoc with matrix repres of ops in irrep
    matrices::Vector{Matrix{T}}  # non-translation assoc with matrix repres of ops in irrep
end
num(sgir::SGIrrep3D) = sgir.sgnum
(sgir::SGIrrep3D)(αβγ=nothing) = sgir.matrices # TODO: add missing αβγ dependence
order(sgir::SGIrrep3D) = sgir.order
iuc(sgir::SGIrrep3D) = sgir.sglabel
operations(sgir::SGIrrep3D) = sgir.ops
isspecial(sgir::SGIrrep3D) = sgir.special
orbit(sgir::SGIrrep3D) = sgir.pmkstar
irdim(sgir::SGIrrep3D) = sgir.irdim
dim(sgir::SGIrrep3D) = 3


# --- Parsing ---

parseisoir(T::Type{Real}; kws...) = parseisoir(Float64; kws...)         # just for being able to call it with Real or Complex
parseisoir(T::Type{Complex}; kws...) = parseisoir(ComplexF64; kws...)   # as input rather than Float64 and ComplexF64

function parseisoir(::Type{T}; datadir = isotropy_datadir()) where T<:Union{Float64,ComplexF64}
    datatag = if T <: Real; "PIR"; elseif T <: Complex; "CIR"; end   
    io = open(joinpath(datadir, datatag*"_data.txt"), "r")

    irreps = Vector{Vector{SGIrrep3D{T}}}()
    while !eof(io)
        # --- READ BASIC INFO, LIKE SG & IR #, NAMES, DIMENSIONALITY, ORDER ---
        irnum = parse(Int, String(read(io, 5))) # read IR# (next 5 characters)
        sgnum = parse(Int, String(read(io, 4))) # read SG# (next 4 characters)
        # read SGlabel, IRlabel, 
        skip(io, 2)
        sglabel = filter(!isequal(' '), readuntil(io, "\""))
        skip(io, 2)
        irlabel = Crystalline.formatirreplabel(
                    Crystalline.roman2greek(filter(!isequal(' '), readuntil(io, "\""))))
        # read irdim, irreality, knum, pmknum, opnum (rest of line; split at spaces)
        #       irdim  : dimension of IR
        #       irreality : reality of IR (see Sec. 5 of Acta Cryst. (2013). A69, 388)
        #                   - 1 : intrinsically real IR
        #                   - 2 : intrinsically complex IR, but equiv. to its 
        #                         own complex conjugate; pseudoreal
        #                   - 3 : intrinsically complex IR, inequivalent to 
        #                         its own complex conjugate
        #       knum   : # of 𝐤-points in 𝐤-star
        #       pmknum : # of 𝐤-points in ±𝐤-star
        #       opnum  : number of symmetry elements in little group of 𝐤-star
        #                (i.e. order of little group of 𝐤-star)
        irdim, irreality′, knum, pmknum, opnum = parsespaced(String(readline(io)))
        
        # --- REMAP REALITY TYPE TO USUAL CONVENTION ---
        # ISOTROPY's convention is: 1 => REAL, 2 => PSEUDOREAL, and 3 => COMPLEX.
        # Our convention (consistent with the Herring and Frobenius-Schur criteria) is
        # instead that REAL = 1, PSEUDOREAL = -1, and COMPLEX = 0, as in encoded in the 
        # Reality enum. We make the swap here, so we get a single consistent convention.
        irreality = fix_isotropy_reality_convention(irreality′)

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
        irtranslation = [zeros(Float64, 3) for _=1:opnum]
        irmatrix = [Matrix{T}(undef, irdim,irdim) for _=1:opnum]
        for i = 1:opnum
            # --- OPERATOR ---
            # matrix form of the symmetry operation (originally in a 4×4 form; the [4,4] idx is a common denominator)
            optempvec   = parsespaced(readline(io))
            opmatrix[i] = rowmajorreshape(optempvec, (4,4))[1:3,:]./optempvec[16] # surprisingly, this is in row-major form..!
            # note the useful convention that the nonsymmorphic translation always ∈[0,1[; in parts of Bilbao, components are 
            # occasionally negative; this makes construction of `MultTable`s unnecessarily cumbersome

            # --- ASSOCIATED IRREP ---
            if !kspecial # if this is a general position, we have to incorporate a translational modulation in the point-part of the irreps
                transtemp = parsespaced(readline(io))
                irtranslation[i] = transtemp[1:3]./transtemp[4]
            else
                irtranslation[i] = zeros(Float64, 3)
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
        end

        # --- STORE DATA IN VECTOR OF IRREPS ---
        irrep = SGIrrep3D{T}(irnum, irlabel, irdim, sgnum, sglabel, irreality,    
                             opnum, knum, pmknum, kspecial, KVec.(k, kabc),
                             SymOperation{3}.(opmatrix), irtranslation, irmatrix)
        length(irreps) < sgnum && push!(irreps, Vector{SGIrrep3D{T}}()) # new space group idx
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
E.g. `parsespaced(Int, "  1  2  5") = [1, 2, 5]`
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
@inline parsespaced(s::AbstractString) = parsespaced(Int, s)



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
    s = readline(io)
    # check that the string indeed is empty; i.e. we didn't miss anything
    @assert isempty(s) "Parsing error; unexpected additional characters after expected EOL"

    return nothing
end

function rowmajorreshape(v::AbstractVector, dims::Tuple)
    return PermutedDimsArray(reshape(v, dims), reverse(ntuple(i->i, length(dims))))
end

function fix_isotropy_reality_convention(iso_reality::Integer)
    iso_reality == 1 && return REAL
    iso_reality == 2 && return PSEUDOREAL
    iso_reality == 3 && return COMPLEX
    throw(DomainError(iso_reality, "unreachable: unexpected ISOTROPY reality value"))
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

function littlegroupirrep(ir::SGIrrep3D{<:Complex})
    lgidx, lgops = littlegroup(operations(ir), orbit(ir)[1], centering(num(ir),3))
    lgirdim′ = irdim(ir)/ir.knum; lgirdim = div(irdim(ir), ir.knum)
    @assert lgirdim′ == lgirdim "The dimension of the little group irrep must be an integer, equaling "*
                                "the dimension of the space group irrep divided by the number of vectors "*
                                "in star{𝐤}"

    kv = orbit(ir)[1] # representative element of the k-star; the k-vector of assoc. w/ this little group   
    if !is_erroneous_lgir(num(ir), label(ir), 3)
        # broadcasting to get all the [1:lgirdim, 1:lgirdim] blocks of every irrep assoc. w/ the lgidx list
        lgirmatrices = getindex.((@view ir()[lgidx]), Ref(Base.OneTo(lgirdim)), Ref(Base.OneTo(lgirdim))) 
        lgirtrans = ir.translations[lgidx]
    else
        #println("Manually swapped out corrected (CDML) LGIrrep for sgnum ", num(ir), ", irrep ", label(ir))
        lgirmatrices, lgirtrans = manually_fixed_lgir(num(ir), label(ir))
    end

    return LGIrrep{3}(label(ir), LittleGroup(num(ir), kv, klabel(ir), collect(lgops)), lgirmatrices, lgirtrans, reality(ir), false)
end

parselittlegroupirreps(; kws...) = parselittlegroupirreps.(parseisoir(Complex; kws...))
function parselittlegroupirreps(irvec::Vector{SGIrrep3D{ComplexF64}})
    lgirsd = Dict{String, Collection{LGIrrep{3}}}()
    curklab = nothing; accidx = Int[]
    for (idx, ir) in enumerate(irvec) # loop over distinct irreps (e.g., Γ1, Γ2, Γ3, Z1, Z2, ..., GP1)
        if curklab == klabel(ir)
            push!(accidx, idx)
        else
            if curklab !== nothing
                lgirs = Vector{LGIrrep{3}}(undef, length(accidx))
                for (pos, kidx) in enumerate(accidx) # write all irreps of a specific k-point to a vector (e.g., Z1, Z2, ...)
                    lgirs[pos] = littlegroupirrep(irvec[kidx])
                end
                push!(lgirsd, curklab=>Collection(lgirs))
            end

            curklab = klabel(ir)
            accidx = [idx,]
        end
    end
    # after the loop finishes, one batch of k-point irreps still needs 
    # incorporation (because we're always _writing_ a new batch, when 
    # we've moved into the next one); for ISOTROPY's default sorting, 
    # this is the GP=Ω=[α,β,γ]ᵀ point)
    lgirs = Vector{LGIrrep{3}}(undef, length(accidx))
    for (pos, kidx) in enumerate(accidx)
        lgirs[pos] = littlegroupirrep(irvec[kidx])
    end
    lastklab = klabel(irvec[last(accidx)])
    @assert lastklab == "Ω"
    push!(lgirsd, lastklab=>Collection(lgirs))

    return lgirsd
end


const ERRONEOUS_LGIRS = (214=>"P₁", 214=>"P₂", 214=>"P₃") # extend to tuple of three-tuples if we ever need D ≠ 3 as well
@inline function is_erroneous_lgir(sgnum::Integer, irlab::String, D::Integer=3)
    D == 1 && return false
    D ≠ 3 && Crystalline._throw_2d_not_yet_implemented(D)
    @simd for ps in ERRONEOUS_LGIRS
        (ps[1]==sgnum && ps[2]==irlab) && return true
    end 
    return false
end

"""
    manually_fixed_lgir(sgnum::Integer, irlab::String)

The small irreps associated with the little group of k-point P ≡ KVec(½,½,½)
of space group 214 are not correct in ISOTROPY's dataset: specifically, while 
they have the correct characters and pass the 1st and 2nd character orthogonality
theorems, they do not pass the grand orthogonality theorem that tests the irrep
matrices themselves. To that end, we manually replace these small irreps with 
those listed by CDML (read off from their tables). 
Those irreps are manually extracted in the scripts/cdml_sg214_P1P2P3.jl file.

The fix is made in littlegroupirrep(ir::SGIrrep3D{<:Complex}), using the check 
in is_erroneous_lgir(...), with the constant "erroneous" tuple ERRONEOUS_LGIRS.

Emailed Stokes & Campton regarding the issue on Sept. 26, 2019; did not yet 
hear back.
"""
function manually_fixed_lgir(sgnum::Integer, irlab::String)
    # TODO: Use their new and corrected dataset (from February 17, 2020) instead of manually
    #       fixing the old dataset.
    #       I already verified that their new dataset is correct (and parses), and that P₁
    #       and P₃ pass tests. Their P₁ and P₃ (and P₂) irreps are not the same as those we
    #       had below, but differ by a transformation R(THEIRS)R⁻¹ = (OURS) with 
    #       R = [1 1; 1 -1]/√2.
    #       Thus, we should remove this method (here and from other callers), update and
    #       commit the new datasets and then finally regenerate/refresh our own saved format
    #       of the irreps (from build/write_littlegroup_irreps_from_ISOTROPY.jl)
    if sgnum == 214
        CP  = cis(π/12)/√2   # C*P       ≈ 0.683013 + 0.183013im
        CQ  = cis(5π/12)/√2  # C*Q       ≈ 0.183013 + 0.683013im
        CcP = cis(-π/12)/√2  # C*conj(P) ≈ 0.683013 - 0.183013im
        CcQ = cis(-5π/12)/√2 # C*conj(Q) ≈ 0.183013 - 0.683013im
        if irlab == "P₁"
            matrices = [[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im],     # x,y,z
                        [0.0+0.0im 1.0+0.0im; 1.0+0.0im 0.0+0.0im],     # x,-y,-z+1/2
                        [0.0+0.0im 0.0-1.0im; 0.0+1.0im 0.0+0.0im],     # -x+1/2,y,-z
                        [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im],    # -x,-y+1/2,z
                        [CP CcQ; CP -CcQ],                              # z,x,y
                        [CcP CcP; CQ -CQ],                              # y,z,x
                        [CcP -CcP; CQ CQ],                              # -y+1/2,z,-x
                        [CP CcQ; -CP CcQ],                              # -z,-x+1/2,y
                        [CcP CcP; -CQ CQ],                              # -y,-z+1/2,x
                        [CP -CcQ; CP CcQ],                              # z,-x,-y+1/2
                        [CQ -CQ; CcP CcP],                              # y,-z,-x+1/2
                        [CcQ CP; -CcQ CP]]                              # -z+1/2,x,-y
        elseif irlab == "P₂" 
            # there is, as far as I can see, nothing wrong with ISOTROPY's (214, P₂)
            # small irrep: but it doesn't agree with the form we extract from CDML 
            # either. To be safe, and to have a consistent set of irreps for the P 
            # point we just swap out this irrep as well.
            matrices =  [[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im],    # x,y,z
                         [0.0+0.0im 0.0+1.0im; 0.0-1.0im 0.0+0.0im],    # x,-y,-z+1/2
                         [0.0+0.0im -1.0+0.0im; -1.0+0.0im 0.0+0.0im],  # -x+1/2,y,-z
                         [-1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im],   # -x,-y+1/2,z
                         [-0.5-0.5im -0.5-0.5im; 0.5-0.5im -0.5+0.5im], # z,x,y
                         [-0.5+0.5im 0.5+0.5im; -0.5+0.5im -0.5-0.5im], # y,z,x
                         [0.5-0.5im 0.5+0.5im; 0.5-0.5im -0.5-0.5im],   # -y+1/2,z,-x
                         [0.5+0.5im 0.5+0.5im; 0.5-0.5im -0.5+0.5im],   # -z,-x+1/2,y
                         [0.5-0.5im -0.5-0.5im; -0.5+0.5im -0.5-0.5im], # -y,-z+1/2,x
                         [0.5+0.5im -0.5-0.5im; -0.5+0.5im -0.5+0.5im], # z,-x,-y+1/2
                         [-0.5-0.5im 0.5-0.5im; 0.5+0.5im 0.5-0.5im],   # y,-z,-x+1/2
                         [-0.5+0.5im 0.5-0.5im; 0.5+0.5im 0.5+0.5im]]   # -z+1/2,x,-y
        elseif irlab == "P₃"
            matrices = [[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im],     # x,y,z
                        [0.0+0.0im 1.0+0.0im; 1.0+0.0im 0.0+0.0im],     # x,-y,-z+1/2
                        [0.0+0.0im 0.0-1.0im; 0.0+1.0im 0.0+0.0im],     # -x+1/2,y,-z
                        [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im],    # -x,-y+1/2,z
                        [-CQ -CcP; -CQ CcP],                            # z,x,y
                        [-CcQ -CcQ; -CP CP],                            # y,z,x
                        [-CcQ CcQ; -CP -CP],                            # -y+1/2,z,-x
                        [-CQ -CcP; CQ -CcP],                            # -z,-x+1/2,y
                        [-CcQ -CcQ; CP -CP],                            # -y,-z+1/2,x
                        [-CQ CcP; -CQ -CcP],                            # z,-x,-y+1/2
                        [-CP CP; -CcQ -CcQ],                            # y,-z,-x+1/2
                        [-CcP -CQ; CcP -CQ]]                            # -z+1/2,x,-y
        else
            throw(DomainError((sgnum, irlab), "should not be called with these input; nothing to fix"))
        end
        translations = [zeros(Float64, 3) for _=Base.OneTo(length(matrices))]

        return matrices, translations
    else
        throw(DomainError((sgnum, irlab), "should not be called with these input; nothing to fix"))
    end
end



# TODO: We need to manually do things for the "special" k-points kᴮ from
# the representation domain Φ that do not exist in the basic domain kᴬ∈Ω
# i.e. for kᴮ∈Φ-Ω (labels like KA, ZA, etc.), which are not included in
# ISOTROPY. We can do this by first finding a transformation R such that
# kᴮ = Rkᴬ, and then transforming the LGIrrep of kᴬ appropriately (see
# CDML Sec. 4.1 or B&C Sec. 5.5.). Parts of the steps are like this:
#
#       kvmaps = ΦnotΩ_kvecs(sgnum, D) # from src/special_representation_domain_kpoints.jl
#       if kvmaps !== nothing # contains KVecs in the representation domain Φ that 
#                           # cannot be mapped to ones in the basic domain Ω 
#           for kvmap in kvmaps # loop over each "new" KVec
#               
#               cdmlᴮ = kvmap.kᴮlab # CDML label of "new" KVec kᴮ∈Φ-Ω
#               cdmlᴬ = kvmap.kᴮlab # CDML label of "old" KVec kᴬ∈Ω
#               R     = kvmap.op    # Mapping from kᴬ to kᴮ: kᴮ = Rkᴬ
#               # pick kᴬ irreps in lgirsd
#               lgirsᴬ = lgirsd[cdmlᴬ]
#           end
#           # ... do stuff to lgirsᴬ to get lgirsᴮ via a transformation {R|v}
#           # derived from a holosymmetric parent group of sgnum and the transformation R
#       end
#
# The only kᴮ included in ISOTROPY is Z′≡ZA for sgs 195, 198, 200, 201, & 205: 
# all other kᴮ∈Φ-Ω points are omitted)
#
# Pa3 (Tₕ⁶), sg 205, cannot be treated by the transformation method and requires  
# manual treatment (B&C p. 415-417); fortunately, the Z′₁=ZA₁ irrep of 205 is  
# already included in ISOTROPY.
#=
function add_special_representation_domain_lgirs(lgirsd::Dict{String, <:AbstractVector{LGIrrep{D}}}) where D
    D ≠ 3 && Crystalline._throw_1d2d_not_yet_implemented(D)

    sgnum = num(first(first(lgirsd)))

    # does this space group contain any nontrivial k-vectors in Φ-Ω?
    
    # what is the method for adding in these missing k-vectors?
            
end
=#


end # module