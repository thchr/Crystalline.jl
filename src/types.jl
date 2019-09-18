import Base: show

# Crystalline lattice
struct Crystal{N}
    R::NTuple{N,Vector{Float64}}
end
Crystal(R,type) = Crystal{length(R)}(R, type)
Crystal(R) = Crystal{length(R)}(R, "")
basis(C::Crystal) = C.R
dim(C::Crystal{N}) where N = N
function show(io::IO, ::MIME"text/plain", C::Crystal)
    print(io, "$(dim(C))D Crystal:")
    print(io, " ($(crystalsystem(C)))");
    for (i,R) in enumerate(basis(C))
        print(io, "\n   R$(i): "); print(io, R); 
    end
end
norms(C::Crystal) = norm.(basis(C))
_angle(rA,rB) = acos(dot(rA,rB)/(norm(rA)*norm(rB)))
function angles(C::Crystal{N}) where N
    R = basis(C)
    γ = _angle(R[1], R[2])
    if N == 3
        α = _angle(R[2], R[3])
        β = _angle(R[3], R[1])
        return α,β,γ
    end
    return γ
end


# Symmetry operations
struct SymOperation
    xyzt::String
    matrix::Matrix{Float64} # TODO: factor this out into a rotation and a translation part (to avoid needless copying)
end
SymOperation(s::AbstractString) = SymOperation(string(s), xyzt2matrix(s))
SymOperation(m::Matrix{<:Real}) = SymOperation(matrix2xyzt(m), float(m))
matrix(op::SymOperation) = op.matrix
xyzt(op::SymOperation) = op.xyzt
dim(op::SymOperation) = size(matrix(op),1)
function show(io::IO, ::MIME"text/plain", op::SymOperation) 
    print(io, seitz(op),":\n   (", xyzt(op), ")\n")
    Base.print_matrix(IOContext(io, :compact=>true), op.matrix, "   ")
end
getindex(op::SymOperation, keys...) = matrix(op)[keys...]   # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
lastindex(op::SymOperation, d::Int64) = size(matrix(op), d) # allows using `end` in indices
rotation(m::Matrix{Float64}) = m[:,1:end-1] # rotational (proper or improper) part of an operation
rotation(op::SymOperation) = matrix(op)[:,1:end-1]        
translation(m::Matrix{Float64}) = m[:,end]  # translation part of an operation
translation(op::SymOperation) = matrix(op)[:,end]   
issymmorph(op::SymOperation) = iszero(translation(op))
(==)(op1::SymOperation, op2::SymOperation) = (xyzt(op1) == xyzt(op2)) && (matrix(op1) == matrix(op2))
isapprox(op1::SymOperation, op2::SymOperation; kwargs...) = isapprox(matrix(op1), matrix(op2); kwargs...)

# Multiplication table
struct MultTable
    operations::Vector{SymOperation}
    indices::Matrix{Int64}
    isgroup::Bool
end
indices(mt::MultTable) = mt.indices
isgroup(mt::MultTable) = mt.isgroup
function show(io::IO, ::MIME"text/plain", mt::MultTable)
    Base.print_matrix(IOContext(io, :compact=>true), mt.indices, "  ")
    print(io, "\nFor operations:\n  ")
    for (i,op) in enumerate(mt.operations)
        print(io, i, " => ", xyzt(op), "\t") # separation could be improved...
        if mod(i,4) == 0; print(io,"\n  "); end
    end
end
getindex(mt::MultTable, keys...) = indices(mt)[keys...]
lastindex(mt::MultTable, d::Int64) = size(indices(mt),d)

# Space group
struct SpaceGroup
    num::Int64
    operations::Vector{SymOperation}
    dim::Int64
end
num(sg::SpaceGroup) = sg.num
operations(sg::SpaceGroup) = sg.operations
dim(sg::SpaceGroup) = sg.dim
order(sg::SpaceGroup) = length(operations(sg))
function show(io::IO, ::MIME"text/plain", sg::SpaceGroup)
    Nops = order(sg)
    groupprefix = dim(sg) == 3 ? "Space" : (dim(sg) == 2 ? "Plane" : nothing)
    println(io, groupprefix, " group #", num(sg))
    for (i,op) in enumerate(operations(sg))
        show(io, "text/plain", op)
        if i < Nops; print(io, "\n\n"); end
    end
end
function show(io::IO, ::MIME"text/plain", sgs::Vector{SpaceGroup})
    Nsgs = length(sgs)
    for (i,sg) in enumerate(sgs); 
        show(io, "text/plain", sg); 
        if i < Nsgs; print(io, '\n'); end
    end
end

# Fourier/plane wave lattice (specified by Gorbits and coefficient interrelations)
struct FourierLattice{N}
    orbits::Vector{Vector{SVector{N, Int64}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}    # Vector of interrelations between coefficients of 𝐆-plane waves within an orbit
end
FourierLattice(orbits, orbitcoefs)   = FourierLattice{length(first(first(orbits)))}(orbits, orbitcoefs)
dim(flat::FourierLattice{N}) where N = N

# K-vectors
# 𝐤-vectors are specified as a pair (k₀, kabc), denoting a 𝐤-vector
#       𝐤 = ∑³ᵢ₌₁ (k₀ᵢ + aᵢα+bᵢβ+cᵢγ)*𝐆ᵢ     (w/ recip. basis vecs. 𝐆ᵢ)
# here the matrix kabc is columns of the vectors (𝐚,𝐛,𝐜) while α,β,γ are free
# parameters ranging over all non-special values (i.e. not coinciding with any 
# high-sym 𝐤)
struct KVec
    k₀::Vector{Float64}
    kabc::Matrix{Float64}
end
KVec(k₀::AbstractVector{<:Real}) = KVec(float.(k₀), zeros(Float64, length(k₀), length(k₀)))
parts(kv::KVec) = (kv.k₀, kv.kabc)
isspecial(kv::KVec) = iszero(kv.kabc)
dim(kv::KVec) = length(kv.k₀)
(kv::KVec)(αβγ::Vector{<:Real}) = begin
    k₀, kabc = parts(kv)
    return k₀ + kabc*αβγ
end
(kv::KVec)(αβγ::Vararg{<:Real, 2}) = kv([αβγ[1], αβγ[2]])
(kv::KVec)(αβγ::Vararg{<:Real, 3}) = kv([αβγ[1], αβγ[2], αβγ[3]])
(kv::KVec)() = kv.k₀
(kv::KVec)(::Nothing) = kv.k₀

function string(kv::KVec)
    k₀, kabc = parts(kv)
    buf = IOBuffer()
    write(buf, "[")
    if isspecial(kv)
        for i in eachindex(k₀) 
            coord = k₀[i] == -0.0 ? 0.0 : k₀[i] # normalize -0.0 to 0.0
            @printf(buf, "%g", coord)
            # prepare for next coordinate/termination
            i == length(k₀) ? write(buf, "]") : write(buf, ", ")
        end
    else
        for i in eachindex(k₀)
            # fixed parts
            if !iszero(k₀[i]) || iszero(@view kabc[i,:]) # don't print zero, if it adds unto anything nonzero
                coord = k₀[i] == -0.0 ? 0.0 : k₀[i] # normalize -0.0 to 0.0
                @printf(buf, "%g", coord)
            end
            # free-parameter parts
            for j in eachindex(k₀) 
                if !iszero(kabc[i,j])
                    sgn = signaschar(kabc[i,j])
                    if !(iszero(k₀[i]) && sgn=='+' && iszero(kabc[i,1:j-1])) # don't print '+' if nothing precedes it
                        write(buf, sgn)
                    end
                    if abs(kabc[i,j]) != oneunit(eltype(kabc)) # don't print prefactors of 1
                        @printf(buf, "%g", abs(kabc[i,j]))
                    end
                    write(buf, j==1 ? 'α' : (j == 2 ? 'β' : 'γ'))
                end
            end
            # prepare for next coordinate/termination
            i == length(k₀) ? write(buf, "]") : write(buf, ", ")
        end
    end
    return String(take!(buf))
end
show(io::IO, ::MIME"text/plain", kv::KVec) = print(io, string(kv))

""" 
    KVec(str::AbstractString) --> KVec

Construct a `KVec` struct from a string representations of a *k*-vector, supplied 
in either of the formats
        `"(\$x,\$y,\$z)"`, `"[\$x,\$y,\$z]"`, `"\$x,\$y,\$z"`,
where the coordinates `x`,`y`, and `z` are strings that can contain fractions,
decimal numbers, and "free" parameters {`'α'`,`'β'`,`'γ'`} (or, alternatively,
{`'u'`,`'v'`,`'w'`}). Returns the associated `KVec`.

Any "fixed"/constant part of a coordinate _must_ precede any free parts, e.g.,
`x="1+α"` is allowable but `x="α+1"` is not.
"""
function KVec(str::AbstractString)
    str = replace(strip(str, ['(',')','[',']']), ' '=>"") # tidy up string (remove parens & spaces)
    xyz = split(str,',')
    dim = length(xyz)
    k₀ = zeros(Float64, dim); kabc = zeros(Float64, dim, dim)
    for (i, coord) in enumerate(xyz)
        # "free" coordinates, kabc[i,:]
        for (j, matchgroup) in enumerate([['α','u'],['β','v'],['γ','w']])
            pos₂ = findfirst(x->any(y->y==x, matchgroup), coord)
            if !isnothing(pos₂)
                match = searchpriornumerals(coord, pos₂)
                kabc[i,j] = parse(Float64, match)
            end
        end

        # "fixed" coordinate, k₀[i]
        if !any(x->x==last(first(split(coord, r"\b(\+|\-)"))), ['α','u','β','v','γ','w']) # check for situations like '±3α' which is not handled by logic below
            nextidx = 0
            while (nextidx=nextind(coord, nextidx)) ≤ lastindex(coord) && !any(x->coord[nextidx]==x, ['α','u','β','v','γ','w'])
                if nextidx != 1 && any(x->coord[nextidx]==x, ['+','-'])
                    break
                else 
                end
            end
            if nextidx != firstindex(coord)
                k₀[i] = parsefraction(coord[firstindex(coord):prevind(coord,nextidx)])
            end
        end
    end
    return KVec(k₀, kabc)
end

# arithmetic with k-vectors
(-)(kv::KVec) = KVec(.- kv.k₀, .- kv.kabc)
(-)(kv1::KVec, kv2::KVec) = KVec(kv1.k₀ .- kv2.k₀, kv1.kabc .- kv2.kabc)
(+)(kv1::KVec, kv2::KVec) = KVec(kv1.k₀ .+ kv2.k₀, kv1.kabc .+ kv2.kabc)

function isapprox(kv1::KVec, kv2::KVec, cntr; kwargs...)
    k₀1, kabc1 = parts(kv1); k₀2, kabc2 = parts(kv2)
    d = length(k₀1)
    # check if k₀1 and k₀2 differ by a _primitive_ reciprocal vector
    diff = k₀1 .- k₀2
    diff = primitivebasismatrix(cntr, d)'*diff
    kbool = all(el -> isapprox(el, round(el); kwargs...), diff) 
    # check if kabc1 == kabc2; no need to check for difference by a reciprocal
    # vector, since kabc contribution is restricted to non-special values
    abcbool = isapprox(kabc1, kabc2; kwargs...)

    return kbool && abcbool
end

# Space group irreps
abstract type AbstractIrrep end
struct SGIrrep{T} <: AbstractIrrep where T
    iridx::Int64    # sequential index assigned to ir by Stokes et al
    cdml::String    # CDML label of irrep (including 𝐤-point label)
    dim::Int64      # dimensionality of irrep (i.e. size)
    sgnum::Int64    # space group number
    sglabel::String # Hermann-Mauguin label of space group
    type::Int64     # real, pseudo-real, or complex (1, 2, or 3)
    order::Int64    # number of operations
    knum::Int64     # number of 𝐤-vecs in star
    pmknum::Int64   # number of ±𝐤-vecs in star
    special::Bool   # whether star{𝐤} describes high-symmetry points
    pmkstar::Vector{KVec}       # star{𝐤} for Complex, star{±𝐤} for Real
    ops::Vector{SymOperation}   # every symmetry operation in space group
    translations::Vector{Vector{Float64}} # translations assoc with matrix repres of symops in irrep
    matrices::Vector{Matrix{T}} # non-translation assoc with matrix repres of symops in irrep
end
irreps(ir::AbstractIrrep) = ir.matrices
characters(ir::AbstractIrrep) = tr.(irreps(ir))
order(ir::AbstractIrrep) = ir.order
label(ir::AbstractIrrep) = ir.cdml
hermannmauguin(ir::AbstractIrrep) = ir.sglabel
operations(ir::AbstractIrrep) = ir.ops
isspecial(ir::AbstractIrrep) = ir.special
kstar(ir::SGIrrep) = ir.pmkstar
num(ir::AbstractIrrep) = ir.sgnum
translations(ir::AbstractIrrep) = ir.translations
type(ir::AbstractIrrep) = ir.type
klabel(ir::AbstractIrrep) = klabel(label(ir))
function klabel(label::String)
    idx = findfirst(!isletter, label)
    return label[firstindex(label):prevind(label,idx)]
end


# Little group Irreps
struct LGIrrep <: AbstractIrrep
    sgnum::Int64 # space group number
    cdml::String # CDML label of irrep (including k-point label)
    kv::KVec
    ops::Vector{SymOperation} # every symmetry operation in little group (modulo _primitive_ 𝐆)
    matrices::Vector{Matrix{ComplexF64}}
    translations::Vector{Vector{Float64}}
    type::Int64 # real, pseudo-real, or complex (⇒ 1, 2, or 3)
end
order(ir::LGIrrep) = length(operations(ir))
function irreps(ir::LGIrrep, αβγ::Union{Vector{<:Real},Nothing})
    P = ir.matrices
    τ = ir.translations
    if !iszero(τ)
        k = kvec(ir)(αβγ)
        P′ = deepcopy(P) # needs deepcopy rather than a copy due to nesting; otherwise we overwrite..!
        for (i,τ′) in enumerate(τ)
            if !iszero(τ′) && !iszero(k)
                P′[i] .*= exp(2π*im*dot(k,τ′)) # This follows the convention in Eq. (11.37) of Inui as well as the 
                                               # Bilbao server; but disagrees (as far as I can tell) with some
                                               # other references (e.g. Herring 1937a, Bilbao's _publications_?!, 
                                               # and Kovalev's book).
                                               # In those other references they have Dᵏ({I|𝐭}) = exp(-i𝐤⋅𝐭), but 
                                               # Inui has Dᵏ({I|𝐭}) = exp(i𝐤⋅𝐭) [cf. (11.36)]. The former choice 
                                               # actually appears more natural, since we usually have symmetry 
                                               # operations acting inversely on functions of spatial coordinates. 
                                               # If we swap the sign here, we probably have to swap t₀ in the check
                                               # for ray-representations in multtable(::MultTable, ::LGIrrep), to 
                                               # account for this difference. It is not enough just to swap the sign
                                               # - I checked (⇒ 112 failures in test/multtable.jl) - you would have 
                                               # to account for the fact that it would be -β⁻¹τ that appears in the 
                                               # inverse operation, not just τ. Same applies here, if you want to 
                                               # adopt the other convention, it should probably not just be a swap 
                                               # to -τ, but to -β⁻¹τ. Probably best to stick with Inui's definition.
                                               # Note that the exp(2πi𝐤⋅τ) is also the convention adopted by Stokes
                                               # et al in Eq. (1) of Acta Cryst. A69, 388 (2013), i.e. in ISOTROPY, 
                                               # so, overall, this is probably the sanest choice for this dataset.
            end
        end
        return P′
    end
    return P
end
irreps(ir::LGIrrep) = irreps(ir, nothing)
kvec(ir::LGIrrep)   = ir.kv
isspecial(ir::LGIrrep) = isspecial(kvec(ir))
issymmorph(ir::LGIrrep) = all(issymmorph.(operations(ir)))

"""
    israyrep(ir::LGIrrep, αβγ=nothing) -> (::Bool, ::Matrix)

Computes whether a given little group irrep `ir` is a ray representation 
by computing the coefficients αᵢⱼ in DᵢDⱼ=αᵢⱼDₖ; if any αᵢⱼ differ 
from unity, we consider the little group irrep a ray representation
(as opposed to the simpler "vector" representations where DᵢDⱼ=Dₖ).
The function returns a boolean (true => ray representation) and the
coefficient matrix αᵢⱼ.
"""
function israyrep(ir::LGIrrep, αβγ::Union{Nothing,Vector{Float64}}=nothing) 
    k = kvec(ir)(αβγ)
    ops = operations(ir)
    Nₒₚ = length(ops)
    α = Matrix{ComplexF64}(undef, Nₒₚ, Nₒₚ)
    mt = multtable(ops, verbose=false)
    for (row, oprow) in enumerate(ops)
        for (col, opcol) in enumerate(ops)
            t₀ = translation(oprow) + rotation(oprow)*translation(opcol) - translation(ops[mt[row,col]])
            ϕ  = 2π*k'*t₀ # include factor of 2π here due to normalized bases
            α[row,col] = exp(1im*ϕ)
        end
    end
    return (any(x->norm(x-1.0)>1e-12, α), α)
end

function show(io::IO, ::MIME"text/plain", lgir::LGIrrep)
    Nₒₚ = order(lgir)
    print(io, label(lgir))
    indent = " "^length(label(lgir))
    for (i,(op,ir)) in enumerate(zip(operations(lgir), irreps(lgir))); 
        if    i == 1; print(io, " ╮ "); 
        else          print(io, indent, " │ "); end
        print(io, xyzt(op), ":\n")
        Base.print_matrix(IOContext(io, :compact=>true), ir, indent*(i == Nₒₚ ? " ╰" : " │")*"    ")
        if i < Nₒₚ; print(io, '\n'); end
    end
end
function show(io::IO, ::MIME"text/plain", lgirvec::Union{AbstractVector{LGIrrep}, NTuple{N,LGIrrep} where N})
    print(io, "LGIrrep(#", num(lgirvec[1]), ") at ", klabel(lgirvec[1]), " = ")
    show(io,"text/plain", kvec(lgirvec[1])); println(io)
    Nᵢᵣ = length(lgirvec)
    for (i,lgir) in enumerate(lgirvec)
        show(io, "text/plain", lgir)
        if i != Nᵢᵣ; println(io); end
    end
end

function findirrep(LGIR, sgnum::Integer, cdml::String)
    kidx = findfirst(x->label(x[1])[1]==cdml[1], LGIR[sgnum])
    irrepidx = findfirst(x->label(x)==cdml, LGIR[sgnum][kidx])
    return LGIR[sgnum][kidx][irrepidx]
end


# band representations
struct BandRep
    wyckpos::String  # Wyckoff position that induces the BR
    sitesym::String  # Site-symmetry point group of Wyckoff pos (IUC notation)
    label::String    # Symbol ρ↑G, with ρ denoting the irrep of the site-symmetry group
    dim::Integer     # Dimension (i.e. # of bands) in band rep
    decomposable::Bool  # Whether a given bandrep can be decomposed further
    spinful::Bool       # Whether a given bandrep involves spinful irreps ("\bar"'ed irreps)
    irrepvec::Vector{Int64}   # Vector the references irreplabs of a parent BandRepSet; 
                              # nonzero entries correspond to an element in the band representation
    irreptags::Vector{String} # vestigial, but quite handy for display'ing; this otherwise 
                              # requires recursive data sharing between BandRep and BandRepSet
end
wyck(BR::BandRep)    = BR.wyckpos
sitesym(BR::BandRep) = BR.sitesym
label(BR::BandRep)   = BR.label
dim(BR::BandRep)     = BR.dim
humanreadable(BR::BandRep) = BR.irreptags
vec(BR::BandRep)     = BR.irrepvec
function show(io::IO, ::MIME"text/plain", BR::BandRep)
    print(label(BR), " (", dim(BR), "): [")
    join(io, map(x->replace(x, '⊕'=>'+'), humanreadable(BR)), ", ") # ⊕ doesn't render well in my terminal; swap for ordinary plus
    print(io, "]")
end

struct BandRepSet
    sgnum::Integer          # space group number, sequential
    bandreps::Vector{BandRep}
    kvs::Vector{KVec}       # Vector of 𝐤-points
    klabs::Vector{String}   # Vector of associated 𝐤-labels (in CDML notation)
    irreplabs::Vector{String} # Vector of (sorted) CDML irrep labels at _all_ 𝐤-points
    allpaths::Bool          # Whether all paths (true) or only maximal 𝐤-points (false) are included
    spinful::Bool           # Whether the band rep set includes (true) or excludes (false) spinful irreps
end
num(BRS::BandRepSet)    = BRS.sgnum
klabels(BRS::BandRepSet) = BRS.klabs
kvecs(BRS::BandRepSet)  = BRS.kvs
hasnonmax(BRS::BandRepSet) = BRS.allpaths
irreplabels(BRS::BandRepSet)   = BRS.irreplabs
isspinful(BRS::BandRepSet) = BRS.spinful
reps(BRS::BandRepSet)   = BRS.bandreps
length(BRS::BandRepSet) = length(reps(BRS))
getindex(BRS::BandRepSet, keys...) = reps(BRS)[keys...]
lastindex(BRS::BandRepSet, d::Int64) = length(BRS)


function show(io::IO, ::MIME"text/plain", BRS::BandRepSet)
    Nirreps = length(irreplabels(BRS))
    println(io, "BandRepSet (#$(num(BRS))):")
    println(io, "k-vecs ($(hasnonmax(BRS) ? "incl. non-maximal" : "maximal only")):")
    for (lab,kv) in zip(klabels(BRS), kvecs(BRS))
        print(io,"   ", lab, ": "); show(io, "text/plain", kv); println(io)
    end

    # prep-work
    maxlen = maximum(x->length(label(x))+ndigits(dim(x)), reps(BRS))+3
    threshold = 30
    if Nirreps > threshold
        toomuch = div((Nirreps-threshold+2),2)
        midpoint = div(Nirreps, 2)
        skiprange = (-toomuch:toomuch) .+ midpoint
        abbreviate = true
    else
        abbreviate = false
    end
    # "title"
    println(io, "$(length(BRS)) band representations", 
                " ($(isspinful(BRS) ? "spinful" : "spinless"))",
                " sampling $(Nirreps) irreps:")
    print(io, " "^(maxlen+5),'║'); # align with spaces
    for (j,lab) in enumerate(irreplabels(BRS)) # irrep labels
        if abbreviate && j∈skiprange
            if j == first(skiprange)
                print(io, "\b  …  ")
            end
        else
            print(io, ' ', lab, j != Nirreps ? " │" : " ║")
        end
    end
    println(io)
    for (i,BR) in enumerate(reps(BRS))
        print(io, "   ", label(BR), " (", dim(BR), "):",                      # bandrep label
                  " "^(maxlen-length(label(BR))-ndigits((dim(BR)))-2), '║')
        for (j,v) in enumerate(vec(BR)) # vector representation of band rep
            if abbreviate && j∈skiprange
                if j == first(skiprange)
                    print(io, mod(i,4) == 0 ? "\b  …  " : "\b     ")
                end
            else
                print(io, "  ")
                !iszero(v) ? print(io, v) : print(io, '·')
                print(io, " "^(length(irreplabels(BRS)[j])-1)) # assumes we will never have ndigit(v) != 1
                print(io, j != Nirreps ? '│' : '║')
            end
        end
        if i != length(BRS); println(io); end
    end
end