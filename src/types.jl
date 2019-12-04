import Base: show

# Crystalline lattice
struct Crystal{D}
    Rs::NTuple{D,Vector{Float64}}
end
Crystal(Rs) = Crystal{length(Rs)}(Rs)
Crystal(Rs...) = Crystal(Rs)
basis(C::Crystal) = C.Rs
dim(C::Crystal{D}) where D = D
function show(io::IO, ::MIME"text/plain", C::Crystal)
    print(io, "$(dim(C))D Crystal:")
    print(io, " ($(crystalsystem(C)))");
    for (i,R) in enumerate(basis(C))
        print(io, "\n   R$(i): "); print(io, R); 
    end
end
norms(C::Crystal) = norm.(basis(C))
_angle(rA,rB) = acos(dot(rA,rB)/(norm(rA)*norm(rB)))
function angles(C::Crystal{D}) where D
    R = basis(C)
    γ = _angle(R[1], R[2])
    if D == 3
        α = _angle(R[2], R[3])
        β = _angle(R[3], R[1])
        return α,β,γ
    end
    return γ
end


# Symmetry operations
struct SymOperation
    xyzt::String
    matrix::Matrix{Float64}
end
SymOperation(s::AbstractString) = SymOperation(string(s), xyzt2matrix(s))
SymOperation(m::Matrix{<:Real}) = SymOperation(matrix2xyzt(m), float(m))
matrix(op::SymOperation) = op.matrix
xyzt(op::SymOperation) = op.xyzt
dim(op::SymOperation) = size(matrix(op),1)
getindex(op::SymOperation, keys...) = matrix(op)[keys...]   # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
lastindex(op::SymOperation, d::Int64) = size(matrix(op), d) # allows using `end` in indices
rotation(m::Matrix{Float64}) = @view m[:,1:end-1] # rotational (proper or improper) part of an operation
rotation(op::SymOperation)   = rotation(matrix(op))
translation(m::Matrix{Float64}) = @view m[:,end]  # translation part of an operation
translation(op::SymOperation)   = translation(matrix(op))
(==)(op1::SymOperation, op2::SymOperation) = (xyzt(op1) == xyzt(op2)) && (matrix(op1) == matrix(op2))
isapprox(op1::SymOperation, op2::SymOperation; kwargs...) = isapprox(matrix(op1), matrix(op2); kwargs...)
function show(io::IO, ::MIME"text/plain", op::SymOperation)
    print(io, seitz(op),":\n   (", xyzt(op), ")\n")
    Base.print_matrix(IOContext(io, :compact=>true), op.matrix, "   ")
end
function show(io::IO, ::MIME"text/plain", ops::AbstractVector{<:SymOperation})
    for (i,op) in enumerate(ops)
        show(io, "text/plain", op)
        if i < length(ops); print(io, "\n"); end
    end
end

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
function (kv::KVec)(αβγ::AbstractVector{<:Real})
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
    write(buf, '[')
    if isspecial(kv)
        for i in eachindex(k₀) 
            coord = k₀[i] == -0.0 ? 0.0 : k₀[i] # normalize -0.0 to 0.0
            print(buf, coord)
            # prepare for next coordinate/termination
            i == length(k₀) ? write(buf, ']') : write(buf, ", ")
        end
    else
        for i in eachindex(k₀)
            # fixed parts
            if !iszero(k₀[i]) || iszero(@view kabc[i,:]) # don't print zero, if it adds unto anything nonzero
                coord = k₀[i] == -0.0 ? 0.0 : k₀[i] # normalize -0.0 to 0.0
                print(buf, coord)
            end
            # free-parameter parts
            for j in eachindex(k₀) 
                if !iszero(kabc[i,j])
                    sgn = signaschar(kabc[i,j])
                    if !(iszero(k₀[i]) && sgn=='+' && iszero(kabc[i,1:j-1])) # don't print '+' if nothing precedes it
                        write(buf, sgn)
                    end
                    if abs(kabc[i,j]) != oneunit(eltype(kabc)) # don't print prefactors of 1
                        print(buf, abs(kabc[i,j]))
                    end
                    write(buf, j==1 ? 'α' : (j == 2 ? 'β' : 'γ'))
                end
            end
            # prepare for next coordinate/termination
            i == length(k₀) ? write(buf, ']') : write(buf, ", ")
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

Fractions such as `1/2` can be parsed: but use of any other special operator
besides `/` will result in faulty operations (e.g. do not use `*`).
"""
function KVec(str::AbstractString)
    str = filter(!isspace, strip(str, ['(',')','[',']'])) # tidy up string (remove parens & spaces)
    xyz = split(str,',')
    dim = length(xyz)
    k₀ = zeros(Float64, dim); kabc = zeros(Float64, dim, dim)
    for (i, coord) in enumerate(xyz)
        # --- "free" coordinates, kabc[i,:] ---
        for (j, matchgroup) in enumerate((('α','u'),('β','v'),('γ','w')))
            pos₂ = findfirst(∈(matchgroup), coord)
            if !isnothing(pos₂)
                match = searchpriornumerals(coord, pos₂)
                kabc[i,j] = parse(Float64, match)
            end
        end
        
        # --- "fixed" coordinate, k₀[i] ---
        sepidx′ = findfirst(r"\b(\+|\-)", coord) # find any +/- separators between fixed and free parts
        # regex matches '+' or '-', except if they are the first character in 
        # string (or if they are preceded by space; but that cannot occur here)   
        if sepidx′===nothing # no separators
            if last(coord) ∈ ('α','u','β','v','γ','w') # free-part only case
                continue # k₀[i] is zero already
            else                                       # constant-part only case
                k₀[i] = parsefraction(coord)
            end
        else # exploit that we require fixed parts to come before free parts
            k₀[i] = parsefraction(coord[firstindex(coord):prevind(coord, first(sepidx′))])
        end
    end
    return KVec(k₀, kabc)
end

# arithmetic with k-vectors
(-)(kv::KVec) = KVec(.- kv.k₀, .- kv.kabc)
(-)(kv1::KVec, kv2::KVec) = KVec(kv1.k₀ .- kv2.k₀, kv1.kabc .- kv2.kabc)
(+)(kv1::KVec, kv2::KVec) = KVec(kv1.k₀ .+ kv2.k₀, kv1.kabc .+ kv2.kabc)

"""
    isapprox(kv1::KVec, kv2::KVec[, cntr::Char]; kwargs...) --> Bool
                                                            
Compute approximate equality of two KVec's `k1` and `k2` modulo any 
primitive G-vectors. To ensure that primitive G-vectors are used, 
the centering type `cntr` (see `centering(cntr, dim)`) must be given
(the dimensionality is inferred from `kv1` and `kv2`).
Optionally, keyword arguments (e.g., `atol` and `rtol`) can be 
provided, to include in calls to `Base.isapprox`.

If `cntr` is not provided, the comparison will not account for equivalence
by primitive G-vectors.
"""
function isapprox(kv1::KVec, kv2::KVec, cntr::Char; kwargs...)
    k₀1, kabc1 = parts(kv1); k₀2, kabc2 = parts(kv2)  # ... unpacking

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
# ... without considerations of G-vectors
function isapprox(kv1::KVec, kv2::KVec; kwargs...) 
    k₀1, kabc1 = parts(kv1); k₀2, kabc2 = parts(kv2)  # ... unpacking
       
    return isapprox(k₀1, k₀2; kwargs...) && isapprox(kabc1, kabc2; kwargs...)
end

function (==)(kv1::KVec, kv2::KVec)   
    k₀1, kabc1 = parts(kv1); k₀2, kabc2 = parts(kv2)  # ... unpacking
       
    return k₀1 == k₀2 && kabc1 == kabc2
end



# Abstract spatial group
abstract type AbstractGroup end
num(g::AbstractGroup) = g.num
operations(g::AbstractGroup) = g.operations
dim(g::AbstractGroup) = g.dim
order(g::AbstractGroup) = length(operations(g))
getindex(g::AbstractGroup, keys...) = operations(g)[keys...]    # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
lastindex(g::AbstractGroup, d::Int64) = size(operations(g), d)  # allows using `end` in indices

function show(io::IO, ::MIME"text/plain", g::T) where T<:AbstractGroup
    if isa(g, SpaceGroup)
        groupprefix = dim(g) == 3 ? "Space group" : (dim(g) == 2 ? "Plane group" : nothing)
    elseif isa(g, PointGroup)
        groupprefix = "Point group"
    else
        groupprefix = string(T)
    end
    println(io, groupprefix, " #", num(g), " (", label(g), ") with ", order(g), " operations:")
    show(io, "text/plain", operations(g))
end
function show(io::IO, ::MIME"text/plain", gs::AbstractVector{<:AbstractGroup})
    Ngs = length(gs)
    for (i,g) in enumerate(gs); 
        show(io, "text/plain", g); 
        if i < Ngs; print(io, '\n'); end
    end
end

# Space group
struct SpaceGroup <: AbstractGroup
    num::Int64
    operations::Vector{SymOperation}
    dim::Int64
end
label(sg::SpaceGroup) = iuc(num(sg), dim(sg))

# Point group
struct PointGroup <: AbstractGroup
    num::Int64
    label::String
    operations::Vector{SymOperation}
    dim::Int64
end
label(pg::PointGroup) = pg.label

# Little group
struct LittleGroup{D} <: AbstractGroup
    num::Int64
    kv::KVec
    klab::String
    operations::Vector{SymOperation}
end
LittleGroup(num::Int64, kv::KVec, klab::String, ops::AbstractVector{SymOperation}) = LittleGroup{dim(kv)}(num, kv, klab, ops)
dim(lg::LittleGroup{D}) where D = D
kvec(lg::LittleGroup) = lg.kv
label(lg::LittleGroup)  = iuc(num(lg), dim(lg))*" at "*string(kvec(lg))


# Abstract group irreps
""" 
    AbstractIrrep (abstract type)

Abstract supertype for irreps: must have fields cdml, matrices, 
and type (and possibly translations). Must implement a function
`irreps` that return the associated irrep matrices.
"""
abstract type AbstractIrrep end
label(ir::AbstractIrrep) = ir.cdml
matrices(ir::AbstractIrrep) = ir.matrices    
type(ir::AbstractIrrep) = ir.type
translations(ir::T) where T<:AbstractIrrep = hasfield(T, :translations) ? ir.translations : nothing
characters(ir::AbstractIrrep) = tr.(irreps(ir))
klabel(ir::AbstractIrrep) = klabel(label(ir))
irdim(ir::AbstractIrrep)  = size(first(matrices(ir)),1)

# 3D Space group irreps
struct SGIrrep{T} <: AbstractIrrep where T
    iridx::Int64    # sequential index assigned to ir by Stokes et al
    cdml::String    # CDML label of irrep (including 𝐤-point label)
    irdim::Int64    # dimensionality of irrep (i.e. size)
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
num(sgir::SGIrrep) = sgir.sgnum
irreps(sgir::SGIrrep) = sgir.matrices
order(sgir::SGIrrep) = sgir.order
hermannmauguin(sgir::SGIrrep) = sgir.sglabel
operations(sgir::SGIrrep) = sgir.ops
isspecial(sgir::SGIrrep) = sgir.special
kstar(sgir::SGIrrep) = sgir.pmkstar
irdim(sgir::SGIrrep) = sgir.irdim
dim(sgir::SGIrrep) = 3

# Little group irreps
struct LGIrrep{D} <: AbstractIrrep
    cdml::String # CDML label of irrep (including k-point label)
    lg::LittleGroup{D} # contains sgnum, kvec, klab, and operations that define the little group (and dimension as type parameter)
    matrices::Vector{Matrix{ComplexF64}}
    translations::Vector{Vector{Float64}}
    type::Int64 # real, pseudo-real, or complex (⇒ 1, 2, or 3)
end
function LGIrrep{D}(cdml::String, lg::LittleGroup{D}, 
                    matrices::Vector{Matrix{ComplexF64}}, 
                    translations_sentinel::Nothing, # sentinel value for all-zero translations
                    type::Int64) where D
    translations = [zeros(Float64,D) for _=Base.OneTo(order(lg))]
    return LGIrrep{D}(cdml, lg, matrices, translations, type)
end
littlegroup(lgir::LGIrrep) = lgir.lg
num(lgir::LGIrrep) = num(littlegroup(lgir))
dim(lgir::LGIrrep{D}) where D = D
operations(lgir::LGIrrep)  = operations(littlegroup(lgir))
order(lgir::LGIrrep) = order(littlegroup(lgir))
kvec(lgir::LGIrrep)  = kvec(littlegroup(lgir))
isspecial(lgir::LGIrrep)  = isspecial(kvec(lgir))
issymmorph(lgir::LGIrrep) = issymmorph(littlegroup(lgir))
kstar(lgir::LGIrrep) = kstar(get_sgops(num(lgir), dim(lgir)), kvec(lgir), centering(num(lgir), dim(lgir)))
function irreps(lgir::LGIrrep, αβγ::Union{Vector{<:Real},Nothing}=nothing)
    P = lgir.matrices
    τ = lgir.translations
    if !iszero(τ)
        k = kvec(lgir)(αβγ)
        P′ = deepcopy(P) # needs deepcopy rather than a copy due to nesting; otherwise we overwrite..!
        for (i,τ′) in enumerate(τ)
            if !iszero(τ′) && !iszero(k)
                P′[i] .*= cis(2π*dot(k,τ′)) # This follows the convention in Eq. (11.37) of Inui as well as the 
                # note cis(x) = exp(ix)     # Bilbao server; but disagrees (as far as I can tell) with some
                                            # other references (e.g. Herring 1937a, Bilbao's _publications_?!, 
                                            # and Kovalev's book).
                                            # In those other references they have Dᵏ({I|𝐭}) = exp(-i𝐤⋅𝐭), but 
                                            # Inui has Dᵏ({I|𝐭}) = exp(i𝐤⋅𝐭) [cf. (11.36)]. The former choice 
                                            # actually appears more natural, since we usually have symmetry 
                                            # operations acting inversely on functions of spatial coordinates. 
                                            # If we swap the sign here, we probably have to swap t₀ in the check
                                            # for ray-representations in multtable(::MultTable, ::LGIrrep), to 
                                            # account for this difference. It is not enough just to swap the sign
                                            # - I checked (⇒ 172 failures in test/multtable.jl) - you would have 
                                            # to account for the fact that it would be -β⁻¹τ that appears in the 
                                            # inverse operation, not just τ. Same applies here, if you want to 
                                            # adopt the other convention, it should probably not just be a swap 
                                            # to -τ, but to -β⁻¹τ. Probably best to stick with Inui's definition.
                                            # Note that the exp(2πi𝐤⋅τ) is also the convention adopted by Stokes
                                            # et al in Eq. (1) of Acta Cryst. A69, 388 (2013), i.e. in ISOTROPY 
                                            # (also expliciated at https://stokes.byu.edu/iso/irtableshelp.php),
                                            # so, overall, this is probably the sanest choice for this dataset.
            end
        end
        return P′
    end
    return P
end

function klabel(cdml::String)
    idx = findfirst(c->isdigit(c) || issubdigit(c), cdml) # look for regular digit or subscript digit
    return cdml[firstindex(cdml):prevind(cdml,idx)]
end


"""
    israyrep(lgir::LGIrrep, αβγ=nothing) -> (::Bool, ::Matrix)

Computes whether a given little group irrep `ir` is a ray representation 
by computing the coefficients αᵢⱼ in DᵢDⱼ=αᵢⱼDₖ; if any αᵢⱼ differ 
from unity, we consider the little group irrep a ray representation
(as opposed to the simpler "vector" representations where DᵢDⱼ=Dₖ).
The function returns a boolean (true => ray representation) and the
coefficient matrix αᵢⱼ.
"""
function israyrep(lgir::LGIrrep, αβγ::Union{Nothing,Vector{Float64}}=nothing) 
    k = kvec(lgir)(αβγ)
    ops = operations(lgir)
    Nₒₚ = length(ops)
    α = Matrix{ComplexF64}(undef, Nₒₚ, Nₒₚ)
    # TODO: Verify that this is OK; not sure if we can just use the primitive basis 
    #       here, given the tricks we then perform subsequently?
    mt = multtable(primitivize.(ops, centering(num(lgir))), verbose=false) 
    for (row, oprow) in enumerate(ops)
        for (col, opcol) in enumerate(ops)
            t₀ = translation(oprow) + rotation(oprow)*translation(opcol) - translation(ops[mt[row,col]])
            ϕ  = 2π*dot(k,t₀) # include factor of 2π here due to normalized bases
            α[row,col] = cis(ϕ)
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