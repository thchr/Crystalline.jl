# --- DirectBasis and ReciprocalBasis for crystalline lattices ---
abstract type Basis{D} <: AbstractVector{SVector{D,Float64}} end
for T in (:DirectBasis, :ReciprocalBasis)
    @eval struct $T{D} <: Basis{D}
              vecs::NTuple{D,SVector{D,Float64}}
          end
    @eval $T(Rs::NTuple{D,AbstractVector{<:Real}}) where D = $T{D}(SVector{D,Float64}.(Rs))
    @eval $T(Rs::NTuple{D,NTuple{D,<:Real}}) where D = $T{D}(SVector{D,Float64}.(Rs))
    @eval $T(Rs::AbstractVector{<:Real}...) = $T(Rs)
    @eval $T(Rs::NTuple{D,<:Real}...) where D = $T{D}(SVector{D,Float64}.(Rs))
end

vecs(Vs::Basis) = Vs.vecs
# define the AbstractArray interface for DirectBasis{D}
getindex(Vs::Basis, i::Int) = vecs(Vs)[i] 
firstindex(::Basis) = 1
lastindex(::Basis{D}) where D = D
setindex!(Vs::Basis, vec::Vector{Float64}, i::Int) = (Vs[i] .= vec)
size(::Basis{D}) where D = (D,)
IndexStyle(::Basis) = IndexLinear()

norms(Rs::Basis) = norm.(Rs)
_angle(rA,rB) = acos(dot(rA,rB)/(norm(rA)*norm(rB)))
function angles(Rs::Basis{D}) where D
    D == 1 && return nothing
    γ = _angle(Rs[1], Rs[2])
    if D == 3
        α = _angle(Rs[2], Rs[3])
        β = _angle(Rs[3], Rs[1])
        return α,β,γ
    end
    return γ
end
"""
    basis2matrix(Vs::Basis{D}) where D

Compute a matrix `[Vs[1] Vs[2] .. Vs[D]]` from `Vs::Basis{D}`, i.e. a matrix whose columns
are the basis vectors in `Vs`. 

Note: Trying to use the iteration interface via `hcat(Vs...)` does not lead to a correctly
      inferred type Matrix::Float64 (and a type-assertion does not improve speed much).
      Instead, we just use the .vec field of `Vs` directly, which achieves good performance.
"""
basis2matrix(Vs::Basis{D}) where D = hcat(vecs(Vs)...)


# --- Symmetry operations ---
struct SymOperation{D} <: AbstractMatrix{Float64}
    matrix::Matrix{Float64}
    # It doesn't seem to be possible to convert `matrix` from ::Matrix{Float64} to
    # ::SMatrix{D,D+1,Float64,D*(D+1)} in a nice way, as it is currently impossible to do
    # computations on type parameters in type definitions, as discussed e.g. in 
    #   https://github.com/JuliaLang/julia/issues/18466 
    #   https://discourse.julialang.org/t/addition-to-parameter-of-parametric-type/20059
    # It doesn't help to split `matrix` into a point `R` and a translation part `τ`, since
    # declaring `R::SMatrix{D,D,Float64,D*D}` also isn't possible; so we'd need to include
    # a useless `L` type, which would be forced to equal `D*D` in the struct constructor.
    # Overall, it doesn't seem worth it at this point: could maybe be done for Julia 2.0.
    # TODO: Splitting `matrix` into a point-group/rotation and a translation part would 
    # probably be worthwhile though, since we only every really deal with them seperately.
end
SymOperation{D}(s::AbstractString) where D = (m=xyzt2matrix(s); SymOperation{D}(m))
# type-unstable convenience constructors; avoid for anything non-REPL related, if possible
SymOperation(m::Matrix{<:Real}) = SymOperation{size(m,1)}(float(m))   
SymOperation(s::AbstractString) = (m=xyzt2matrix(s); SymOperation(m)) 

matrix(op::SymOperation) = op.matrix
xyzt(op::SymOperation) = matrix2xyzt(matrix(op))
dim(::SymOperation{D}) where D = D
# define the AbstractArray interface for SymOperation
getindex(op::SymOperation, keys...) = matrix(op)[keys...]
firstindex(::SymOperation) = 1
lastindex(op::SymOperation{D}) where D = D*(D+1)
lastindex(op::SymOperation{D}, d::Int64) where D = d == 1 ? D : (d == 2 ? D+1 : 1)
IndexStyle(::SymOperation) = IndexLinear()
size(::SymOperation{D}) where D = (D,D+1)
eltype(::SymOperation) = Float64

rotation(m::Matrix{<:Real}) = @view m[:,1:end-1] # rotational (proper or improper) part of an operation
rotation(op::SymOperation)  = rotation(matrix(op))
translation(m::Matrix{<:Real}) = @view m[:,end]  # translation part of an operation
translation(op::SymOperation)  = translation(matrix(op))
(==)(op1::SymOperation, op2::SymOperation) = (dim(op1) == dim(op2) && xyzt(op1) == xyzt(op2)) && (matrix(op1) == matrix(op2))
isapprox(op1::SymOperation, op2::SymOperation; kwargs...)= (dim(op1) == dim(op2) && isapprox(matrix(op1), matrix(op2); kwargs...))
unpack(op::SymOperation) = (rotation(op), translation(op))

# --- Multiplication table ---
struct MultTable{D} <: AbstractMatrix{Int64}
    operations::Vector{SymOperation{D}}
    indices::Matrix{Int64}
    isgroup::Bool
end
indices(mt::MultTable) = mt.indices
isgroup(mt::MultTable) = mt.isgroup
getindex(mt::MultTable, keys...) = indices(mt)[keys...]
firstindex(mt::MultTable, d) = 1
lastindex(mt::MultTable, d::Int64) = size(indices(mt),d)

# --- 𝐤-vectors ---
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
KVec(k₀s::T...) where T<:Real = KVec([float.(k₀s)...])
parts(kv::KVec) = (kv.k₀, kv.kabc)
dim(kv::KVec) = length(kv.k₀)
isspecial(kv::KVec) = iszero(kv.kabc)
# returns a vector whose entries are true (false) if α,β,γ, respectively, are free parameters (not featured) in `kv`
freeparams(kv::KVec)  = map(j->!iszero(@view kv.kabc[:,j]), Base.OneTo(dim(kv))) 
nfreeparams(kv::KVec) = count(j->!iszero(@view kv.kabc[:,j]), Base.OneTo(dim(kv))) # total number of free parameters in `kv`
function (kv::KVec)(αβγ::AbstractVector{<:Real})
    k₀, kabc = parts(kv)
    return k₀ + kabc*αβγ
end
(kv::KVec)(αβγ::Vararg{<:Real, 2}) = kv([αβγ[1], αβγ[2]])
(kv::KVec)(αβγ::Vararg{<:Real, 3}) = kv([αβγ[1], αβγ[2], αβγ[3]])
(kv::KVec)() = kv.k₀
(kv::KVec)(::Nothing) = kv.k₀

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
                # TODO: Seems broken for k-planes like KVec("α+β,β,0") due to the α+β term
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
zero(kv::KVec) = KVec(zero(kv.k₀))

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

# --- Abstract spatial group ---
abstract type AbstractGroup{D} <: AbstractVector{SymOperation{D}} end
num(g::AbstractGroup) = g.num
operations(g::AbstractGroup) = g.operations
dim(g::AbstractGroup{D}) where D = D
# define the AbstractArray interface for AbstractGroup
getindex(g::AbstractGroup, keys...) = operations(g)[keys...]    # allows direct indexing into an op::SymOperation like op[1,2] to get matrix(op)[1,2]
firstindex(::AbstractGroup) = 1
lastindex(g::AbstractGroup, d::Int64) = size(operations(g), d)  # allows using `end` in indices
setindex!(g::AbstractGroup, op::SymOperation, i::Int) = (operations(g)[i] .= op)
size(g::AbstractGroup) = (length(operations(g)),)
IndexStyle(::AbstractGroup) = IndexLinear()
eltype(::AbstractGroup{D}) where D = SymOperation{D}
order(g::AbstractGroup) = length(g)

# --- Space group ---
struct SpaceGroup{D} <: AbstractGroup{D}
    num::Int64
    operations::Vector{SymOperation{D}}
end
label(sg::SpaceGroup) = iuc(sg)

# --- Point group ---
struct PointGroup{D} <: AbstractGroup{D}
    num::Int64
    label::String
    operations::Vector{SymOperation{D}}
end
label(pg::PointGroup) = pg.label
iuc(pg::PointGroup) = label(pg)

# --- Little group ---
struct LittleGroup{D} <: AbstractGroup{D}
    num::Int64
    kv::KVec
    klab::String
    operations::Vector{SymOperation{D}}
end
LittleGroup(num::Int64, kv::KVec, klab::String, ops::AbstractVector{SymOperation{D}}) where D = LittleGroup{D}(num, kv, klab, ops)
LittleGroup(num::Int64, kv::KVec, ops::AbstractVector{SymOperation{D}}) where D = LittleGroup{D}(num, kv, "", ops)
kvec(lg::LittleGroup) = lg.kv
klabel(lg::LittleGroup) = lg.klab
label(lg::LittleGroup)  = iuc(num(lg), dim(lg))*" at "*klabel(lg)*" = "*string(kvec(lg))

# --- Abstract group irreps ---
""" 
    AbstractIrrep{D} (abstract type)

Abstract supertype for irreps of dimensionality `D`: must have fields `cdml`, `matrices`,
and `type` (and possibly `translations`). Must implement a function `irreps` that returns
the associated irrep matrices.
"""
abstract type AbstractIrrep{D} end
label(ir::AbstractIrrep) = ir.cdml
matrices(ir::AbstractIrrep) = ir.matrices    
type(ir::AbstractIrrep) = ir.type
translations(ir::T) where T<:AbstractIrrep = hasfield(T, :translations) ? ir.translations : nothing
characters(ir::AbstractIrrep, αβγ::Union{AbstractVector{<:Real},Nothing}=nothing) = tr.(irreps(ir, αβγ))
irdim(ir::AbstractIrrep)  = size(first(matrices(ir)),1)
klabel(ir::AbstractIrrep) = klabel(label(ir))
order(ir::AbstractIrrep)  = order(group(ir))
operations(ir::AbstractIrrep) = operations(group(ir))
num(ir::AbstractIrrep) = num(group(ir))
dim(ir::AbstractIrrep{D}) where D = D
function klabel(cdml::String)
    idx = findfirst(c->isdigit(c) || issubdigit(c), cdml) # look for regular digit or subscript digit
    previdx = idx !== nothing ? prevind(cdml, idx) : lastindex(cdml)
    return cdml[firstindex(cdml):previdx]
end

# --- Point group irreps ---
struct PGIrrep{D} <: AbstractIrrep{D}
    cdml::String
    pg::PointGroup{D}
    matrices::Vector{Matrix{ComplexF64}}
    type::Int64
end
irreps(pgir::PGIrrep, αβγ::Nothing=nothing) = pgir.matrices
group(pgir::PGIrrep) = pgir.pg

# printing
function prettyprint_irrep_matrix(io::IO, pgir::PGIrrep, i::Integer, prefix::AbstractString)
    P = pgir.matrices[i]
    prettyprint_scalar_or_matrix(io, P, prefix, false)
end

# --- Little group irreps ---
struct LGIrrep{D} <: AbstractIrrep{D}
    cdml::String # CDML label of irrep (including k-point label)
    lg::LittleGroup{D} # contains sgnum, kvec, klab, and operations that define the little group (and dimension as type parameter)
    matrices::Vector{Matrix{ComplexF64}}
    translations::Vector{Vector{Float64}}
    type::Int64 # real, pseudo-real, or complex (⇒ 1, 2, or 3)
    iscorep::Bool # Whether this irrep really represents a corep (only relevant for `type`s 2 and 3; leads to special handling for `irreps(..)` and printing)
end
function LGIrrep{D}(cdml::String, lg::LittleGroup{D}, 
                    matrices::Vector{Matrix{ComplexF64}}, 
                    translations::Vector{Vector{Float64}},
                    type::Int64) where D
    return LGIrrep{D}(cdml, lg, matrices, translations, type, false)
end
function LGIrrep{D}(cdml::String, lg::LittleGroup{D}, 
                    matrices::Vector{Matrix{ComplexF64}}, 
                    translations_sentinel::Nothing, # sentinel value for all-zero translations
                    type::Int64) where D
    translations = [zeros(Float64,D) for _=Base.OneTo(order(lg))]
    return LGIrrep{D}(cdml, lg, matrices, translations, type)
end
group(lgir::LGIrrep) = lgir.lg
iscorep(lgir::LGIrrep) = lgir.iscorep
kvec(lgir::LGIrrep)  = kvec(group(lgir))
isspecial(lgir::LGIrrep)  = isspecial(kvec(lgir))
issymmorph(lgir::LGIrrep) = issymmorph(group(lgir))
kstar(lgir::LGIrrep) = kstar(spacegroup(num(lgir), dim(lgir)), 
                             kvec(lgir), centering(num(lgir), dim(lgir)))
function irreps(lgir::LGIrrep, αβγ::Union{Vector{<:Real},Nothing}=nothing)
    P = lgir.matrices
    τ = lgir.translations
    if !iszero(τ)
        k = kvec(lgir)(αβγ)
        P = deepcopy(P) # needs deepcopy rather than a copy due to nesting; otherwise we overwrite..!
        for (i,τ′) in enumerate(τ)
            if !iszero(τ′) && !iszero(k)
                P[i] .*= cis(2π*dot(k,τ′)) # This follows the convention in Eq. (11.37) of Inui as well as the 
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
    end

    if iscorep(lgir)
        t = type(lgir) 
        if t == 2 # Pseudo-real (doubles)
            return _blockdiag2x2.(P)
        elseif t == 3 # Complex (conj-doubles)
            return _blockdiag2x2_conj.(P)
        else
            throw(DomainError(type, "Unexpected combination of iscorep=true and type≠{2,3}"))
        end
    else
        return P
    end
    return P
end

function _blockdiag2x2(A::Matrix{T}) where T
    n = LinearAlgebra.checksquare(A)
    B = zeros(T, 2*n, 2*n)
    @inbounds for I in 0:1
        I′ = I*n
        for i in Base.OneTo(n)
            i′ = I′+i
            for j in Base.OneTo(n)
                B[i′,I′+j] = A[i,j]
            end
        end
    end
    return B
end
function _blockdiag2x2_conj(A::Matrix{T}) where T
    n = LinearAlgebra.checksquare(A)
    B = zeros(T, 2*n, 2*n)
    @inbounds for i in Base.OneTo(n) # upper left block
        for j in Base.OneTo(n)
            B[i,j] = A[i,j]
        end
    end
    @inbounds for i in Base.OneTo(n) # lower right block
        i′ = n+i
        for j in Base.OneTo(n)
            B[i′,n+j] = conj(A[i,j])
        end
    end
    return B
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
    return any(x->norm(x-1.0)>DEFAULT_ATOL, α), α
end

function find_lgirreps(lgirsvec::AbstractVector{<:AbstractVector{<:LGIrrep}}, klab::String, verbose::Bool=false)
    kidx = findfirst(x->klabel(first(x))==klab, lgirsvec)
    if kidx === nothing
        if verbose
            println("Didn't find any matching k-label in lgirsvec: "*
                    "the label may be specified incorrectly, or the irrep is missing "*
                    "(e.g. the irrep could be a axes-dependent irrep)")
            @info klab klabel.(first.(lgirsvec))
        end
        return nothing 
    else
        return lgirsvec[kidx] # return an "lgirs" (vector of `LGIrrep`s)
    end
end
find_lgirreps(sgnum::Integer, klab::String, Dᵛ::Val{D}) where D = find_lgirreps(get_lgirreps(sgnum, Dᵛ), klab)
find_lgirreps(sgnum::Integer, klab::String, D::Integer=3) = find_lgirreps(sgnum, klab, Val(D))



# --- Character table ---
struct CharacterTable{D}
    ops::Vector{SymOperation{D}}
    irlabs::Vector{String}
    chartable::Matrix{ComplexF64} # Stored as irreps-along-columns & operations-along-rows
    # TODO: for LGIrreps, it might be nice to keep this more versatile and include the 
    #       translations and kvec as well; then we could print a result that doesn't  
    #       specialize on a given αβγ choice (see also CharacterTable(::LGirrep))
    tag::String
end
CharacterTable{D}(ops::AbstractVector{SymOperation{D}}, 
                  irlabs::Vector{String}, 
                  chartable::Matrix{ComplexF64}) where D = CharacterTable{D}(ops, irlabs, chartable, "")
operations(ct::CharacterTable) = ct.ops
labels(ct::CharacterTable) = ct.irlabs
characters(ct::CharacterTable) = ct.chartable
tag(ct::CharacterTable) = ct.tag

"""
    CharacterTable(irs::AbstractVector{<:AbstractIrrep}, αβγ=nothing)

Returns a `CharacterTable` associated with vector of `AbstractIrrep`s `irs`. 

Optionally, an `αβγ::AbstractVector{<:Real}` variable can be passed to evaluate the irrep
(and associated characters) with concrete free parameters (e.g., for `LGIrrep`s, a concrete
k-vector sampled from a "line-irrep"). Defaults to `nothing`, indicating it being either 
irrelevant (e.g., for `PGIrrep`s) or all free parameters implicitly set to zero.
"""
function CharacterTable(irs::AbstractVector{<:AbstractIrrep{D}},
                        αβγ::Union{AbstractVector{<:Real}, Nothing}=nothing) where D
    table = Array{ComplexF64}(undef, order(first(irs)), length(irs))
    for (i,col) in enumerate(eachcol(table))
        col .= characters(irs[i], αβγ)
    end
    g = group(first(irs))
    tag = "#"*string(num(g))*" ("*label(g)*")"
    return CharacterTable{D}(operations(first(irs)), label.(irs), table, tag)
end

# --- Band representations ---
struct BandRep <: AbstractVector{Int64}
    wyckpos::String  # Wyckoff position that induces the BR
    sitesym::String  # Site-symmetry point group of Wyckoff pos (IUC notation)
    label::String    # Symbol ρ↑G, with ρ denoting the irrep of the site-symmetry group
    dim::Integer     # Dimension (i.e. # of bands) in band rep
    decomposable::Bool  # Whether a given bandrep can be decomposed further
    spinful::Bool       # Whether a given bandrep involves spinful irreps ("\bar"'ed irreps)
    irvec::Vector{Int64}      # Vector that references irlabs of a parent BandRepSet; 
                              # nonzero entries correspond to an element in the band representation
    irreptags::Vector{String} # Vestigial, but quite handy for display'ing; this otherwise 
                              # requires recursive data sharing between BandRep and BandRepSet
end
wyck(BR::BandRep)    = BR.wyckpos
sitesym(BR::BandRep) = BR.sitesym
label(BR::BandRep)   = BR.label
humanreadable(BR::BandRep) = BR.irreptags
vec(BR::BandRep)     = BR.irvec

"""
    dim(BR::BandRep) --> Int64

Get the number of bands included in a single BandRep `BR`; i.e. the "band filling"
ν discussed in Po's papers.
"""
dim(BR::BandRep)     = BR.dim

# define the AbstractArray interface for BandRep
size(BR::BandRep)    = (length(vec(BR)),) # number of irreps samplable by BandRep
getindex(BR::BandRep, keys...) = vec(BR)[keys...]
firstindex(::BandRep) = 1
lastindex(BR::BandRep) = length(vec(BR))
IndexStyle(::BandRep) = IndexLinear()
eltype(::BandRep) = Int64

struct BandRepSet <: AbstractVector{BandRep}
    sgnum::Integer          # space group number, sequential
    bandreps::Vector{BandRep}
    kvs::Vector{KVec}       # Vector of 𝐤-points
    klabs::Vector{String}   # Vector of associated 𝐤-labels (in CDML notation)
    irlabs::Vector{String}  # Vector of (sorted) CDML irrep labels at _all_ 𝐤-points
    allpaths::Bool          # Whether all paths (true) or only maximal 𝐤-points (false) are included
    spinful::Bool           # Whether the band rep set includes (true) or excludes (false) spinful irreps
    timeinvar::Bool         # Whether the band rep set assumes time-reversal symmetry (true) or not (false) 
end
num(BRS::BandRepSet)         = BRS.sgnum
klabels(BRS::BandRepSet)     = BRS.klabs
kvecs(BRS::BandRepSet)       = BRS.kvs
hasnonmax(BRS::BandRepSet)   = BRS.allpaths
irreplabels(BRS::BandRepSet) = BRS.irlabs
isspinful(BRS::BandRepSet)   = BRS.spinful
istimeinvar(BRS::BandRepSet) = BRS.timeinvar
reps(BRS::BandRepSet)        = BRS.bandreps

# define the AbstractArray interface for BandRepSet
size(BRS::BandRepSet) = (length(reps(BRS)),) # number of distinct band representations
getindex(BRS::BandRepSet, keys...) = reps(BRS)[keys...]
firstindex(::BandRepSet) = 1
lastindex(BRS::BandRepSet) = length(reps(BRS))
IndexStyle(::BandRepSet) = IndexLinear()
eltype(::BandRepSet) = BandRep

"""
    matrix(BRS::BandRepSet[, includedim::Bool=false])

Return a matrix representation of `BRS::BandRepSet`, with band representations as columns 
and irreps over rows.

For `includedim=true` the band filling (i.e. `dim.(BRS)`) is included as the last row.
"""
function matrix(BRS::BandRepSet, includedim::Bool=false)
    Nⁱʳʳ, Nᵉᵇʳ = length(BRS[1]), length(BRS)
    M = Matrix{Int64}(undef, Nⁱʳʳ+includedim, Nᵉᵇʳ)
    @inbounds for (j, BR) in enumerate(BRS)
        for (i, v) in enumerate(vec(BR)) # bit over-explicit, but faster this way than with 
            M[i,j] = v                   # broadcasting/iterator interface (why!?)
        end
        if includedim
            M[Nⁱʳʳ+1,j] = dim(BR)
        end
    end

    return M
end 