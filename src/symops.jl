""" 
    read_sgops_xyzt(sgnum::Integer, dim::Integer=3)

Obtains the symmetry operations in xyzt format for a given space group
number `sgnum` by reading from json files; see `get_sgops` for additional
details. Much faster than crawling; generally preferred.
"""
function read_sgops_xyzt(sgnum::Integer, dim::Integer=3)
    if dim ∉ (1,2,3); throw(DomainError(dim, "dim must be 1, 2, or 3")); end
    if sgnum < 1 || dim == 3 && sgnum > 230 || dim == 2 && sgnum > 17 || dim == 1 && sgnum > 2
        throw(DomainError(sgnum, "sgnum must be in range 1:2 in 1D, 1:17 in 2D, and in 1:230 in 3D")) 
    end

    filepath = (@__DIR__)*"/../data/symops/"*string(dim)*"d/"*string(sgnum)*".json"
    sgops_str::Vector{String} = open(filepath) do io
        JSON2.read(io)
    end

    return sgops_str
end

""" 
    get_sgops(sgnum::Integer, dim::Integer=3) --> SpaceGroup

Obtains the space group symmetry operations in xyzt and matrix format
for a given space group number (`= sgnum`). The symmetry operations  
are specified relative to the conventional basis vector choices, i.e.
not necessarily primitive. 
If desired, operations on a primitive unit cell can be subsequently 
generated using `primitivize(...)` and `reduce_ops(...)`.

The default choices for basis vectors are specified in Bilbao as:
- Unique axis b (cell choice 1) for space groups within the
    monoclinic system.
- Obverse triple hexagonal unit cell for R space groups.
- Origin choice 2 - inversion center at (0,0,0) - for the
    centrosymmetric space groups for which there are two origin
    choices, within the orthorhombic, tetragonal and cubic systems.
"""
function get_sgops(sgnum::Integer, dim::Integer=3)
    sgops_str = read_sgops_xyzt(sgnum, dim)
    sgops = SymOperation.(sgops_str)

    return SpaceGroup(sgnum, sgops, dim)
end

function xyzt2matrix(s::String)
    ssub = split(s, ",")
    dim = length(ssub)
    xyzt2matrix!(zeros(Float64, dim, dim+1), ssub)
end

function xyzt2matrix!(O::Matrix{Float64}, s::Union{T, Array{T}} where T<:AbstractString)
    if s isa AbstractString
        itr = split(s,",")
    elseif s isa Array
        itr = s
    end

    for (i,op) in enumerate(itr)
        # rotation/inversion/reflection part
        nextidx = 1
        while true
            idx = findnext(r"x|y|z", op, nextidx);
            if !isnothing(idx)
                opchar = op[idx]
                if     opchar == "x"; j = 1; 
                elseif opchar == "y"; j = 2;
                elseif opchar == "z"; j = 3; end
                
                if idx[1] == 1 || op[prevind(op, idx[1])] == '+'
                    O[i,j] = 1.0
                elseif op[prevind(op, idx[1])] == '-'
                    O[i,j] = -1.0
                end
                nextidx = nextind(op, idx[end])
            else
                break
            end
        end
        
        # nonsymmorphic part/fractional translation part
        nonsymmorph = op[nextidx:end]
        if !isempty(nonsymmorph)
            slashidx = findfirst(x->x=='/',nonsymmorph)
            num=nonsymmorph[1:prevind(nonsymmorph, slashidx)]
            den=nonsymmorph[nextind(nonsymmorph, slashidx):end]
            O[i,end] = parse(Int64, num)/parse(Int64, den)
        end
    end
        
    return O
end

signaschar(x::Number) = signbit(x) ? '-' : '+'
const IDX2XYZ = ('x', 'y', 'z')

function matrix2xyzt(O::Matrix{T}) where T<:Real
    dim = size(O,1)
    buf = IOBuffer()
    # rotation/inversion/reflection part
    for (i, row) in enumerate(eachrow(O))
        # rotation/inversion/reflection part
        firstchar = true
        for j = 1:dim
            if !iszero(row[j])
                if !firstchar || signbit(row[j])
                    write(buf, signaschar(row[j]))
                end
                write(buf, IDX2XYZ[j]) 
                firstchar = false
            end
        end

        # nonsymmorphic/fractional translation part
        if size(O,2) == dim+1 # for size(O) = dim×dim+1, interpret as a space-group operation and check for nonsymmorphic parts; otherwise, assume a point-group operation
            if !iszero(row[end])
                write(buf, signaschar(row[end]))
                t = rationalize(float(row[end]), tol=1e-2) # convert to "minimal" Rational fraction (within nearest 1e-2 neighborhood)
                write(buf, string(abs(numerator(t)), '/', denominator(t)))
            end
        end
        if i != dim; write(buf, ','); end
    end

    return String(take!(buf))
end


"""
    issymmorph(op::SymOperation, cntr::Char) --> Bool

Checks whether a given symmetry operation `op` is symmorphic (true) or
nonsymmorphic (false). The operation is assumed to be given in a 
conventional basis; but the check requires that the translation is zero 
in a primitive basis. Accordingly, the centering `cntr` must provided.
"""
@inline function issymmorph(op::SymOperation, cntr::Char)
    P = primitivebasismatrix(cntr, dim(op))
    w_primitive = transform_translation(op, P, nothing) # translation in a primitive basis
    return iszero(w_primitive)
end
"""
    issymmorph(sg::AbstractGroup) --> Bool

Checks whether a given space group `sg` is symmorphic (true) or
nonsymmorphic (false).
"""
issymmorph(g::AbstractGroup) = all(op->issymmorph(op, centering(num(g), dim(g))), operations(g))

"""
    issymmorph(sgnum::Integer, dim::Integer=3) --> Bool

Checks whether a given space group `sgnum` is symmorphic (true) or
nonsymmorphic (false).
"""
issymmorph(sgnum::Integer, dim::Integer=3) = issymmorph(get_sgops(sgnum, dim))

# ----- POINT GROUP ASSOCIATED WITH SPACE/PLANE GROUP (FULL OR LITTLE) ---
"""
    pointgroup(ops:AbstractVector{SymOperation})

Computes the point group associated with a space group SG (characterized by
a set of operators `ops`, which, jointly with lattice translations generate 
the space group), obtained by "taking away" any translational parts and 
then reducing to the resulting unique rotational operations.
(technically, in the language of Bradley & Cracknell, this is the so-called
isogonal point group of SG; see Sec. 1.5).
"""
function pointgroup(ops::AbstractVector{SymOperation})
    # find SymOperations that are unique with respect to their rotational parts
    unique_rotation_ops = unique(rotation, ops) 
    # return rotation-only SymOperations from the above unique set
    return SymOperation.(hcat.(rotation.(unique_rotation_ops), Ref(zeros(Float64, dim(first(ops))))))
end
pointgroup(sg::AbstractGroup) = pointgroup(operations(sg))
pointgroup(pg::PointGroup) = operations(pg)
pointgroup(sgnum::Integer, dim::Integer=3) = pointgroup(get_sgops(sgnum, dim))

# ----- GROUP ELEMENT COMPOSITION -----
""" 
    (∘)(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation

Compose two symmetry operations `op1`={W₁|w₁} and `op2`={W₂|w₂}
using the composition rule (in Seitz notation)

    {W₁|w₁}{W₂|w₂} = {W₁*W₂|w₁+W₁*w₂}

for symmetry operations opᵢ = {Wᵢ|wᵢ}. By default, the translation part of
the {W₁*W₂|w₁+W₁*w₂} is reduced to the range [0,1], i.e. computed modulo 1.
This can be toggled off (or on) by the Boolean flag `modτ` (enabled, i.e. 
`true` by default). Returns another `SymOperation`.
"""
(∘)(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation = SymOperation((∘)(matrix(op1), matrix(op2), modτ))
function (∘)(op1::T, op2::T, modτ::Bool=true) where T<:Matrix{Float64}
    W′ = rotation(op1)*rotation(op2)
    w′ = translation(op1) .+ rotation(op1)*translation(op2)
    if modτ; w′ .= mod.(w′, 1.0); end

    return [W′ w′]
end
const compose = ∘



"""
    (⊚)(op1::T, op2::T) where T<:SymOperation -->  Vector{Float64}

Compose two symmetry operations `op1`={W₁|w₁} and `op2`={W₂|w₂} and
return the quotient of w₁+W₁*w₂ and 1. This functionality complements
`op1∘op2`, which yields the translation modulo 1; accordingly, 
`translation(op1∘op2) + op1⊚op2` yields the translation component
of the composition `op1` and `op2` **without** taking it modulo 1,
i.e. including any "trivial" lattice translation.

Note that ⊚ can be auto-completed in Julia via \\circledcirc+[tab]
""" 
function (⊚)(op1::T, op2::T) where T<:SymOperation
    # Translation result _without_ taking `mod`
    w′ = translation(op1) .+ rotation(op1)*translation(op2)  
    # Below, we combine `mod` and `rem` to ensure correctness in 
    # case any component `τ[i] < 0` (since `mod`, as used in ∘, 
    # is not the "partner" of `div`; `rem` is, in the sense 
    # `div(x,1) + rem(x,1) = x`, while `div(x,1) + mod(x,1) = x`
    # is only true for x ≥ 0).
    w′_lattice = div.(w′, 1.0) + rem.(w′, 1.0) .- mod.(w′, 1.0) 

    return w′_lattice
end

"""
    inv(op::SymOperation) --> SymOperation

Compute the inverse {W|w}⁻¹ of an operator `op`≡{W|w}.
"""
function inv(op::SymOperation)
    W = rotation(op)
    w = translation(op)

    W⁻¹ = inv(W)
    w⁻¹ = -W⁻¹*w

    return SymOperation([W⁻¹ w⁻¹])
end


"""
    multtable(ops::T) where T<:Union{Vector{SymOperation}, SpaceGroup}

Computes the multiplication table of a set of symmetry operations.
A MultTable is returned, which contains symmetry operations 
resulting from composition of `row ∘ col` operators; the table of 
indices give the symmetry operators relative to the ordering of 
`ops`.
"""
function multtable(ops::AbstractVector{SymOperation}; verbose::Bool=false)
    havewarned = false
    N = length(ops)
    indices = Matrix{Int64}(undef, N,N)
    for (row,oprow) in enumerate(ops)
        for (col,opcol) in enumerate(ops)
            op′ = matrix(oprow) ∘ matrix(opcol)
            match = findfirst(op′′ -> op′≈matrix(op′′), ops)
            if isnothing(match)
                if !havewarned
                    if verbose; @warn "The given operations do not form a group!"; end
                    havewarned = true
                end
                match = 0
            end
            @inbounds indices[row,col] = match
        end
    end
    return MultTable(ops, indices, !havewarned)
end
multtable(g::AbstractGroup) = multtable(operations(g))


checkmulttable(lgir::LGIrrep, αβγ=nothing; verbose::Bool=false) = begin
    ops = operations(lgir)
    sgnum = num(lgir); cntr = centering(sgnum, dim(first(ops)))
    primitive_ops = primitivize.(ops, cntr) # must do multiplication table in primitive basis, cf. choices for composition/∘
    checkmulttable(multtable(primitive_ops), lgir, αβγ; verbose=verbose)
end
function checkmulttable(mt::MultTable, lgir::LGIrrep, αβγ=nothing; verbose::Bool=false)
    havewarned = false
    irs = irreps(lgir, αβγ)
    ops = operations(lgir)
    k = kvec(lgir)(αβγ)
    N = length(ops)
    mtindices = indices(mt)
    checked = trues(N, N)
    for (row,irrow) in enumerate(irs)
        for (col,ircol) in enumerate(irs)
            @inbounds mtidx = mtindices[row,col]
            if iszero(mtidx) && !havewarned
                @warn "Provided multtable is not a group; cannot compare with irreps"
                checked[row,col] = false
                havewarned = true
            end
            ir′ = irrow*ircol
            # If 𝐤 is on the BZ boundary and if the little group is nonsymmorphic
            # the representation could be a ray representation (see Inui, p. 89),
            # such that DᵢDⱼ = αᵢⱼᵏDₖ with a phase factor αᵢⱼᵏ = exp(i*𝐤⋅𝐭₀) where
            # 𝐭₀ is a lattice vector 𝐭₀ = τᵢ + βᵢτⱼ - τₖ, for symmetry operations
            # {βᵢ|τᵢ}. To ensure we capture this, we include this phase here.
            # See Inui et al. Eq. (5.29) for explanation.
            # Note that the phase's sign is opposite to that used in many other 
            # conventions (e.g. Bradley & Cracknell, 1972, Eq. 3.7.7 & 3.7.8), 
            # but consistent with that used in Stokes' paper (see irreps(::LGIrrep)).
            # It is still a puzzle to me why I cannot successfully flip the sign 
            # of `ϕ` here and in `irreps(::LGIrrep)`.
            t₀ = translation(ops[row]) .+ rotation(ops[row])*translation(ops[col]) .- translation(ops[mtidx])
            ϕ =  2π*dot(k, t₀) # accumulated ray-phase
            match = ir′ ≈ cis(ϕ)*irs[mtidx] # cis(x) = exp(ix)
            if !match
                checked[row,col] = false
                if !havewarned
                    if verbose
                        println("""Provided irreps do not match group multiplication table for sg $(num(lgir)) in irrep $(label(lgir)):
                                 First failure at (row,col) = ($(row),$(col));
                                 Expected idx $(mtidx), got idx $(findall(ir′′ -> ir′′≈ir′, irs))
                                 Expected irrep = $(cis(ϕ)*irs[mtidx])
                                 Got irrep      = $(ir′)""")
                    end
                    havewarned = true
                end
            end
        end
    end
    return checked
end


# ----- LITTLE GROUP OF 𝐤 -----
# A symmetry operation g acts on a wave vector as (𝐤′)ᵀ = 𝐤ᵀg⁻¹ since we 
# generically operate with g on functions f(𝐫) via gf(𝐫) = f(g⁻¹𝐫), such that 
# the operation on a plane wave creates exp(i𝐤⋅g⁻¹𝐫); invariant plane waves 
# then define the little group elements {g}ₖ associated with wave vector 𝐤. 
# The plane waves are evidently invariant if 𝐤ᵀg⁻¹ = 𝐤ᵀ, or since g⁻¹ = gᵀ 
# (orthogonal transformations), if (𝐤ᵀg⁻¹)ᵀ = 𝐤 = (g⁻¹)ᵀ𝐤 = g𝐤; corresponding
# to the requirement that 𝐤 = g𝐤). Because we have g and 𝐤 in different bases
# (in the direct {𝐑} and reciprocal {𝐆} bases, respectively), we have to take 
# a little extra care here. Consider each side of the equation 𝐤ᵀ = 𝐤ᵀg⁻¹, 
# originally written in Cartesian coordinates, and rewrite each Cartesian term
# through basis-transformation to a representation we know (w/ P(𝐗) denoting 
# a matrix with columns of 𝐗m that facilitates this transformation):
#   𝐤ᵀ = [P(𝐆)𝐤(𝐆)]ᵀ = 𝐤(𝐆)ᵀP(𝐆)ᵀ                    (1)
#   𝐤ᵀg⁻¹ = [P(𝐆)𝐤(𝐆)]ᵀ[P(𝐑)g(𝐑)P(𝐑)⁻¹]⁻¹
#         = 𝐤(𝐆)ᵀP(𝐆)ᵀ[P(𝐑)⁻¹]⁻¹g(𝐑)⁻¹P(𝐑)⁻¹
#         = 𝐤(𝐆)ᵀ2πg(𝐑)⁻¹P(𝐑)⁻¹                       (2)
# (1+2): 𝐤′(𝐆)ᵀP(𝐆)ᵀ = 𝐤(𝐆)ᵀ2πg(𝐑)⁻¹P(𝐑)⁻¹
#     ⇔ 𝐤′(𝐆)ᵀ = 𝐤(𝐆)ᵀ2πg(𝐑)⁻¹P(𝐑)⁻¹[P(𝐆)ᵀ]⁻¹ 
#               = 𝐤(𝐆)ᵀ2πg(𝐑)⁻¹P(𝐑)⁻¹[2πP(𝐑)⁻¹]⁻¹
#               = 𝐤(𝐆)ᵀg(𝐑)⁻¹
#     ⇔  𝐤′(𝐆) = [g(𝐑)⁻¹]ᵀ𝐤(𝐆) = [g(𝐑)ᵀ]⁻¹𝐤(𝐆) 
# where we have used that P(𝐆)ᵀ = 2πP(𝐑)⁻¹ several times. Importantly, this
# essentially shows that we can consider g(𝐆) and g(𝐑) mutually interchangeable
# in practice.
# By similar means, one can show that 
#   [g(𝐑)⁻¹]ᵀ = P(𝐑)ᵀP(𝐑)g(𝐑)[P(𝐑)ᵀP(𝐑)]⁻¹
#             = [P(𝐆)ᵀP(𝐆)]⁻¹g(𝐑)[P(𝐆)ᵀP(𝐆)],
# by using that g(C)ᵀ = g(C)⁻¹ is an orthogonal matrix in the Cartesian basis.
# [ *) We transform from a Cartesian basis to an arbitrary 𝐗ⱼ basis via a 
# [    transformation matrix P(𝐗) = [𝐗₁ 𝐗₂ 𝐗₃] with columns of 𝐗ⱼ; a vector 
# [    v(𝐗) in the 𝐗-representation corresponds to a Cartesian vector v(C)≡v via
# [      v(C) = P(𝐗)v(𝐗)
# [    while an operator O(𝐗) corresponds to a Cartesian operator O(C)≡O via
# [      O(C) = P(𝐗)O(𝐗)P(𝐗)⁻¹
function littlegroup(ops::AbstractVector{SymOperation}, kv::KVec, cntr::Char='P')
    k₀, kabc = parts(kv)
    checkabc = !iszero(kabc)
    idxlist = [1]
    dim = length(k₀)
    for (idx, op) in enumerate(@view ops[2:end]) # note: `idx` is offset by 1 relative to position of op in ops
        k₀′, kabc′ = parts(compose(op, kv, checkabc)) # this is k₀(𝐆)′ = [g(𝐑)ᵀ]⁻¹k₀(𝐆)  
        diff = k₀′ .- k₀
        diff = primitivebasismatrix(cntr, dim)'*diff 
        kbool = all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), diff) # check if k₀ and k₀′ differ by a _primitive_ reciprocal vector
        abcbool = checkabc ? isapprox(kabc′, kabc, atol=DEFAULT_ATOL) : true # check if kabc == kabc′; no need to check for difference by a reciprocal vec, since kabc is in interior of BZ

        if kbool && abcbool # ⇒ part of little group
            push!(idxlist, idx+1) # `idx+1` is due to previously noted `idx` offset 
        end
    end
    return idxlist, view(ops, idxlist)
end
function littlegroup(sg::SpaceGroup, kv::KVec) 
    _, lgops = littlegroup(operations(sg), kv, centering(num(sg), dim(sg)))
    return LittleGroup{dim(sg)}(num(sg), kv, lgops)
end

function kstar(ops::Vector{SymOperation}, kv::KVec, cntr::Char)
    # we refer to kv by its parts (k₀, kabc) in the comments below
    kstar = [kv] 
    checkabc = !iszero(kv.kabc)
    d = dim(kv)
    for op in (@view ops[2:end])
        k₀′, kabc′ = parts(compose(op, kv, checkabc))

        newkbool = true
        for kv′′ in kstar
            k₀′′, kabc′′ = parts(kv′′)
            diff = k₀′ .- k₀′′
            diff = primitivebasismatrix(cntr, d)'*diff
            kbool = all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), diff)    # check if k₀ and k₀′ differ by a _primitive_ G-vector
            abcbool = checkabc ? isapprox(kabc′, kabc′′, atol=DEFAULT_ATOL) : true # check if kabc == kabc′ (no need to check for difference by G-vectors, since kabc ∈ interior of BZ)

            if kbool && abcbool # ⇒ we've already seen this KVec for (mod 𝐆) - we can skip it and go to next operator
                newkbool = false
                break # no need to check the rest of the kvecs currently in kstar; already found a match
            end
        end

        if newkbool
            push!(kstar, KVec(k₀′, kabc′))
        end
    end
    return kstar
end
kstar(sg::SpaceGroup, kv::KVec) = kstar(operations(sg), kv, centering(num(sg), dim(sg)))

"""
    (∘)(op::SymOperation, kv::KVec, checkabc::Bool=true) --> KVec

Computes the action of the SymOperation `op`=g on a KVec `kv`=k
using that g acts on k-vectors as k(G)′ = [g(R)ᵀ]⁻¹k(G), with g 
in an R-basis and k in a G-basis. Returns a new KVec, that is 
possibly distinct from its original only by a reciprocal lattice
vector (i.e. multiple of integers).

If `checkabc` = false, the free part of KVec is not transformed
(can be useful in situation where `kabc` is zero, and several 
transformations are requested).
"""
@inline function (∘)(op::SymOperation, kv::KVec, checkabc::Bool=true)
    k₀, kabc = parts(kv)
    k₀′ = rotation(op)'\k₀      
    kabc′ = checkabc ? rotation(op)'\kabc : kabc
    return KVec(k₀′, kabc′)
end



"""
    primitivize(op::SymOperation, cntr::Char) --> SymOperation

Transforms a symmetry operation `op`={W|w} from a conventional cell 
to a primitive cell (specified by its centering character `cntr`), 
then denoted {W′|w′}; i.e. performs a basis change 
    {W′|w′} = {P|p}⁻¹{W|w}{P|p}
where P and p describe basis change and origin shifts, respectively,
associated with the coordinate transformation. 

For additional details, see ITA6 Sec. 1.5.2.3, p. 84.
"""
function primitivize(op::SymOperation, cntr::Char)
    if cntr === 'P' || cntr === 'p' # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, dim(op))
        return transform(op, P, nothing)
    end
end

function conventionalize(op::SymOperation, cntr::Char)
    if cntr === 'P' || cntr === 'p' # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, dim(op))
        return transform(op, inv(P), nothing)
    end
end

""" 
    transform(op::SymOperation, P::Matrix{<:Real}, 
              p::Union{Vector{<:Real}, Nothing}=nothing) --> SymOperation

Transforms a symmetry operation `op = {W|w}` by a rotation matrix `P` and 
a translation vector `p` (can be `nothing` for zero-translations), producing
a new symmetry operation `op′ = {W′|w′}`: (see ITA6, Sec. 1.5.2.3.)
    {W′|w′} = {P|p}⁻¹{W|w}{P|p}
    with   W′ =  P⁻¹WP
           w′ = P⁻¹(w+Wp-p)
with the translation `w′` reduced to the range [0, 1). 

See also `primivitze` and `conventionalize`.
"""
# translation (usually zero; can then be given as `nothing`)
function transform(op::SymOperation, P::Matrix{<:Real}, 
                   p::Union{Vector{<:Real}, Nothing}=nothing)    
    W′ = transform_rotation(op, P)       # = P⁻¹WP       (+ rounding)
    w′ = transform_translation(op, P, p) # = P⁻¹(w+Wp-p)
                                         # with W ≡ rotation(op) and w ≡ translation(op)

    return SymOperation([W′ w′])
end

function transform_rotation(op::SymOperation, P::Matrix{<:Real})
    W = rotation(op)
    W′ = P\(W*P)        # = P⁻¹WP
    # clean up rounding-errors introduced by transformation (e.g. 
    # occassionally produces -0.0). The rotational part should 
    # always have integer coefficients in a valid lattice basis.
    @inbounds for (idx, el) in enumerate(W′) 
        rel = round(el)
        if !isapprox(el, rel, atol=DEFAULT_ATOL)
            throw(ErrorException("The transformed operator must have integer coefficients in its rotational part; got $(W′)"))
        end
        # since round(x) takes positive values x∈[0,0.5] to 0.0 and negative
        # values x∈[-0.5,-0.0] to -0.0 -- and since it is bad for us to have
        # both 0.0 and -0.0 -- we convert -0.0 to 0.0 here
        if rel===-zero(Float64); rel = zero(Float64); end

        W′[idx] = rel
    end
    return W′
end

function transform_translation(op::SymOperation, P::Matrix{<:Real}, 
                               p::Union{Vector{<:Real}, Nothing}=nothing)
    w = translation(op)

    if !isnothing(p)
        w′ = P\(w+rotation(op)*p-p)  # = P⁻¹(w+Wp-p)
    else
        w′ = P\w                     # = P⁻¹w  [with p = zero(dim(op))]
    end
    w′ .= mod.(w′, 1.0)
    return w′
end

function reduce_ops(ops::AbstractVector{SymOperation}, cntr::Char, conv_or_prim::Bool=true)
    P = primitivebasismatrix(cntr, dim(first(ops)))
    ops′ = transform.(ops, Ref(P), nothing)         # equiv. to `primitivize.(ops, cntr)` [but avoids loading P anew for each SymOperation]
    # remove equivalent operations
    ops′_reduced = SymOperation.(uniquetol(matrix.(ops′), atol=SGOps.DEFAULT_ATOL))

    if conv_or_prim # (true) return in conventional basis
        return transform.(ops′_reduced, Ref(inv(P))) # equiv. to conventionalize.(ops′_reduced, cntr)
    else            # (false) return in primitive basis
        return ops′_reduced
    end
end
reduce_ops(sg::SpaceGroup, conv_or_prim::Bool=true) = reduce_ops(operations(sg), centering(num(sg), dim(sg)), conv_or_prim)
reduce_ops(sgnum::Int64, dim::Int64=3, conv_or_prim::Bool=true) = reduce_ops(get_sgops(sgnum, dim), conv_or_prim)


"""
    findequiv(op::SymOperation, ops::AbstractVector{SymOperation}, cntr::Char) 
                                                --> Tuple{Int, Vector{Float64}}

Search for an operator `op′` in `ops` which is equivalent, modulo differences
by **primitive** lattice translations `Δw`, to `op`. Return the index of `op′` in 
`ops`, as well as the primitive translation difference `Δw`. If no match is found
returns `(nothing, nothing)`.

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


