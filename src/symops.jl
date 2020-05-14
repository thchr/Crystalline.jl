""" 
    read_sgops_xyzt(sgnum::Integer, dim::Integer=3)

Obtains the symmetry operations in xyzt format for a given space group
number `sgnum` by reading from json files; see `spacegroup` for additional
details. Much faster than crawling; generally preferred.
"""
function read_sgops_xyzt(sgnum::Integer, D::Integer=3)
    D ∉ (1,2,3) && _throw_invaliddim(D)
    if sgnum < 1 || D == 3 && sgnum > 230 || D == 2 && sgnum > 17 || D == 1 && sgnum > 2
        throw(DomainError(sgnum, "sgnum must be in range 1:2 in 1D, 1:17 in 2D, and in 1:230 in 3D")) 
    end

    filepath = (@__DIR__)*"/../data/sgops/"*string(D)*"d/"*string(sgnum)*".json"
    sgops_str::Vector{String} = open(filepath) do io
        JSON2.read(io)
    end

    return sgops_str
end

""" 
    spacegroup(sgnum::Integer, D::Integer=3) --> SpaceGroup{D}

Obtains the space group symmetry operations in xyzt and matrix format
for a given space group number (`= sgnum`) and dimensionality `D`.
The symmetry operations are specified relative to the conventional basis
vector choices, i.e. not necessarily primitive. 
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
@inline function spacegroup(sgnum::Integer, ::Val{D}=Val(3)) where D
    sgops_str = read_sgops_xyzt(sgnum, D)
    sgops = SymOperation{D}.(sgops_str)

    return SpaceGroup{D}(sgnum, sgops)
end
@inline spacegroup(sgnum::Integer, D::Integer) = spacegroup(sgnum, Val(D)) # behind a function barrier for type-inference's sake

function xyzt2matrix(s::String)
    ssub = split(s, ',')
    D = length(ssub)
    xyzt2matrix!(zeros(Float64, D, D+1), ssub)
end

function xyzt2matrix!(O::Matrix{Float64}, s::Union{T, AbstractVector{T}} where T<:AbstractString)
    if s isa AbstractString
        itr = split(s, ',')
    elseif s isa Array
        itr = s
    end

    @inbounds for (i, op) in enumerate(itr)
        # rotation/inversion/reflection part
        firstidx = nextidx = firstindex(op)
        while true
            idx = findnext(c -> c==='x' || c==='y' || c==='z', op, nextidx)
            if idx !== nothing
                opchar = op[idx]
                if      opchar === 'x';   j = 1; 
                elseif  opchar === 'y';   j = 2;
                else #= opchar === 'z' =# j = 3; end # opchar can only be 'z' at this point; no need to check
                
                previdx = prevind(op, idx)
                if idx == firstidx || op[previdx] === '+'
                    O[i,j] = 1.0
                elseif op[previdx] === '-'
                    O[i,j] = -1.0
                end
                nextidx = nextind(op, idx)
            else
                break
            end
        end
        
        # nonsymmorphic part/fractional translation part
        lastidx = lastindex(op)
        if nextidx ≤ lastidx # ... then there's stuff "remaining" in op; a nonsymmorphic part
            slashidx = findnext(==('/'), op, nextidx)
            if slashidx !== nothing # interpret as integer fraction
                num = SubString(op, nextidx, prevind(op, slashidx))
                den = SubString(op, nextind(op, slashidx), lastidx)
                O[i,end] = parse(Int64, num)/parse(Int64, den)
            else                    # interpret at floating point number
                O[i,end] = parse(Float64, SubString(op, nextidx, lastidx))
            end
        end
    end
        
    return O
end

signaschar(x::Number) = signbit(x) ? '-' : '+'
const IDX2XYZ = ('x', 'y', 'z')

function matrix2xyzt(O::AbstractMatrix{<:Real})
    D = size(O,1)
    buf = IOBuffer()
    @inbounds for i in Base.OneTo(D)
        # point group part
        firstchar = true
        for j in Base.OneTo(D)
            Oᵢⱼ = O[i,j]
            if !iszero(Oᵢⱼ)
                if !firstchar || signbit(Oᵢⱼ)
                    print(buf, signaschar(Oᵢⱼ))
                end
                print(buf, IDX2XYZ[j]) 
                firstchar = false
            end
        end

        # nonsymmorphic/fractional translation part
        if size(O,2) == D+1 # for size(O) = D×D+1, interpret as a space-group operation and 
                            # check for nonsymmorphic parts
            if !iszero(O[i,D+1])
                fractionify!(buf, O[i,D+1])
            end
        end
        i != D && print(buf, ',')
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
issymmorph(g::Union{SpaceGroup,LittleGroup}) = all(op->issymmorph(op, centering(g)), operations(g))

"""
    issymmorph(sgnum::Integer, D::Integer=3) --> Bool

Checks whether a given space group `sgnum` (of dimensionality `D`)
is symmorphic (true) or nonsymmorphic (false).
"""
issymmorph(sgnum::Integer, D::Integer=3) = issymmorph(spacegroup(sgnum, D))

# ----- POINT GROUP ASSOCIATED WITH SPACE/PLANE GROUP (FULL OR LITTLE) ---
"""
    pointgroup(ops:AbstractVector{SymOperation{D}})
    pointgroup(sg::AbstractGroup)

Computes the point group associated with a space group `sg` (characterized by
a set of operators `ops`, which, jointly with lattice translations generate 
the space group), obtained by "taking away" any translational parts and 
then reducing to the resulting unique rotational operations.
(technically, in the language of Bradley & Cracknell, this is the so-called
isogonal point group of `sg`; see Sec. 1.5).

Returns a `Vector` of `SymOperation`s.
"""
function pointgroup(ops::AbstractVector{SymOperation{D}}) where D
    # find SymOperations that are unique with respect to their rotational parts
    unique_rotation_ops = unique(rotation, ops) 
    # return rotation-only SymOperations from the above unique set
    return SymOperation{D}.(hcat.(rotation.(unique_rotation_ops), Ref(zeros(Float64, D))))
    # TODO: Return a PointGroup?
end
pointgroup(sg::Union{SpaceGroup,LittleGroup}) = pointgroup(operations(sg))

# ----- GROUP ELEMENT COMPOSITION -----
""" 
    (∘)(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation

Compose two symmetry operations `op1`={W₁|w₁} and `op2`={W₂|w₂}
using the composition rule (in Seitz notation)

    {W₁|w₁}{W₂|w₂} = {W₁*W₂|w₁+W₁*w₂}

for symmetry operations opᵢ = {Wᵢ|wᵢ}. By default, the translation part of
the {W₁*W₂|w₁+W₁*w₂} is reduced to the range [0,1], i.e. computed modulo 1.
This can be toggled off (or on) by the Boolean flag `modτ` (enabled, i.e. 
`true`, by default). Returns another `SymOperation`.
"""
function(∘)(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation
    T((∘)(matrix(op1), matrix(op2), modτ))
end
function (∘)(op1::T, op2::T, modτ::Bool=true) where T<:AbstractMatrix{Float64}
    W′ = rotation(op1)*rotation(op2)
    w′ = translation(op1) .+ rotation(op1)*translation(op2)

    if modτ
        reduce_translation_to_unitrange!(w′)
    end

    return [W′ w′]
end
const compose = ∘

function reduce_translation_to_unitrange!(w::AbstractVector{Float64}) # mutates w; reduces components to range [0.0, 1.0[
    # naïve approach to achieve semi-robust reduction of integer-translation
    # via a slightly awful "approximate" modulo approach; basically just the
    # equivalent of w′ .= mod.(w′,1.0), but reducing in a range DEFAULT_ATOL 
    # around each integer.
    w .= mod.(w, 1.0)
    # sometimes, mod(w′, 1.0) can omit reducing values that are very nearly 1.0
    # due to floating point errors: we use a tolerance here to round everything 
    # close to 0.0 or 1.0 exactly to 0.0
    @simd for i in eachindex(w)
        if isapprox(round(w[i]), w[i], atol=DEFAULT_ATOL)
            w[i] = zero(eltype(w))
        end
    end
    return w
end
reduce_translation_to_unitrange(w::AbstractVector{Float64}) = reduce_translation_to_unitrange!(copy(w)) # non-mutating variant

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
    # Then we take w′ modulo lattice vectors
    w′′ = reduce_translation_to_unitrange(w′)
    # Then we subtract the two (reuse w′′ to avoid additional allocations)
    w′′ .= w′ .- w′′
    return w′′
end

"""
    inv(op::SymOperation{D}) --> SymOperation{D}

Compute the inverse {W|w}⁻¹ ≡ {W⁻¹|-W⁻¹w} of an operator `op` ≡ {W|w}.
"""
function inv(op::T) where T<:SymOperation
    W = rotation(op)
    w = translation(op)

    W⁻¹ = inv(W)
    w⁻¹ = -W⁻¹*w

    return T([W⁻¹ w⁻¹])
end


"""
    multtable(ops::AbstractVector{<:SymOperation{D}})

Compute the multiplication (or Cayley) table of a set of symmetry operations.
A MultTable is returned, which contains symmetry operations 
resulting from composition of `row ∘ col` operators; the table of 
indices give the symmetry operators relative to the ordering of 
`ops`.
"""
function multtable(ops::AbstractVector{SymOperation{D}}; verbose::Bool=false) where D
    havewarned = false
    N = length(ops)
    indices = Matrix{Int64}(undef, N,N)
    for (row,oprow) in enumerate(ops)
        for (col,opcol) in enumerate(ops)
            op′ = oprow ∘ opcol
            match = findfirst(op′′ -> op′≈op′′, ops)
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
    return MultTable{D}(ops, indices, !havewarned)
end


function checkmulttable(lgir::LGIrrep{D}, αβγ=nothing; verbose::Bool=false) where D
    ops = operations(lgir)
    sgnum = num(lgir); cntr = centering(sgnum, D)
    primitive_ops = primitivize.(ops, cntr) # must do multiplication table in primitive basis, cf. choices for composition/∘
    checkmulttable(multtable(primitive_ops), lgir, αβγ; verbose=verbose)
end
function checkmulttable(mt::MultTable, ir::AbstractIrrep, αβγ=nothing; verbose::Bool=false)
    havewarned = false
    Ds = irreps(ir, αβγ)
    ops = operations(ir)
    if ir isa LGIrrep
        k = kvec(ir)(αβγ)
    end
    N = length(ops)
    mtindices = indices(mt)
    checked = trues(N, N)
    for (i,Dⁱ) in enumerate(Ds)     # rows
        for (j,Dʲ) in enumerate(Ds) # cols
            @inbounds mtidx = mtindices[i,j]
            if iszero(mtidx) && !havewarned
                @warn "Provided multtable is not a group; cannot compare with irreps"
                checked[i,j] = false
                havewarned = true
            end
            Dⁱʲ = Dⁱ*Dʲ
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
            if ir isa LGIrrep
                t₀ = translation(ops[i]) .+ rotation(ops[i])*translation(ops[j]) .- 
                     translation(ops[mtidx])
                ϕ =  2π*dot(k, t₀) # accumulated ray-phase
                match = Dⁱʲ ≈ cis(ϕ)*Ds[mtidx] # cis(x) = exp(ix)
            else
                match = Dⁱʲ ≈ Ds[mtidx]
            end
            if !match
                checked[i,j] = false
                if !havewarned
                    if verbose
                        println("""Provided irreps do not match group multiplication table for group $(num(ir)) in irrep $(label(ir)):
                                 First failure at (row,col) = ($(i),$(j));
                                 Expected idx $(mtidx), got idx $(findall(D′ -> D′≈Dⁱʲ, Ds))""")
                        print("Expected irrep = ")
                        if ir isa LGIrrep
                            println(cis(ϕ)*Ds[mtidx])
                        else
                            println(Dⁱʲ)
                        end
                        println("Got irrep      = $(Dⁱʲ)")
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
function littlegroup(ops::AbstractVector{SymOperation{D}}, kv::KVec, cntr::Char='P') where D
    k₀, kabc = parts(kv)
    checkabc = !iszero(kabc)
    idxlist = [1]
    for (idx, op) in enumerate(@view ops[2:end]) # note: `idx` is offset by 1 relative to position of op in ops
        k₀′, kabc′ = parts(compose(op, kv, checkabc)) # this is k₀(𝐆)′ = [g(𝐑)ᵀ]⁻¹k₀(𝐆)  
        diff = k₀′ .- k₀
        diff = primitivebasismatrix(cntr, D)'*diff 
        kbool = all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), diff) # check if k₀ and k₀′ differ by a _primitive_ reciprocal vector
        abcbool = checkabc ? isapprox(kabc′, kabc, atol=DEFAULT_ATOL) : true # check if kabc == kabc′; no need to check for difference by a reciprocal vec, since kabc is in interior of BZ

        if kbool && abcbool # ⇒ part of little group
            push!(idxlist, idx+1) # `idx+1` is due to previously noted `idx` offset 
        end
    end
    return idxlist, view(ops, idxlist)
end
function littlegroup(sg::SpaceGroup, kv::KVec) 
    _, lgops = littlegroup(operations(sg), kv, centering(sg))
    return LittleGroup{dim(sg)}(num(sg), kv, "", lgops)
end

function kstar(ops::AbstractVector{SymOperation{D}}, kv::KVec, cntr::Char) where D
    # we refer to kv by its parts (k₀, kabc) in the comments below
    kstar = [kv] 
    checkabc = !iszero(kv.kabc)
    for op in (@view ops[2:end])
        k₀′, kabc′ = parts(compose(op, kv, checkabc))

        newkbool = true
        for kv′′ in kstar
            k₀′′, kabc′′ = parts(kv′′)
            diff = k₀′ .- k₀′′
            diff = primitivebasismatrix(cntr, D)'*diff
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
kstar(sg::SpaceGroup, kv::KVec) = kstar(sg, kv, centering(sg))

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
function primitivize(op::SymOperation{D}, cntr::Char) where D
    if (D == 3 && cntr === 'P') || (D == 2 && cntr === 'p')
        # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, D)
        return transform(op, P, nothing)
    end
end

function conventionalize(op::SymOperation{D}, cntr::Char) where D
    if (D == 3 && cntr === 'P') || (D == 2 && cntr === 'p')
        # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, D)
        return transform(op, inv(P), nothing)
    end
end

function primitivize(kv::KVec, cntr::Char)
    k₀, kabc = parts(kv)
    P = primitivebasismatrix(cntr, dim(kv))
    # P transforms reciprocal coordinates as 𝐤′ = Pᵀ𝐤
    # while P transforms direct coordinates as r′=P⁻¹r,
    # see ITA7 Sec. 1.5.1.2 and 1.5.2.1 (and note the 
    # distinction between transforming the basis and
    # the coordinates of a vector!).
    k₀′ = P'*k₀
    kabc′ = P'*kabc
    return KVec(k₀′, kabc′)
end

""" 
    transform(op::SymOperation, P::Matrix{<:Real}, 
              p::Union{Vector{<:Real}, Nothing}=nothing,
              modw::Bool=true)                          --> SymOperation

Transforms a symmetry operation `op = {W|w}` by a rotation matrix `P` and 
a translation vector `p` (can be `nothing` for zero-translations), producing
a new symmetry operation `op′ = {W′|w′}`: (see ITA6, Sec. 1.5.2.3.)

        {W′|w′} = {P|p}⁻¹{W|w}{P|p}
        with   W′ = P⁻¹WP and w′ = P⁻¹(w+Wp-p)

By default, the translation part of `op′`, i.e. `w′`, is reduced to the range 
[0,1), i.e. computed modulo 1 (corresponding to `modw=true`). This can be 
disabled by setting `modw=false`.

See also `primitivize` and `conventionalize`. 
"""
function transform(op::SymOperation{D}, P::AbstractMatrix{<:Real}, 
                   p::Union{AbstractVector{<:Real}, Nothing}=nothing,
                   modw::Bool=true) where D
    W′ = transform_rotation(op, P)             # = P⁻¹WP       (+ rounding)
    w′ = transform_translation(op, P, p, modw) # = P⁻¹(w+Wp-p)
                                               # with W ≡ rotation(op) and w ≡ translation(op)

    return SymOperation{D}([W′ w′])
end

function transform_rotation(op::SymOperation, P::AbstractMatrix{<:Real})
    W = rotation(op)
    W′ = P\(W*P)        # = P⁻¹WP
    # clean up rounding-errors introduced by transformation (e.g. 
    # occassionally produces -0.0). The rotational part will 
    # always have integer coefficients if it is in the conventional
    # or primitive basis of its lattice; if transformed to a nonstandard
    # lattice, it might not have that though.
    @inbounds for (idx, el) in enumerate(W′)
        rel = round(el)
        if !isapprox(el, rel, atol=DEFAULT_ATOL)
            rel = el # non-standard lattice transformation; fractional elements 
                     # (this is why we need Float64 in SymOperation{D})
        end
        # since round(x) takes positive values x∈[0,0.5] to 0.0 and negative
        # values x∈[-0.5,-0.0] to -0.0 -- and since it is bad for us to have
        # both 0.0 and -0.0 -- we convert -0.0 to 0.0 here
        if rel===-zero(Float64); rel = zero(Float64); end

        W′[idx] = rel
    end
    return W′
end

function transform_translation(op::SymOperation, P::AbstractMatrix{<:Real}, 
                               p::Union{AbstractVector{<:Real}, Nothing}=nothing,
                               modw::Bool=true)
    w = translation(op)

    if !isnothing(p)
        w′ = P\(w+rotation(op)*p-p)  # = P⁻¹(w+Wp-p)
    else
        w′ = P\w                     # = P⁻¹w  [with p = zero(dim(op))]
    end
    
    if modw
        return mod.(w′, 1.0)
    else
        return w′
    end
end

# TODO: Maybe implement this in mutating form; lots of unnecessary allocations below in many usecases
function reduce_ops(ops::AbstractVector{SymOperation{D}}, cntr::Char, conv_or_prim::Bool=true) where D
    P = primitivebasismatrix(cntr, D)
    ops′ = transform.(ops, Ref(P))         # equiv. to `primitivize.(ops, cntr)` [but avoids loading P anew for each SymOperation]
    # remove equivalent operations
    ops′_reduced = SymOperation{D}.(uniquetol(matrix.(ops′), atol=Crystalline.DEFAULT_ATOL))

    if conv_or_prim # (true) return in conventional basis
        return transform.(ops′_reduced, Ref(inv(P))) # equiv. to conventionalize.(ops′_reduced, cntr)
    else            # (false) return in primitive basis
        return ops′_reduced
    end
end
@inline function reduce_ops(slg::Union{<:SpaceGroup, <:LittleGroup}, conv_or_prim::Bool=true)
    return reduce_ops(operations(slg), centering(slg), conv_or_prim)
end
primitivize(sg::T) where T<:SpaceGroup = T(num(sg), reduce_ops(sg, false))
function primitivize(lg::T) where T<:LittleGroup 
    cntr = centering(lg)
    kv′  = primitivize(kvec(lg), cntr)              # transform both k-point and operations
    ops′ = reduce_ops(operations(lg), cntr, false)
    return T(num(lg), kv′, klabel(lg), ops′)
end

"""
    cartesianize(op::SymOperation{D}, Rs::DirectBasis{D}) --> SymOperation{D}

Convert a `SymOperation` `opˡ` from the lattice basis to a Cartesian basis, by computing the
transformed operators `opᶜ = 𝐑*opˡ*𝐑⁻¹` via the Cartesian basis matrix 𝐑 (whose columns are
the `DirectBasis` vectors `Rs[i]`). 

# Note 1
The matrix 𝐑 maps vectors coefficients in a lattice basis 𝐯ˡ to coefficients in a Cartesian
basis 𝐯ᶜ as 𝐯ˡ = 𝐑⁻¹𝐯ᶜ and vice versa as 𝐯ᶜ = 𝐑𝐯ˡ. Since a general transformation P 
transforms an "original" vectors with coefficients 𝐯 to new coefficients 𝐯′ via 𝐯′ = P⁻¹𝐯
and since we here here consider the lattice basis as the "original" bais we have P = 𝐑⁻¹. 
As such, the transformation of the operator `op` transforms as `opᶜ = P⁻¹*opˡ*P`, i.e.
`opᶜ = transform(opˡ,P) = transform(opˡ,𝐑⁻¹)`.

# Note 2
The display (e.g. Seitz and xyzt notation) of `SymOperation`s e.g. in the REPL implicitly
assumes integer coefficients for its point-group matrix: as a consequence, displaying 
`SymOperation`s in a Cartesian basis may yield undefine behavior. The matrix representation
remains valid, however.
"""
function cartesianize(op::SymOperation{D}, Rs::DirectBasis{D}) where D
    𝐑 = basis2matrix(Rs)
    # avoids inv(𝐑) by not calling out to transform(opˡ, inv(𝐑))
    op′ = SymOperation{D}([𝐑*rotation(op)/𝐑 𝐑\translation(op)]) 
    return op′
end
cartesianize(sg::SpaceGroup{D}, Rs::DirectBasis{D}) where D = SpaceGroup{D}(num(sg), cartesianize.(operations(sg), Ref(Rs)))

"""
    findequiv(op::SymOperation, ops::AbstractVector{SymOperation{D}}, cntr::Char) 
                                                --> Tuple{Int, Vector{Float64}}

Search for an operator `op′` in `ops` which is equivalent, modulo differences
by **primitive** lattice translations `Δw`, to `op`. Return the index of `op′` in 
`ops`, as well as the primitive translation difference `Δw`. If no match is found
returns `(nothing, nothing)`.

The small irreps of `op` at wavevector k, Dⱼᵏ[`op`], can be computed from 
the small irreps of `op′`, Dⱼᵏ[`op′`], via Dⱼᵏ[`op`] = exp(2πik⋅`Δw`)Dⱼᵏ[`op′`]
"""
function findequiv(op::SymOperation{D}, ops::AbstractVector{SymOperation{D}}, cntr::Char) where D
    W = rotation(op)
    w = translation(op)

    P = primitivebasismatrix(cntr, D)
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
    _findsubgroup(opsᴳ, opsᴴ) --> Bool, Vector{Int64}

Determine whether the group H (with operators `opsᴴ`) is a subgroup
of the group G (with operators `opsᴳ`), i.e. whether H<G, and returns
an indexing vector `idxs` of `opsᴳ` into `opsᴴ` (empty if `false`), such
that `opsᴳ[idxs] ≡ H`. 
The first return argument is a Boolean (whether H<G); the second is `idxs`.

"""
function _findsubgroup(opsᴳ::T, opsᴴ::T) where T<:AbstractVector{<:SymOperation{<:Any}}
    idxsᴳ²ᴴ = Vector{Int64}(undef, length(opsᴴ))
    @inbounds for (idxᴴ, opᴴ) in enumerate(opsᴴ)
        idxᴳ = findfirst(==(opᴴ), opsᴳ)
        if idxᴳ !== nothing
            idxsᴳ²ᴴ[idxᴴ] = idxᴳ
        else
            return false, Int64[]
        end
    end
    return true, idxsᴳ²ᴴ
end
_findsubgroup(G::T, H::T) where T<:SpaceGroup = _findsubgroup(operations(G), operations(H))

"""
    issubgroup(opsᴳ::T, opsᴴ::T) --> Bool

Determine whether the operations in group H are a subgroup of the group G (each with 
operations `opsᴳ` and `opsᴴ`, respectively, of type `T::AbstractVector{SymOperation{D}}`),
i.e. whether H<G. Specifically, this requires that G and H are both groups and that 
for every h∈H there exists an element g∈G such that h=g.

Returns a Boolean answer (true if normal, false if not).

**Note:** This compares space groups rather than space group types, i.e. the 
comparison assumes a matching setting choice between H and G. To compare space 
group types with different conventional settings, they must first be transformed
to a shared setting.
"""
function issubgroup(opsᴳ::T, opsᴴ::T) where T<:AbstractVector{SymOperation{D}} where D
    ΔW = Matrix{Float64}(undef, D, D) # work matrices
    Δw = Vector{Float64}(undef, D)
    for h in opsᴴ
        found = false
        for g in opsᴳ
            ΔW .= rotation(h) .- rotation(g)
            Δw .= translation(h) .- translation(g)

            @inbounds @simd for i in Base.OneTo(D) # consider two operations identical if they differ by a near-integer translation
                rΔwᵢ = round(Δw[i])
                if isapprox(Δw[i], rΔwᵢ, atol=DEFAULT_ATOL)
                    Δw[i] = zero(Float64)
                end
            end
            
            if norm(ΔW) < DEFAULT_ATOL && norm(Δw) < DEFAULT_ATOL
                found = true
                continue
            end
        end
        if !found
            return false
        end
    end
    return true
end
issubgroup(G::T, H::T) where T<:SpaceGroup = issubgroup(operations(G), operations(H))


"""
    isnormal(opsᴳ::T, opsᴴ::T; verbose::Bool=false) --> Bool

Determine whether the operations in group H are normal in the group G (each with 
operations `opsᴳ` and `opsᴴ`, respectively, of type `T::AbstractVector{SymOperation{D}}`),
in the sense that 
    
    ghg⁻¹ ∈ H ∀ g∈G, h∈H

Returns a Boolean answer (true if normal, false if not).

**Note:** This that this compares space groups rather than space group types, i.e. the 
comparison assumes a matching setting choice between H and G. To compare space 
group types with different conventional settings, they must first be transformed
to a shared setting.
"""
function isnormal(opsᴳ::T, opsᴴ::T; verbose::Bool=false) where T<:AbstractVector{<:SymOperation{<:Any}}
    for g in opsᴳ
        g⁻¹ = inv(g)
        for h in opsᴴ
            # check if ghg⁻¹ ∉ G
            h′ = g∘h∘g⁻¹
            if !isapproxin(h′, opsᴴ, atol=Crystalline.DEFAULT_ATOL)
                if verbose
                    println("\nNormality-check failure:\n",
                            "Found h′ = ", seitz(h′), "\n",
                            "But h′ should be an element of the group: ", 
                            join(seitz.(opsᴴ), ", "))
                end
                return false
            end
        end
    end
    
    return true
end
isnormal(G::T, H::T) where T<:SpaceGroup = isnormal(operations(G), operations(H))