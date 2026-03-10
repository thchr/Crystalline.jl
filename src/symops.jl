"""
    S"..." --> SymOperation

Construct a `SymOperation` from a triplet form given as a string.

The expressions specify coordinate mappings from initial positions `x,y,z` (reduced 
coordinates relative to lattice basis) with constant terms as translations.

## Examples

```jldoctest
julia> S"-y,x"
4⁺ ──────────────────────────────── (-y,x)
 ┌ 0 -1 ╷ 0 ┐
 └ 1  0 ╵ 0 ┘

julia> S"x,-y,z+1/2"
{m₀₁₀|0,0,½} ──────────────── (x,-y,z+1/2)
 ┌ 1  0  0 ╷   0 ┐
 │ 0 -1  0 ┆   0 │
 └ 0  0  1 ╵ 1/2 ┘
```

See also [`SymOperation`](@ref).
"""
macro S_str(s)
    SymOperation(s)
end

function xyzt2matrix(s::AbstractString, Dᵛ::Val{D}) where D
    W, w = xyzt2components(s, Dᵛ)
    return hcat(SMatrix(W), SVector(w))
end

if !isdefined(Base, :eachsplit) # compat (https://github.com/JuliaLang/julia/pull/39245)
    eachsplit(str::T, splitter; limit::Integer=0, keepempty::Bool=true) where T = 
                                                    split(str, splitter; limit, keepempty)
end
function xyzt2components(s::AbstractString, ::Val{D}) where D
    # initialize zero'd MArrays for rotation/translation parts (allocation will be elided)
    W = zero(MMatrix{D, D, Float64}) # rotation
    w = zero(MVector{D, Float64})    # translation
    
    # "fill in" `W` and `w` according to content of `xyzts`
    xyzts = eachsplit(s, ',')
    xyzt2components!(W, w, xyzts)

    # convert to SArrays (elides allocation since `xyzt2components!` is inlined)
    return SMatrix(W), SVector(w)
end

@inline function xyzt2components!(W::MMatrix{D,D,T}, w::MVector{D,T},
                                  xyzts #= iterator of AbstractStrings =#) where {D,T<:Real}

    chars = D == 3 ? ('x','y','z') : D == 2 ? ('x','y') : ('x',)
    d = 0 # dimension-consistency check

    @inbounds for (i,s) in enumerate(xyzts)
        # rotation/inversion/reflection part
        firstidx = nextidx = firstindex(s)
        while (idx = findnext(c -> c ∈ chars, s, nextidx)) !== nothing
            idx = something(idx)
            c = s[idx]
            j = c=='x' ? 1 : (c=='y' ? 2 : 3)
            
            if idx == firstidx
                W[i,j] = one(T)
            else
                previdx = prevind(s, idx)
                while (c′=s[previdx]; isspace(s[previdx]))
                    previdx = prevind(s, previdx)
                    previdx ≤ firstidx && break
                end
                if c′ == '+' || isspace(c′)
                    W[i,j] = one(T)
                elseif c′ == '-'
                    W[i,j] = -one(T)
                else
                    throw(ArgumentError("failed to parse provided string representation"))
                end
            end
            nextidx = nextind(s, idx)
        end

        # nonsymmorphic part/fractional translation part
        lastidx = lastindex(s)
        if nextidx ≤ lastidx # ⇒ stuff "remaining" in `s`; a nonsymmorphic part
            slashidx = findnext(==('/'), s, nextidx)
            if slashidx !== nothing # interpret as integer fraction
                num = SubString(s, nextidx, prevind(s, slashidx))
                den = SubString(s, nextind(s, slashidx), lastidx)
                w[i] = convert(T, parse(Int, num)/parse(Int, den))
            else                    # interpret at number of type `T`
                w[i] = parse(T, SubString(s, nextidx, lastidx))
            end
        end
        d = i
    end
    d == D || throw(DimensionMismatch("incompatible matrix size and string format"))

    return W, w
end

const IDX2XYZ = ('x', 'y', 'z')
function matrix2xyzt(O::AbstractMatrix{<:Real})
    D = size(O,1)
    buf = IOBuffer()
    @inbounds for i in OneTo(D)
        # point group part
        firstchar = true
        for j in OneTo(D)
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

Return whether a given symmetry operation `op` is symmorphic (`true`) or nonsymmorphic
(`false`). 

The operation is assumed provided in conventional basis with centering type `cntr`: 
checking symmorphism is then equivalent to checking whether the operation's translation
part is zero or a lattice vector in the associated primitive basis.
"""
@inline function issymmorph(op::SymOperation{D}, cntr::Char) where D
    P = primitivebasismatrix(cntr, Val(D))
    w_primitive = transform_translation(op, P, nothing) # translation in a primitive basis
    return iszero(w_primitive)
end

"""
    issymmorph(sg::Union{SpaceGroup, LittleGroup}) --> Bool

Return whether a given space group `sg` is symmorphic (`true`) or nonsymmorphic (`false`).
"""
function issymmorph(g::AbstractGroup)
    all(op->issymmorph(op, centering(g)), operations(g))
end
issymmorph(::PointGroup) = true

"""
    issymmorph(sgnum::Integer, D::Integer=3) --> Bool

Return whether the space group with number `sgnum` and dimensionality `D` is symmorphic 
(`true`) or nonsymmorphic (`false`).

Equivalent to `issymmorph(spacegroup(sgnum, D))` but uses memoization for performance.
"""
function issymmorph(sgnum::Integer, D::Integer=3)
    @boundscheck _check_valid_sgnum_and_dim(sgnum, D)

    # check against memoized look-up of calling `issymmorph(spacegroup(sgnum, D))`
    if D == 3
        return sgnum ∈ (1, 2, 3, 5, 6, 8, 10, 12, 16, 21, 22, 23, 25, 35, 38, 42, 44, 47, 
                        65, 69, 71, 75, 79, 81, 82, 83, 87, 89, 97, 99, 107, 111, 115, 119, 
                        121, 123, 139, 143, 146, 147, 148, 149, 150, 155, 156, 157, 160, 
                        162, 164, 166, 168, 174, 175, 177, 183, 187, 189, 191, 195, 196, 
                        197, 200, 202, 204, 207, 209, 211, 215, 216, 217, 221, 225, 229)
    elseif D == 2
        return sgnum ∉ (4, 7, 8, 12)
    elseif D == 1
        return true
    else
        error("unreachable: file a bug report")
    end
end

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
    pgops = unique(rotation, ops) 
    # return rotations only from the above unique set (set translations to zero)
    pgops .= SymOperation.(rotation.(pgops))
    return pgops
end
pointgroup(sg::Union{SpaceGroup,LittleGroup}) = pointgroup(operations(sg))

# ----- GROUP ELEMENT COMPOSITION -----
""" 
    compose(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation

Compose two symmetry operations `op1` ``= \\{W₁|w₁\\}`` and `op2` ``= \\{W₂|w₂\\}``
using the composition rule (in Seitz notation)

``\\{W₁|w₁\\}\\{W₂|w₂\\} = \\{W₁W₂|w₁+W₁w₂\\}``

By default, the translation part of the ``\\{W₁W₂|w₁+W₁w₂\\}`` is reduced to the range
``[0,1[``, i.e. computed modulo 1. This can be toggled off (or on) by the Boolean flag
`modτ` (enabled, i.e. `true`, by default). Returns another `SymOperation`.

The multiplication operator `*` is overloaded for `SymOperation`s to call `compose`, in the
manner `op1 * op2 = compose(op1, op2, modτ=true)`.
"""
function compose(op1::T, op2::T, modτ::Bool=true) where T<:SymOperation
    T(compose(unpack(op1)..., unpack(op2)..., modτ)...)
end
function compose(W₁::T, w₁::R, W₂::T, w₂::R, modτ::Bool=true) where T<:SMatrix{D,D,<:Real} where R<:SVector{D,<:Real} where D
    W′ = W₁*W₂
    w′ = w₁ + W₁*w₂

    if modτ
        w′ = reduce_translation_to_unitrange(w′)
    end

    return W′, w′
end
(*)(op1::T, op2::T) where T<:SymOperation = compose(op1, op2)

"""
    reduce_translation_to_unitrange_q(w::AbstractVector{<:Real}) --> :: typeof(w)

Reduce the components of the vector `w` to range [0.0, 1.0[, incorporating a tolerance
`atol` in the reduction.
"""
function reduce_translation_to_unitrange(w::AbstractVector{<:Real}; atol=DEFAULT_ATOL)
    # naïve approach to achieve semi-robust reduction of integer-translation via a slightly
    # awful "approximate" modulo approach; basically just the equivalent of mod.(w′,1.0),
    # but reducing in a range DEFAULT_ATOL around each value.
    w′ = similar(w)
    for (i, wᵢ) in enumerate(w)
        w′ᵢ = mod(wᵢ, one(eltype(w)))
        if _fastpath_atol_isapprox(round(w′ᵢ), w′ᵢ; atol)
            # mod.(w, 1.0) occassionally does not reduce values that are very nearly 1.0
            # due to floating point errors: we use a tolerance here to round everything 
            # close to 0.0 or 1.0 exactly to 0.0
            w′ᵢ = zero(eltype(w))
        end
        @inbounds w′[i] = w′ᵢ
    end
    return typeof(w)(w′)
end

"""
    (⊚)(op1::T, op2::T) where T<:SymOperation -->  Vector{Float64}

Compose two symmetry operations `op1` ``= \\{W₁|w₁\\}`` and `op2` ``= \\{W₂|w₂\\}`` and
return the quotient of ``w₁+W₁w₂`` and 1. This functionality complements
`op1*op2`, which yields the translation modulo 1; accordingly, 
`translation(op1*op2) + op1⊚op2` yields the translation component
of the composition `op1` and `op2` *without* taking it modulo 1,
i.e. including any "trivial" lattice translation.

Note that ⊚ can be auto-completed in Julia via \\circledcirc+[tab]
""" 
function (⊚)(op1::T, op2::T) where T<:SymOperation
    # Translation result _without_ taking `mod`
    w′ = translation(op1) + rotation(op1)*translation(op2)  
    # Then we take w′ modulo lattice vectors
    w′′ = reduce_translation_to_unitrange(w′)
    # Then we subtract the two
    w′′′ = w′ - w′′
    return w′′′
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

    return T(W⁻¹, w⁻¹)
end


"""
    MultTable(ops::AbstractVector, modτ=true)

Compute the multiplication (or Cayley) table of `ops`, an iterable of `SymOperation`s.

A `MultTable` is returned, which contains symmetry operations resulting from composition 
of `row` and `col` operators; the table of indices give the symmetry operators relative to
the ordering of `ops`.

## Keyword arguments
- `modτ` (default: `true`): whether composition of operations is taken modulo lattice
vectors (`true`) or not (`false`).
"""
function MultTable(ops; modτ::Bool=true)
    N = length(ops)
    table = Matrix{Int}(undef, N,N)
    for (row,oprow) in enumerate(ops)
        for (col,opcol) in enumerate(ops)
            op′ = compose(oprow, opcol, modτ)
            match = findfirst(≈(op′), ops)
            if isnothing(match)
                throw(DomainError(ops, "provided operations do not form a group; group does not contain $op′"))
            end
            @inbounds table[row,col] = match
        end
    end
    return MultTable(ops, table)
end


function check_multtable_vs_ir(lgir::LGIrrep{D}, αβγ=nothing) where D
    ops = operations(lgir)
    sgnum = num(lgir); cntr = centering(sgnum, D)
    primitive_ops = primitivize.(ops, cntr) # must do multiplication table in primitive basis, cf. choices in `compose`
    check_multtable_vs_ir(MultTable(primitive_ops), lgir, αβγ)
end
function check_multtable_vs_ir(mt::MultTable, ir::AbstractIrrep, αβγ=nothing;
            verbose::Bool=false)
    havewarned = false
    Ds = ir(αβγ)
    ops = operations(ir)
    if ir isa LGIrrep
        k = position(ir)(αβγ)
    end
    N = length(ops)

    checked = trues(N, N)
    for (i,Dⁱ) in enumerate(Ds)     # rows
        for (j,Dʲ) in enumerate(Ds) # cols
            @inbounds mtidx = mt.table[i,j]
            if iszero(mtidx) && !havewarned
                @warn "Provided MultTable is not a group; cannot compare with irreps"
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
            # of `ϕ` here and in `(lgir::LGIrrep)(αβγ)`.
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
                if !havewarned && verbose
                    println("""Provided irreps do not match group multiplication table for group $(num(ir)) in irrep $(label(ir)):
                                First failure at (row,col) = ($(i),$(j));
                                Expected idx $(mtidx), got idx $(findall(≈(Dⁱʲ), Ds))""")
                    print("Expected irrep = ")
                    if ir isa LGIrrep
                        println(cis(ϕ)*Ds[mtidx])
                    else
                        println(Dⁱʲ)
                    end
                    println("Got irrep      = $(Dⁱʲ)")
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
function littlegroup(ops::AbstractVector{SymOperation{D}}, kv::KVec{D},
                cntr::Char='P') where D
    k₀, kabc = parts(kv)
    checkabc = !iszero(kabc)
    idxlist = Int[]
    for (idx, op) in enumerate(ops)
        k₀′, kabc′ = parts(compose(op, kv, checkabc)) # this is k₀(𝐆)′ = [g(𝐑)ᵀ]⁻¹k₀(𝐆)  
        diff = k₀′ .- k₀
        diff = primitivebasismatrix(cntr, Val(D))'*diff 
        kbool = all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), diff) # check if k₀ and k₀′ differ by a _primitive_ reciprocal vector
        abcbool = checkabc ? isapprox(kabc′, kabc, atol=DEFAULT_ATOL) : true # check if kabc == kabc′; no need to check for difference by a reciprocal vec, since kabc is in interior of BZ

        if kbool && abcbool # ⇒ part of little group
            push!(idxlist, idx)
        end
    end
    return idxlist, view(ops, idxlist)
end

"""
$(TYPEDSIGNATURES)

Return the little group associated with space group `sg` at the **k**-vector `kv`.

Optionally, an associated **k**-vector label `klab` can be provided; if not provided, the
empty string is used as label.
"""
function littlegroup(sg::SpaceGroup, kv::KVec, klab::String="") 
    _, lgops = littlegroup(operations(sg), kv, centering(sg))
    return LittleGroup{dim(sg)}(num(sg), kv, klab, lgops)
end


# TODO: Generalize docstring to actual (broader) signature
"""
    orbit(g::AbstractVector{<:SymOperation}, kv::KVec, cntr::Char)  -->  Vector{KVec{D}}
    orbit(lg::LittleGroup)
    orbit(lgir::LGIrrep)

Return the orbit of of the reciprocal-space vector `kv` under the action of the group
`g`, also known as the star of **k**.

The orbit of `kv` in `g` is the set of inequivalent **k**-points obtained by composition
of all the symmetry operations of `g` with `kv`. Two reciprocal vectors ``\\mathbf{k}``
and ``\\mathbf{k}'`` are equivalent if they differ by a primitive reciprocal lattice vector.

If `kv` and `g` are specified in a conventional basis but refer to a non-primitive
lattice, the centering type `cntr` must be provided to ensure that only equivalence by
primitive (not conventional) reciprocal lattice vectors are considered.
If the centering type of the group `g` can be inferred from `g` (e.g., if `g` is a
`SpaceGroup`), `orbit` will assume a conventional setting and use the inferred centering
type; otherwise, if `cntr` is neither explicitly set nor inferrable, a primitive setting is
assumed.
"""
function orbit(g::AbstractVector{SymOperation{D}},
               v::Union{AbstractVec{D}, Bravais.AbstractPoint{D}},
               P::Union{Nothing, AbstractMatrix{<:Real}} = nothing) where D
    vs = [v]
    for op in g
        v′ = op*v
        if !isapproxin(v′, vs, P, #=modw=#true)
            push!(vs, v′)
        end
    end
    return vs
end
function orbit(g::AbstractVector{SymOperation{D}},
               v::Union{AbstractVec{D}, Bravais.AbstractPoint{D}},
               cntr::Char) where D
    return orbit(g, v, primitivebasismatrix(cntr, Val(D)))
end
orbit(sg::SpaceGroup{D}, kv::KVec{D}) where D = orbit(sg, kv, centering(sg))

"""
    cosets(G, H)

For a subgroup `H` ``=H`` of `G` = ``G``, find a set of (left) coset representatives 
``\\{g_i\\}`` of ``H`` in ``G``, such that (see e.g., Inui et al., Section 2.7)
```math
    G = \\bigcup_i g_i H.
```
The identity operation ``1`` is always included in ``\\{g_i\\}``.

## Example
```jldoctest cosets
julia> G = pointgroup("6mm");

julia> H = pointgroup("3");

julia> Q = cosets(G, H)
4-element Vector{SymOperation{3}}:
 1
 2₀₀₁
 m₁₁₀
 m₁₋₁₀
```
To generate the associated cosets, simply compose the representatives with `H`:
```jldoctest cosets
julia> [compose.(Ref(q), H) for q in Q]
4-element Vector{Vector{SymOperation{3}}}:
 [1, 3₀₀₁⁺, 3₀₀₁⁻]
 [2₀₀₁, 6₀₀₁⁻, 6₀₀₁⁺]
 [m₁₁₀, m₁₀₀, m₀₁₀]
 [m₁₋₁₀, m₁₂₀, m₂₁₀]
```
"""
function cosets(
            G::AbstractVector{T},
            H::AbstractVector{T}
            ) where T<:AbstractOperation

    iszero(rem(length(G), length(H))) || error("H must be a subgroup of G: failed Lagrange's theorem")
    ind = div(length(G), length(H)) # [H:G]
    
    representatives = [one(T)]
    _cosets         = Vector{T}[operations(H)]
    sizehint!(representatives, ind)
    sizehint!(_cosets, ind)
    for g in G
        any(_coset -> isapproxin(g, _coset), _cosets) && continue
        
        push!(representatives, g)
        push!(_cosets, compose.(Ref(g), H))

        length(representatives) == ind && break
    end
    length(representatives) == ind || error("failed to find a set of coset representatives")
    return representatives
end

@doc raw"""
    compose(op::SymOperation, kv::KVec[, checkabc::Bool=true])  -->  KVec

Return the composition `op` ``= \{\mathbf{W}|\mathbf{w}\}`` and a reciprocal-space vector
`kv` ``= \mathbf{k}``.

The operation is assumed to be specified the direct lattice basis, and `kv` in the
reciprocal lattice basis. If both were specified in the Cartesian basis, `op` would act
directly via its rotation part. However, because of the different bases, it now acts as:
```math
    \mathbf{k}' = (\mathbf{W}^{\text{T}})^{-1}\mathbf{k}.
```
Note the transposition _and_ inverse of ``\mathbf{W}``, arising as a result of the implicit
real-space basis of ``\{\mathbf{W}|\mathbf{w}\}`` versus the reciprocal-space basis of
``\mathbf{k}`` and the requirement that 
``\mathbf{k} \cdot \mathbf{r} = \mathbf{k}' \cdot \mathbf{r}'`` for rotated (no translation,
``\mathbf{w} = \mathbf{0}``) vectors ``\mathbf{k}'`` and ``\mathbf{r}'``.

Note also that the composition of ``\{\mathbf{W}|\mathbf{w}\}`` with ``\mathbf{k}`` is
independent of ``\mathbf{w}``, i.e., translations do not act in reciprocal space.

## Extended help

If `checkabc = false`, the free part of `KVec` is not transformed (can be improve 
performance in situations when `kabc` is zero, and several transformations are requested).
"""
@inline function compose(op::SymOperation{D}, kv::KVec{D}, checkabc::Bool=true) where D
    k₀, kabc = parts(kv)
    W⁻¹ᵀ = transpose(inv(rotation(op)))
    k₀′ = W⁻¹ᵀ*k₀
    kabc′ = checkabc ? W⁻¹ᵀ * kabc : kabc
    return KVec{D}(k₀′, kabc′)
end
function compose(op::SymOperation{D}, k::ReciprocalPoint{D}) where D
    W⁻¹ᵀ = transpose(inv(rotation(op)))
    return ReciprocalPoint{D}(W⁻¹ᵀ*parent(k))
end

@doc raw"""
    compose(op::SymOperation, rv::RVec)  -->  RVec

Return the composition of `op` ``= \{\mathbf{W}|\mathbf{w}\}`` and a real-space vector `rv`
``= \mathbf{r}``.

The operation is taken to act directly, returning
```math
    \mathbf{r}' = \{\mathbf{W}|\mathbf{w}\}\mathbf{r} = \mathbf{W}\mathbf{r} + \mathbf{w}.
```
The corresponding inverse action ``\{\mathbf{W}|\mathbf{w}\}^{-1}\mathbf{r} = 
\mathbf{W}^{-1}\mathbf{r} - \mathbf{W}^{-1}\mathbf{w}`` can be obtained via 
`compose(inv(op), rv)`.
"""
function compose(op::SymOperation{D}, rv::RVec{D}) where D
    cnst, free = parts(rv)
    W, w       = unpack(op)

    cnst′ = W*cnst + w
    free′ = W*free

    return RVec{D}(cnst′, free′)
end
function compose(op::SymOperation{D}, r::DirectPoint{D}) where D
    W, w = unpack(op)
    return DirectPoint{D}(W*parent(r) + w)
end

(*)(op::SymOperation{D}, v::Union{AbstractVec{D}, AbstractPoint{D}}) where D = compose(op, v)

"""
    primitivize(op::SymOperation, cntr::Char, modw::Bool=true) --> typeof(op)

Return a symmetry operation `op′` ``≡ \\{W'|w'\\}`` in a primitive setting, transformed
from an input symmetry operation `op` ``= \\{W|w\\}`` in a conventional setting.
The operations ``\\{W'|w'\\}`` and ``\\{W|w\\}`` are related by a transformation
``\\{P|p\\}`` via (cf. Section 1.5.2.3 of [^ITA6]):

```math
\\{W'|w'\\} = \\{P|p\\}⁻¹\\{W|w\\}\\{P|p\\}.
```

where ``P`` and ``p`` are the basis change matrix and origin shifts, respectively.
The relevant transformation ``\\{P|p\\}`` is inferred from the centering type, as provided
by `cntr` (see also [`Bravais.centering`](@ref)).

By default, translation parts of `op′`, i.e. ``w'`` are reduced modulo 1 (`modw = true`); to
disable this, set `modw = false`.

[^ITA6]: M.I. Aroyo, International Tables of Crystallography, Vol. A, 6th ed. (2016).
"""
function primitivize(op::SymOperation{D}, cntr::Char, modw::Bool=true) where D
    if (D == 3 && cntr === 'P') || (D == 2 && cntr === 'p')
        # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(op, P, nothing, modw)
    end
end

"""
    conventionalize(op′::SymOperation, cntr::Char, modw::Bool=true) --> typeof(op)

Return a symmetry operation `op` ``= \\{W|w\\}`` in a conventional setting, transformed
from an input symmetry operation `op′` ``≡ \\{W'|w'\\}`` in a primitive setting.

See [`primitivize(::SymOperation, ::Char, ::Bool)`](@ref) for description of the centering
argument `cntr` and optional argument `modw`.
"""
function conventionalize(op::SymOperation{D}, cntr::Char, modw::Bool=true) where D
    if (D == 3 && cntr === 'P') || (D == 2 && cntr === 'p')
        # primitive basis: identity-transform, short circuit
        return op
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(op, inv(P), nothing, modw)
    end
end

@doc raw""" 
    transform(op::SymOperation, P::AbstractMatrix{<:Real}, 
              p::Union{AbstractVector{<:Real}, Nothing}=nothing,
              modw::Bool=true)                           --> SymOperation

Transforms a `op` ``= \{\mathbf{W}|\mathbf{w}\}`` by a rotation matrix `P` and a translation
vector `p` (can be `nothing` for zero-translations), producing a new symmetry operation 
`op′` ``= \{\mathbf{W}'|\mathbf{w}'\}`` (cf. Section 1.5.2.3 of [^ITA6])

```math
\{\mathbf{W}'|\mathbf{w}'\} = \{\mathbf{P}|\mathbf{p}\}^{-1}\{\mathbf{W}|\mathbf{w}\}
\{\mathbf{P}|\mathbf{p}\}
```

with

```math
\mathbf{W}' = \mathbf{P}^{-1}\mathbf{W}\mathbf{P}
\text{ and }
\mathbf{w}' = \mathbf{P}^{-1}(\mathbf{w}+\mathbf{W}\mathbf{p}-\mathbf{p})
```

By default, the translation part of `op′`, i.e. ``\mathbf{w}'``, is reduced to the range
``[0,1)``, i.e. computed modulo 1. This can be disabled by setting `modw = false` (default,
`modw = true`).

See also [`Bravais.primitivize(::SymOperation, ::Char, ::Bool)`](@ref) and
[`Bravais.conventionalize(::SymOperation, ::Char, ::Bool)`](@ref).

[^ITA6]: M.I. Aroyo, International Tables of Crystallography, Vol. A, 6th ed. (2016).
"""
function transform(op::SymOperation{D}, P::AbstractMatrix{<:Real}, 
                   p::Union{AbstractVector{<:Real}, Nothing}=nothing,
                   modw::Bool=true) where D
    W′ = transform_rotation(op, P)             # = P⁻¹WP       (+ rounding)
    w′ = transform_translation(op, P, p, modw) # = P⁻¹(w+Wp-p)
                                               # with W ≡ rotation(op) and w ≡ translation(op)

    return SymOperation{D}(W′, w′)
end

function transform_rotation(op::SymOperation{D}, P::AbstractMatrix{<:Real}) where D
    W = rotation(op)
    W′ = P\(W*P)        # = P⁻¹WP

    # clean up rounding-errors introduced by transformation (e.g. 
    # occassionally produces -0.0). The rotational part will 
    # always have integer coefficients if it is in the conventional
    # or primitive basis of its lattice; if transformed to a nonstandard
    # lattice, it might not have that though.
    W′_cleanup = ntuple(Val(D*D)) do i
        @inbounds W′ᵢ = W′[i]
        rW′ᵢ = round(W′ᵢ)
        if !isapprox(W′ᵢ, rW′ᵢ, atol=DEFAULT_ATOL)
            rW′ᵢ = W′ᵢ # non-standard lattice transformation; fractional elements 
                       # (this is why we need Float64 in SymOperation{D})
        end
        # since round(x) takes positive values x∈[0,0.5] to 0.0 and negative
        # values x∈[-0.5,-0.0] to -0.0 -- and since it is bad for us to have
        # both 0.0 and -0.0 -- we convert -0.0 to 0.0 here
        rW′ᵢ === -zero(Float64) && (rW′ᵢ = zero(Float64))

        return W′ᵢ
    end

    if W′ isa SMatrix{D,D,Float64,D*D}
        return SMatrix{D,D,Float64,D*D}(W′_cleanup)
    else # P was not an SMatrix, so output isn't either
        return copyto!(W′, W′_cleanup)
    end
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
        return reduce_translation_to_unitrange(w′)
    else
        return w′
    end
end

# TODO: Maybe implement this in mutating form; lots of unnecessary allocations below in
#       many usecases
"""
    reduce_ops(ops::AbstractVector{SymOperation{D}},
               cntr::Char,
               conv_or_prim::Bool=true,
               modw::Bool=true) --> Vector{SymOperation{D}}

Reduce the operations `ops`, removing operations that are identical in the primitive basis
associated with the centering `cntr`. 

If `conv_or_prim = false`, the reduced operations are returned in the primitive basis
associated with `cntr`; otherwise, in the conventional.
If `modw = true`, the comparison in the primitive basis is done modulo unit primitive
lattice vectors; otherwise not.
A final argument of type `::Val{P}` can be specified to indicate a subperiodic group of
periodicity dimension `P`, different from the spatial embedding dimension `D`.
"""
function reduce_ops(ops::AbstractVector{SymOperation{D}}, cntr::Char, 
                    conv_or_prim::Bool=true, modw::Bool=true,
                    ::Val{Pdim}=Val(D) #= to allow subperiodic groups =#) where {D,Pdim}
    
    P = primitivebasismatrix(cntr, Val(D), Val(Pdim))
    # transform ops (equiv. to `primitivize.(ops, cntr, modw)` but avoids loading `P` anew
    # for each SymOperation
    ops′ = transform.(ops, Ref(P), nothing, modw)

    # remove equivalent operations
    ops′_reduced = SymOperation{D}.(uniquetol(matrix.(ops′), atol=Crystalline.DEFAULT_ATOL))

    if conv_or_prim # `true`: return in conventional basis
        return transform.(ops′_reduced, Ref(inv(P)), nothing, modw)
    else            # `false`: return in primitive basis
        return ops′_reduced
    end
end
@inline function reduce_ops(slg::Union{SpaceGroup, LittleGroup}, 
                            conv_or_prim::Bool=true, modw::Bool=true)
    return reduce_ops(operations(slg), centering(slg), conv_or_prim, modw)
end
function primitivize(sg::SpaceGroup, modw::Bool=true)
    return typeof(sg)(num(sg), reduce_ops(sg, false, modw))
end
function primitivize(lg::LittleGroup, modw::Bool=true)
    cntr = centering(lg)
    # transform both k-point and operations
    kv′  = primitivize(position(lg), cntr)
    ops′ = reduce_ops(operations(lg), cntr, false, modw)
    return typeof(lg)(num(lg), kv′, klabel(lg), ops′)
end

"""
    cartesianize(op::SymOperation{D}, Rs::DirectBasis{D}) --> SymOperation{D}

Converts `opˡ` from a lattice basis to a Cartesian basis, by computing the
transformed operators `opᶜ = 𝐑*opˡ*𝐑⁻¹` via the Cartesian basis matrix 𝐑 (whose columns are
the `DirectBasis` vectors `Rs[i]`). 

# Note 1
The matrix 𝐑 maps vectors coefficients in a lattice basis 𝐯ˡ to coefficients in a Cartesian
basis 𝐯ᶜ as 𝐯ˡ = 𝐑⁻¹𝐯ᶜ and vice versa as 𝐯ᶜ = 𝐑𝐯ˡ. Since a general transformation P 
transforms an "original" vectors with coefficients 𝐯 to new coefficients 𝐯′ via 𝐯′ = P⁻¹𝐯
and since we here here consider the lattice basis as the "original" basis we have P = 𝐑⁻¹. 
As such, the transformation of the operator `op` transforms as `opᶜ = P⁻¹*opˡ*P`, i.e.
`opᶜ = transform(opˡ,P) = transform(opˡ,𝐑⁻¹)`.

# Note 2
The display (e.g. Seitz and xyzt notation) of `SymOperation`s e.g. in the REPL implicitly
assumes integer coefficients for its point-group matrix: as a consequence, displaying 
`SymOperation`s in a Cartesian basis may produce undefined behavior. The matrix
representation remains valid, however.
"""
function cartesianize(op::SymOperation{D}, Rs::DirectBasis{D}) where D
    𝐑 = stack(Rs)
    # avoids computing inv(𝐑) by _not_ calling out to transform(opˡ, inv(𝐑))
    op′ = SymOperation{D}([𝐑*rotation(op)/𝐑 𝐑*translation(op)])
    return op′
end
function cartesianize(sg::SpaceGroup{D}, Rs::DirectBasis{D}) where D
    return SpaceGroup{D}(num(sg), cartesianize.(operations(sg), Ref(Rs)))
end

"""
    findequiv(op::SymOperation, ops::AbstractVector{SymOperation{D}}, cntr::Char) 
                                                --> Tuple{Int, Vector{Float64}}

Search for an operator `op′` in `ops` which is equivalent, modulo differences
by *primitive* lattice translations `Δw`, to `op`. Return the index of `op′` in 
`ops`, as well as the primitive translation difference `Δw`. If no match is found
returns `nothing`.

The small irreps of `op` at wavevector k, Dⱼᵏ[`op`], can be computed from 
the small irreps of `op′`, Dⱼᵏ[`op′`], via Dⱼᵏ[`op`] = exp(2πik⋅`Δw`)Dⱼᵏ[`op′`]
"""
function findequiv(op::SymOperation{D}, ops::AbstractVector{SymOperation{D}},
            cntr::Char) where D
    W = rotation(op)
    w = translation(op)

    P = primitivebasismatrix(cntr, Val(D))
    w′ = P\w # `w` in its primitive basis

    for (j, opⱼ) in enumerate(ops)
        Wⱼ = rotation(opⱼ)
        wⱼ = translation(opⱼ)
        wⱼ′ = P\wⱼ

        if W == Wⱼ # rotation-part of op and opⱼ is identical
            # check if translation-part of op and opⱼ is equivalent, modulo a primitive lattice translation
            if all(el -> isapprox(el, round(el), atol=DEFAULT_ATOL), w′.-wⱼ′)
                return j, w.-wⱼ
            end
        end
    end
    return nothing # didn't find any match
end


"""
    _findsubgroup(opsᴳ::T, opsᴴ::T′[, cntr]) where T⁽′⁾<:AbstractVector{SymOperation}
                                                    --> Tuple{Bool, Vector{Int}}

Returns a 2-tuple with elements:

1. A boolean, indicating whether the group ``H`` (with operators `opsᴴ`) is a subgroup of
   the group ``G`` (with operators `opsᴳ`), i.e. whether ``H < G``.
2. An indexing vector `idxsᴳ²ᴴ` of `opsᴳ` into `opsᴴ` (empty if `H` is not a subgroup of 
   `G`), such that `all(isapprox(opsᴳ[idxsᴳ²ᴴ], opsᴴ, cntr)`.

If `cntr` is provided, comparison of operations in ``G`` and ``H`` is done in the associated
primitive basis (i.e., `cntr` is passed to `isapprox` which compares operations in ``G`` and
``H``).
"""
function _findsubgroup(opsᴳ::AbstractVector{SymOperation{D}},
                       opsᴴ::AbstractVector{SymOperation{D}},
                       cntr::Union{Char, Nothing}=nothing) where D
    idxsᴳ²ᴴ = Vector{Int}(undef, length(opsᴴ))
    @inbounds for (idxᴴ, opᴴ) in enumerate(opsᴴ)
        idxᴳ = findfirst(opᴳ -> isapprox(opᴳ, opᴴ, cntr), opsᴳ)
        if idxᴳ !== nothing
            idxsᴳ²ᴴ[idxᴴ] = idxᴳ
        else
            return false, Int[]
        end
    end
    return true, idxsᴳ²ᴴ
end

"""
    issubgroup(opsᴳ::T, opsᴴ::T′) where T⁽′⁾<:AbstractVector{SymOperation} --> Bool

Determine whether the operations in group ``H`` are a subgroup of the group ``G`` (each with 
operations `opsᴳ` and `opsᴴ`, respectively), i.e. whether ``H<G``.
Specifically, this requires that ``G`` and ``H`` are both groups and that for every ``h∈H``
there exists an element ``g∈G`` such that ``h=g``.

Returns a Boolean answer (`true` if normal, `false` if not).

## Note
This compares space groups rather than space group types, i.e. the comparison assumes a
matching setting choice between ``H`` and ``G``. To compare space group types with different
conventional settings, they must first be transformed to a shared setting.
"""
function issubgroup(opsᴳ::AbstractVector{SymOperation{D}},
                    opsᴴ::AbstractVector{SymOperation{D}}) where D

    _subgroup_fastpath_checks(length(opsᴳ), length(opsᴴ)) || return false
    return _is_setting_matched_subgroup(opsᴳ, opsᴴ)
end

function _is_setting_matched_subgroup(opsᴳ::AbstractVector{SymOperation{D}},
                                      opsᴴ::AbstractVector{SymOperation{D}}) where D
    for h in opsᴴ
        found = false
        for g in opsᴳ
            # consider two operations identical if they differ by an integer translation         
            Δw = translation(h) - translation(g)
            if isapprox(rotation(h), rotation(g), atol=DEFAULT_ATOL) &&
               isapprox(Δw, round.(Δw), atol=DEFAULT_ATOL)
                found = true
                continue
            end
        end
        found || return false
    end
    return true
end

function _subgroup_fastpath_checks(Nᴳ::Integer, Nᴴ::Integer)
    # "fast-path" checks for H !< G (returns false if violated, true as sentinel otherwise)
    Nᴳ > Nᴴ || return false   # 1. Must be a proper supergroup (i.e. higher order)
    return rem(Nᴳ, Nᴴ) == 0   # 2. Must obey Lagrange's theorem (index is an integer)
end

"""
    isnormal(opsᴳ::AbstractVector{<:SymOperation},
             opsᴴ::AbstractVector{<:SymOperation};
             verbose::Bool=false)                    --> Bool

Determine whether the operations in group ``H`` are normal in the group ``G`` (each with 
operations `opsᴳ` and `opsᴴ`), in the sense that 
    
```math
ghg⁻¹ ∈ H, ∀ g∈G, ∀ h∈H
```

Returns a Boolean answer (`true` if normal, `false` if not).

## Note 
This compares space groups rather than space group types, i.e. the comparison assumes a
matching setting choice between ``H`` and ``G``. To compare space group types with
different conventional settings, they must first be transformed to a shared setting.
"""
function isnormal(opsᴳ::AbstractVector{SymOperation{D}},
                  opsᴴ::AbstractVector{SymOperation{D}};
                  verbose::Bool=false) where D
    for g in opsᴳ
        g⁻¹ = inv(g)
        for h in opsᴴ
            # check if ghg⁻¹ ∉ G
            h′ = g*h*g⁻¹
            if !isapproxin(h′, opsᴴ, nothing; atol=Crystalline.DEFAULT_ATOL)
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

"""
$(TYPEDSIGNATURES)

Return the group generated from a finite set of generators `gens`.

## Keyword arguments
- `cntr` (default, `nothing`): check equivalence of operations modulo primitive lattice
  vectors (see also `isapprox(::SymOperation, ::SymOperation, cntr::Union{Nothing, Char})`; 
  only nonequivalent operations are included in the returned group.
- `modτ` (default, `true`): the group composition operation can either be taken modulo
  lattice vectors (`true`) or not (`false`, useful e.g. for site symmetry groups). In this
  case, the provided generators will also be taken modulo integer lattice translations.
- `Nmax` (default, `256`): the maximum size of the generated group. This is essentially
  a cutoff set to ensure halting of execution in case the provided set of generators do not
  define a *finite* group (especially relevant if `modτ=false`). If more operations than
  `Nmax` are generated, the method throws an overflow error.
"""
function generate(gens::AbstractVector{SymOperation{D}};
                  cntr::Union{Nothing,Char}=nothing,
                  modτ::Bool = true,
                  Nmax::Integer = 256) where D
    ops = if modτ
        [SymOperation{D}(op.rotation,
                         reduce_translation_to_unitrange(translation(op))) for op in gens]
    else
        collect(gens)
    end
    unique!(ops)
    Ngens = length(ops)
    
    # FIXME: there's probably a more efficient way to do this?
    # Ideas: apparently, the group order can be inferred from the generators without
    #        needing to compute the group itself (Schreier-Sims algorithm, see GAP; &
    #        also the Orbit-Stabilizer theorem; https://math.stackexchange.com/a/1706006).
    while true
        Nₒₚ = length(ops)
        for opᵢ in (@view ops[OneTo(Nₒₚ)])
            for opⱼ in (@view ops[OneTo(Ngens)]) # every op can be written as products w/ original generators
                opᵢⱼ = compose(opᵢ, opⱼ, modτ)
                if !isapproxin(opᵢⱼ, ops, cntr, modτ)
                    push!(ops, opᵢⱼ)
                    # early out if generators don't seem to form a closed group ...
                    length(ops) > Nmax && _throw_overflowed_generation()
                end
            end
        end
        Nₒₚ == length(ops) && return GenericGroup{D}(ops)
    end
end

_throw_overflowed_generation() = 
    throw(OverflowError("The provided set of generators overflowed Nmax distinct "*
                        "operations: generators may not form a finite group; "*
                        "otherwise, try increasing Nmax"))