"""
    schoenflies(sgnum::Integer) --> String

Returns the Schoenflies notation for a given space group number
`sgnum`. Schoenflies notation only applies to point groups and 
space groups, not plane groups, so this notation is only relevant
in three dimensions.
"""
schoenflies(sgnum::Integer) = SCHOENFLIES_TABLE[sgnum]
schoenflies(sg::SpaceGroup{3}) = schoenflies(num(sg))

"""
    iuc(sgnum::Integer, D::Integer=3) --> String

Returns the IUC (International Union of Crystallography) notation for space group number
`sgnum` and dimensionality `D`, as used in the International Tables of Crystallography. 
The notation is sometimes also known as the Hermann-Mauguin notation; the functionality
is consequently aliased by `hermannmauguin(sgnum, D)`. 
IUC/Hermann-Mauguin notation applies in one, two, and three-dimensions.

For more information, see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.
"""
@inline iuc(sgnum::Integer, D::Integer=3) = SGS_IUC_NOTATION[D][sgnum]
@inline iuc(sg::Union{SpaceGroup{D},LittleGroup{D}}) where D = iuc(num(sg), D)
const hermannmauguin = iuc # alias

""" 
    centering(sg::SpaceGroup) --> Char
    centering(sgnum::Integer, D::Integer=3) --> Char

Determines the conventional centering type of a given space/plane group `sg` (alternatively
specified by its conventional number `sgnum` and dimensionality `D` by comparison with the
Hermann-Mauguin notation's first letter. 

Possible output values, depending on dimensionality `D`, are (see ITA Sec. 9.1.4):

    D=2 ┌ 'p': no centring (primitive)
        └ 'c': face centered

    D=3 ┌ 'P': no centring (primitive)
        ├ 'I': body centred (innenzentriert)
        ├ 'F': all-face centred
        ├ 'A', 'B', 'C': one-face centred, (b,c) or (c,a) or (a,b)
        └ 'R': hexagonal cell rhombohedrally centred
"""
centering(sgnum::Integer, D::Integer=3) = first(iuc(sgnum, D))
centering(sg::Union{SpaceGroup{D},LittleGroup{D}}) where D = first(centering(num(sg), D))

# Schoenflies notation, ordered relative to space group number
# [from https://bruceravel.github.io/demeter/artug/atoms/space.html]
const SCHOENFLIES_TABLE = (
# triclinic
"C₁¹",    "Cᵢ¹",
# monoclinic
"C₂¹",    "C₂²",    "C₂³",    "Cₛ¹",    "Cₛ²",    "Cₛ³",
"Cₛ⁴",    "C₂ₕ¹",   "C₂ₕ²",   "C₂ₕ³",   "C₂ₕ⁴",   "C₂ₕ⁵",
"C₂ₕ⁶",
# orthorhombic
"D₂¹",    "D₂²",    "D₂³",    "D₂⁴",    "D₂⁵",    "D₂⁶",
"D₂⁷",    "D₂⁸",    "D₂⁹",    "C₂ᵥ¹",   "C₂ᵥ²",   "C₂ᵥ³",
"C₂ᵥ⁴",   "C₂ᵥ⁵",   "C₂ᵥ⁶",   "C₂ᵥ⁷",   "C₂ᵥ⁸",   "C₂ᵥ⁹",
"C₂ᵥ¹⁰",  "C₂ᵥ¹¹",  "C₂ᵥ¹²",  "C₂ᵥ¹³",  "C₂ᵥ¹⁴",  "C₂ᵥ¹⁵",
"C₂ᵥ¹⁶",  "C₂ᵥ¹⁷",  "C₂ᵥ¹⁸",  "C₂ᵥ¹⁹",  "C₂ᵥ²⁰",  "C₂ᵥ²¹",
"C₂ᵥ²²",  "D₂ₕ¹",   "D₂ₕ²",   "D₂ₕ³",   "D₂ₕ⁴",   "D₂ₕ⁵",
"D₂ₕ⁶",   "D₂ₕ⁷",   "D₂ₕ⁸",   "D₂ₕ⁹",   "D₂ₕ¹⁰",  "D₂ₕ¹¹",
"D₂ₕ¹²",  "D₂ₕ¹³",  "D₂ₕ¹⁴",  "D₂ₕ¹⁵",  "D₂ₕ¹⁶",  "D₂ₕ¹⁷",
"D₂ₕ¹⁸",  "D₂ₕ¹⁹",  "D₂ₕ²⁰",  "D₂ₕ²¹",  "D₂ₕ²²",  "D₂ₕ²³",
"D₂ₕ²⁴",  "D₂ₕ²⁵",  "D₂ₕ²⁶",  "D₂ₕ²⁷",  "D₂ₕ²⁸",
# tetragonal
"C₄¹",    "C₄²",    "C₄³",    "C₄⁴",    "C₄⁵",    "C₄⁶",
"S₄¹",    "S₄²",    "C₄ₕ¹",   "C₄ₕ²",   "C₄ₕ³",   "C₄ₕ⁴",
"C₄ₕ⁵",   "C₄ₕ⁶",   "D₄¹",    "D₄²",    "D₄³",    "D₄⁴",
"D₄⁵",    "D₄⁶",    "D₄⁷",    "D₄⁸",    "D₄⁹",    "D₄¹⁰",
"C₄ᵥ¹",   "C₄ᵥ²",   "C₄ᵥ³",   "C₄ᵥ⁴",   "C₄ᵥ⁵",   "C₄ᵥ⁶",
"C₄ᵥ⁷",   "C₄ᵥ⁸",   "C₄ᵥ⁹",   "C₄ᵥ¹⁰",  "C₄ᵥ¹¹",  "C₄ᵥ¹²",
"D₂d¹",   "D₂d²",   "D₂d³",   "D₂d⁴",   "D₂d⁵",   "D₂d⁶",
"D₂d⁷",   "D₂d⁸",   "D₂d⁹",   "D₂d¹⁰",  "D₂d¹¹",  "D₂d¹²",
"D₄ₕ¹",   "D₄ₕ²",   "D₄ₕ³",   "D₄ₕ⁴",   "D₄ₕ⁵",   "D₄ₕ⁶",
"D₄ₕ⁷",   "D₄ₕ⁸",   "D₄ₕ⁹",   "D₄ₕ¹⁰",  "D₄ₕ¹¹",  "D₄ₕ¹²",
"D₄ₕ¹³",  "D₄ₕ¹⁴",  "D₄ₕ¹⁵",  "D₄ₕ¹⁶",  "D₄ₕ¹⁷",  "D₄ₕ¹⁸",
"D₄ₕ¹⁹",  "D₄ₕ²⁰",
# trigonal
"C₃¹",    "C₃²",    "C₃³",    "C₃⁴",    "C₃ᵢ¹",   "C₃ᵢ²",
"D₃¹",    "D₃²",    "D₃³",    "D₃⁴",    "D₃⁵",    "D₃⁶",
"D₃⁷",    "C₃ᵥ¹",   "C₃ᵥ²",   "C₃ᵥ³",   "C₃ᵥ⁴",   "C₃ᵥ⁵",
"C₃ᵥ⁶",   "D₃d¹",   "D₃d²",   "D₃d³",   "D₃d⁴",   "D₃d⁵",
"D₃d⁶",
# hexagonal
"C₆¹",    "C₆²",    "C₆³",    "C₆⁴",    "C₆⁵",    "C₆⁶",
"C₃ₕ¹",   "C₆ₕ¹",   "C₆ₕ²",   "D₆¹",    "D₆²",    "D₆³",
"D₆⁴",    "D₆⁵",    "D₆⁶",    "C₆ᵥ¹",   "C₆ᵥ²",   "C₆ᵥ³",
"C₆ᵥ⁴",   "D₃ₕ¹",   "D₃ₕ²",   "D₃ₕ³",   "D₃ₕ⁴",   "D₆ₕ¹",
"D₆ₕ²",   "D₆ₕ³",   "D₆ₕ⁴",
# cubic
"T¹",      "T²",    "T³",     "T⁴",    "T⁵",     "Tₕ¹",
"Tₕ²",     "Tₕ³",    "Tₕ⁴",    "Tₕ⁵",    "Tₕ⁶",    "Tₕ⁷",
"O¹",      "O²",    "O³",     "O⁴",    "O⁵",     "O⁶",
"O⁷",      "O⁸",    "Td¹",    "Td²",   "Td³",    "Td⁴",
"Td⁵",     "Td⁶",   "Oₕ¹",    "Oₕ²",    "Oₕ³",    "Oₕ⁴",
"Oₕ⁵",     "Oₕ⁶",    "Oₕ⁷",    "Oₕ⁸",    "Oₕ⁹",    "Oₕ¹⁰"
)

# IUC/Hermann-Mauguin notation, ordered relative to space/plane group number
const SGS_IUC_NOTATION = (
# ------------------------------------------------------------------------------------------
# line-group notation (one dimension) [see https://en.wikipedia.org/wiki/Line_group]
# ------------------------------------------------------------------------------------------    
("p1", "p1m"),
# ------------------------------------------------------------------------------------------
# plane-group notation (two dimensions) [see e.g. Table 19 of Cracknell, Adv. Phys. 1974, or
# https://www.cryst.ehu.es/cgi-bin/plane/programs/nph-plane_getgen?from=getwp]
# ------------------------------------------------------------------------------------------
(
# oblique
"p1",   "p2",
# rectangular ('p' or 'c' centering; c-centered lattices are rhombic in their primitive cell)
"p1m1", "p1g1", "c1m1", "p2mm", "p2mg", "p2gg", "c2mm",   
# square
"p4",   "p4mm", "p4gm",
# hexagonal
"p3",   "p3m1", "p31m", "p6",   "p6mm"
),
# ------------------------------------------------------------------------------------------
# space-group notation (three dimensions) following the conventions of ITA and Bilbao:
# https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen
# ------------------------------------------------------------------------------------------
(
# triclinic
"P1",      "P-1",
# monoclinic
"P2",      "P2₁",     "C2",      "Pm",      "Pc",      "Cm",
"Cc",      "P2/m",    "P2₁/m",   "C2/m",    "P2/c",    "P2₁/c",
"C2/c",
# orthorhombic
"P222",    "P222₁",   "P2₁2₁2",  "P2₁2₁2₁", "C222₁",   "C222",
"F222",    "I222",    "I2₁2₁2₁", "Pmm2",    "Pmc2₁",   "Pcc2",
"Pma2",    "Pca2₁",   "Pnc2",    "Pmn2₁",   "Pba2",    "Pna2₁",
"Pnn2",    "Cmm2",    "Cmc2₁",   "Ccc2",    "Amm2",    "Aem2",
"Ama2",    "Aea2",    "Fmm2",    "Fdd2",    "Imm2",    "Iba2",
"Ima2",    "Pmmm",    "Pnnn",    "Pccm",    "Pban",    "Pmma",
"Pnna",    "Pmna",    "Pcca",    "Pbam",    "Pccn",    "Pbcm",
"Pnnm",    "Pmmn",    "Pbcn",    "Pbca",    "Pnma",    "Cmcm",
"Cmce",    "Cmmm",    "Cccm",    "Cmme",    "Ccce",    "Fmmm",
"Fddd",    "Immm",    "Ibam",    "Ibca",    "Imma",
# tetragonal
"P4",      "P4₁",     "P4₂",     "P4₃",     "I4",      "I4₁",
"P-4",     "I-4",     "P4/m",    "P4₂/m",   "P4/n",    "P4₂/n",
"I4/m",    "I4₁/a",   "P422",    "P42₁2",   "P4₁22",   "P4₁2₁2",
"P4₂22",   "P4₂2₁2",  "P4₃22",   "P4₃2₁2",  "I422",    "I4₁22",
"P4mm",    "P4bm",    "P4₂cm",   "P4₂nm",   "P4cc",    "P4nc",
"P4₂mc",   "P4₂bc",   "I4mm",    "I4cm",    "I4₁md",   "I4₁cd",
"P-42m",   "P-42c",   "P-42₁m",  "P-42₁c",  "P-4m2",   "P-4c2",
"P-4b2",   "P-4n2",   "I-4m2",   "I-4c2",   "I-42m",   "I-42d",
"P4/mmm",  "P4/mcc",  "P4/nbm",  "P4/nnc",  "P4/mbm",  "P4/mnc",
"P4/nmm",  "P4/ncc",  "P4₂/mmc", "P4₂/mcm", "P4₂/nbc", "P4₂/nnm",
"P4₂/mbc", "P4₂/mnm", "P4₂/nmc", "P4₂/ncm", "I4/mmm",  "I4/mcm",
"I4₁/amd", "I4₁/acd",
# trigonal
"P3",      "P3₁",     "P3₂",     "R3",      "P-3",     "R-3",
"P312",    "P321",    "P3₁12",   "P3₁21",   "P3₂12",   "P3₂21",
"R32",     "P3m1",    "P31m",    "P3c1",    "P31c",    "R3m",
"R3c",     "P-31m",   "P-31c",   "P-3m1",   "P-3c1",   "R-3m",
"R-3c",
# hexagonal
"P6",      "P6₁",     "P6₅",     "P6₂",     "P6₄",     "P6₃",
"P-6",     "P6/m",    "P6₃/m",   "P622",    "P6₁22",   "P6₅22",
"P6₂22",   "P6₄22",   "P6₃22",   "P6mm",    "P6cc",    "P6₃cm",
"P6₃mc",   "P-6m2",   "P-6c2",   "P-62m",   "P-62c",   "P6/mmm",
"P6/mcc",  "P6₃/mcm", "P6₃/mmc",
# cubic
"P23",     "F23",     "I23",     "P2₁3",    "I2₁3",    "Pm-3",
"Pn-3",    "Fm-3",    "Fd-3",    "Im-3",    "Pa-3",    "Ia-3",
"P432",    "P4₂32",   "F432",    "F4₁32",   "I432",    "P4₃32",
"P4₁32",   "I4₁32",   "P-43m",   "F-43m",   "I-43m",   "P-43n",
"F-43c",   "I-43d",   "Pm-3m",   "Pn-3n",   "Pm-3n",   "Pn-3m",
"Fm-3m",   "Fm-3c",   "Fd-3m",   "Fd-3c",   "Im-3m",   "Ia-3d"
)
)



""" 
    seitz(op::SymOperation) --> String

Computes the correponding Seitz notation for a symmetry operation in triplet/xyzt form.

Implementation based on ITA5 Table 11.2.1.1, with 3D point group parts inferred from
the trace and determinant of the matrix ``W`` in the triplet ``{W|w}``.


| detW\\trW | -3 | -2 | -1 | 0  | 1 | 2 | 3 |
|-----------|----|----|----|----|---|---|---|
|  1        |    |    |  2 | 3  | 4 | 6 | 1 |
|  -1       | -1 | -6 | -4 | -3 | m |   |   |

with the elements of the table giving the type of symmetry operation in
in Hermann-Mauguin notation. The rotation axis and the rotation sense are 
computed following the rules in ITA6 Sec. 1.2.2.4(1)(b-c).
The implementation has been checked against the Tables 1.4.2.1-5 of ITA6.

Note that the orientation of axis (i.e. its sign) is not necessarily equal
to the orientation picked in those tables; it is a matter of convention,
and the conventions have not been explicated in ITA6.

2D operations are treated by the same procedure, by elevation in a third dimension; 1D
operations by a simple inspection of sign.
"""
function seitz(op::SymOperation{D}) where D
    w = translation(op)
    if D == 3
        W = rotation(op)
    elseif D == 2 # augment 2D case by "adding" an invariant z dimension
        W′ = rotation(op)
        W  = @inbounds SMatrix{3,3,Float64,9}( # build by column (= [W′ zeros(2); 0 0 1])
                W′[1], W′[2], 0.0, W′[3], W′[4], 0.0, 0.0, 0.0, 1.0 )
    elseif D == 1
        W = rotation(op)
        isone(abs(W[1])) || throw(DomainError((W,w), "not a valid 1D symmetry operation"))
        W_str = signbit(W[1]) ? "-1" : "1"
        if iszero(w[1])
            return W_str
        else
            w_str = unicode_frac(w[1])
            return '{'*W_str*'|'*w_str*'}'
        end
    else
        throw(DomainError(D, "dimension different from 1, 2, or 3 is not supported"))
    end

    detW′ = det(W); detW = round(Int64, detW′) # det, then trunc & check
    isapprox(detW′, detW, atol=DEFAULT_ATOL) || throw(DomainError(detW′, "det W must be an integer for a SymOperation {W|w}"))
    trW′  = tr(W);  trW  = round(Int64, trW′)   # tr, then trunc & check
    isapprox(trW′,  trW,  atol=DEFAULT_ATOL) || throw(DomainError(trW′,  "tr W must be an integer for a SymOperation {W|w}"))

    io_pgop = IOBuffer()
    iszero(w) || print(io_pgop, '{')

    # --- rotation order (and proper/improper determination) ---
    rot = rotation_order_3d(detW, trW) # works for 2D also, since we augmented W above
    order = abs(rot)
    if rot == -2
        print(io_pgop, 'm')
    else
        print(io_pgop, rot)
    end
    
    if order ≠ 1
        # --- rotation axis (for order ≠ 1)---
        u = if D == 2 && rot ≠ -2   # only need orientation in 2D for mirrors 
            SVector{3,Int}(0, 0, 1)
        else
            uv = rotation_axis_3d(W, detW, order)
            SVector{3,Int}((uv[1], uv[2], uv[3]))
        end

        if !(D == 2 && rot ≠ -2)
            # (for 2D, ignore z-component)
            join(io_pgop, (subscriptify(string(uᵢ)) for uᵢ in u[SOneTo(D)]))
        end
        
        # --- rotation sense (for order > 2}) ---
        # ±-rotation sense is determined from sign of det(𝐙) where
        # 𝐙 ≡ [𝐮|𝐱|det(𝐖)𝐖𝐱] where 𝐱 is an arbitrary vector that 
        # is not parallel to 𝐮. [ITA6  Vol. A, p. 16, Sec. 1.2.2.4(1)(c)]
        if order > 2
            x = rand(-1:1, SVector{3, Int})
            while iszero(x×u) # check that generated 𝐱 is not parallel to 𝐮 (if it is, 𝐱×𝐮 = 0)
                x = rand(-1:1, SVector{3, Int}) 
            end
            Z = hcat(u, x, detW*(W*x))
            print(io_pgop, signbit(det(Z)) ? '⁻' : '⁺')
        end
    end

    # --- add translation for nonsymorphic operations ---
    if !iszero(w)
        print(io_pgop, '|')
        join(io_pgop, (unicode_frac(wᵢ) for wᵢ in w), ',')
        print(io_pgop, '}')
    end
    return String(take!(io_pgop))
end
seitz(str::String) = seitz(SymOperation(str))

"""
    rotation_order_3d(detW::Real, trW::Real) --> Int

Determine the integer rotation order of a 3D point group operation with a 3×3 matrix 
representation `W` (alternatively specified by its determinant `detW` and its trace `trW`).

The rotation order of
- Proper rotations is positive.
- Improper (mirrors, inversion, roto-inversions) is negative.
"""
function rotation_order_3d(detW::Real, trW::Real)
    if detW == 1 # proper rotations
        if -1 ≤ trW ≤ 1 # 2-, 3-, or 4-fold rotation
            rot = convert(Int, trW) + 3
        elseif trW == 2 # 6-fold rotation
            rot = 6
        elseif trW == 3 # identity operation
            rot = 1
        else 
            _throw_seitzerror(trW, detW)
        end
    elseif detW == -1 # improper rotations (rotoinversions)
        if trW == -3     # inversion
            rot = -1
        elseif trW == -2 # 6-fold rotoinversion
            rot = -6
        elseif -1 ≤ trW ≤ 0 # 4- and 3-fold rotoinversion
            rot = convert(Int, trW) - 3
        elseif trW == 1  # mirror, note that "m" == "-2" conceptually
            rot = -2
        else
            _throw_seitzerror(trW, detW)
        end
    else
        _throw_seitzerror(trW, detW)
    end
    
    return rot
end

"""
    rotation_order(W::Matrix{<:Real}) --> Int
    rotation_order(op::SymOperation)  --> Int

Determine the integer rotation order of a point group operation, input either as a matrix
`W` or `op::SymOperation`.

The rotation order of
- Proper rotations is positive.
- Improper (mirrors, inversion, roto-inversions) is negative.
"""
function rotation_order(W::AbstractMatrix{<:Real})
    if size(W) == (1,1)
        return convert(Int, W[1,1])
    elseif size(W) == (2,2)
        # we augment 2D case by effectively "adding" an invariant z dimension, extending
        # into 3D: then we just shortcut to what det and tr of that "extended" matrix is
        return rotation_order_3d(det(W), tr(W) + one(eltype(W)))
    elseif size(W) ≠ (3,3)
        throw(DomainError(size(W), "Point group operation must have a dimension ≤3"))
    end

    return rotation_order_3d(det(W), tr(W))
end
rotation_order(op::SymOperation) = rotation_order(rotation(op))

function rotation_axis_3d(W::AbstractMatrix{<:Real}, detW::Real, order::Integer)
    # the rotation axis 𝐮 of a 3D rotation 𝐖 of order k is determined from the product of 
    #       𝐘ₖ(𝐖) ≡ (d𝐖)ᵏ⁻¹+(d𝐖)ᵏ⁻² + ... + (d𝐖) + 𝐈 where d ≡ det(𝐖) 
    # with an arbitrary vector 𝐯 that is not perpendicular to 𝐮 [cf. ITA6  Vol. A, p. 16,
    # Sec. 1.2.2.4(1)(b)]

    order ≤ 0 && throw(DomainError(order, "order must be positive (i.e. not include sign)"))
    # if W is the identity or inversion, the notion of an axis doesn't make sense
    isone(order) && throw(DomainError(order, "order must be non-unity (i.e. operation must not be identity or inversion)"))

    Yₖ   = SMatrix{3,3,Float64}(I) # calculate Yₖ by iteration
    term = SMatrix{3,3,eltype(W)}(I)
    for j in OneTo(order-1)
        term = term*W # iteratively computes Wʲ
        if detW^j == -1;
            Yₖ = Yₖ - term 
        else
            Yₖ = Yₖ + term
        end
    end
    u′ = Yₖ*rand(SVector{3, Float64})
    while LinearAlgebra.norm(u′) < 1e-6
        # there is near-infinitesimal chance that u′ is zero for random v, but check anyway
        u′ = Yₖ*rand(SVector{3, Float64})
    end
    norm = minimum(Base.Filter(x->abs(x)>DEFAULT_ATOL, u′)) # minimum nonzero element
    u′ = u′/norm # normalize
    u  = round.(Int64, u′) # convert from float to integer and check validity of conversion
    if !isapprox(u′, u, atol=DEFAULT_ATOL)
        throw(DomainError(u′, "the rotation axis must be equivalent to an integer vector by appropriate normalization"))
    end
    # the sign of u is arbitrary: we adopt the convention of '-' elements coming "before"
    # '+' elements; e.g. [-1 -1 1] is picked over [1 1 -1] and [-1 1 -1] is picked over
    # [1 -1 1]; note that this impacts the sense of rotation which depends on the sign of
    # the rotation axis; finally, if all elements have the same sign (or zero), we pick a
    # positive overall sign ('+')
    if all(≤(0), u)
        u = -u
    else
        negidx = findfirst(signbit, u)
        firstnonzero = findfirst(≠(0), u) # don't need to bother taking abs, as -0 = 0 for integers (and floats)
        if negidx ≠ nothing && (negidx ≠ firstnonzero || negidx === firstnonzero === 3)
            u = -u 
        end
    end

    return u

end
rotation_axis_3d(W::AbstractMatrix)   = rotation_axis_3d(W, det(W), rotation_order(W))
rotation_axis_3d(op::SymOperation{3}) = (W=rotation(op); rotation_axis_3d(W, det(W), abs(rotation_order(W))))

_throw_seitzerror(trW, detW) = throw(DomainError((trW, detW), "trW = $(trW) for detW = $(detW) is not a valid symmetry operation; see ITA5 Vol A, Table 11.2.1.1"))