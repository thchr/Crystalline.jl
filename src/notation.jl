"""
    schoenflies(sgnum::Integer) --> String

Return the [Schoenflies notation](https://en.wikipedia.org/wiki/Schoenflies_notation) for
space group number `sgnum` in dimension 3.

Note that Schoenflies notation is well-defined only for 3D point and space groups.
"""
schoenflies(sgnum::Integer) = SCHOENFLIES_TABLE[sgnum]
schoenflies(sg::SpaceGroup{3}) = schoenflies(num(sg))

"""
    iuc(sgnum::Integer, D::Integer=3) --> String

Return the IUC (International Union of Crystallography) notation for space group number
`sgnum` in dimension `D` (1, 2, or 3), as used in the International Tables of
Crystallography.

The notation is sometimes also known as the
[Hermann-Mauguin notation](https://en.wikipedia.org/wiki/Hermann–Mauguin_notation).
"""
@inline function iuc(sgnum::Integer, D::Integer=3)
    if D==3
        @boundscheck (sgnum ∈ 1:230) || _throw_invalid_sgnum(sgnum, 3)
        return @inbounds SGS_IUC_NOTATION[3][sgnum]
    elseif D==2
        @boundscheck (sgnum ∈ 1:17) || _throw_invalid_sgnum(sgnum, 2)
        return @inbounds SGS_IUC_NOTATION[2][sgnum]
    elseif D==1
        @boundscheck (sgnum ∈ 1:2) || _throw_invalid_sgnum(sgnum, 1)
        return @inbounds SGS_IUC_NOTATION[1][sgnum]
    else
        _throw_invaliddim(D)
    end
end
@inline iuc(sg::Union{SpaceGroup{D},LittleGroup{D}}) where D = iuc(num(sg), D)

""" 
    centering(g::AbstractGroup) --> Char

Return the conventional centering type of a group. 

For groups without lattice structure (e.g., point groups), return `nothing`.
"""
centering(sg_or_lg::Union{SpaceGroup{D},LittleGroup{D}}) where D = centering(num(sg_or_lg), D)

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



@doc raw"""
    seitz(op::SymOperation) --> String

Computes the correponding Seitz notation for a symmetry operation in triplet/xyzt form.

Implementation based on ITA5 Table 11.2.1.1, with 3D point group parts inferred from
the trace and determinant of the matrix ``\mathb{W}`` in the triplet
``\{\mathbf{W}|\mathbf{w}\}``.

| detW/trW | -3 | -2 | -1 | 0  | 1 | 2 | 3 |
|:---------|----|----|----|----|---|---|---|
| **1**    |    |    |  2 | 3  | 4 | 6 | 1 |
| **-1**   | -1 | -6 | -4 | -3 | m |   |   |

with the elements of the table giving the type of symmetry operation in in Hermann-Mauguin
notation. The rotation axis and the rotation sense are computed following the rules in ITA6
Sec. 1.2.2.4(1)(b-c). See also .

Note that the orientation of the axis (i.e. its sign) does not necessarily match the
orientation picked in Tables 1.4.2.1-5 of ITA6; it is a matter of (arbitrary) convention,
and the conventions have not been explicated in ITA.

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

    detW′ = det(W); detW = round(Int, detW′) # det, then trunc & check
    isapprox(detW′, detW, atol=DEFAULT_ATOL) || throw(DomainError(detW′, "det W must be an integer for a SymOperation {W|w}"))
    trW′  = tr(W);  trW  = round(Int, trW′)   # tr, then trunc & check
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
            rotation_axis_3d(W, detW, order)
        end

        if !(D == 2 && rot ≠ -2)
            # (for 2D, ignore z-component)
            join(io_pgop, (subscriptify(string(u[i])) for i in SOneTo(D)))
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
    u  = round.(Int, u′) # convert from float to integer and check validity of conversion
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



# -----------------------------------------------------------------------------------------
# MULLIKEN NOTATION FOR POINT GROUP IRREPS

const PGIRLABS_CDML2MULLIKEN_3D = ImmutableDict(
    # sorted in ascending order wrt. Γᵢ CDML sorting; i.e. as 
    #       Γ₁, Γ₂, ... 
    #   or  Γ₁⁺, Γ₁⁻, Γ₂⁺, Γ₂⁻, ...
    # the association between CDMl and Mulliken labels are obtained obtained from
    # https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=spg
    # note that e.g., https://www.cryst.ehu.es/rep/point.html cannot be used, because the 
    # Γ-labels there do not always refer to the CDML convention; more likely, the B&C 
    # convention. For "setting = 2" cases, we used the `bilbao_pgs_url(..)` from the 
    # point group irrep crawl script
    # includes all labels in PGS_IUCs[3]
    "1"     => ImmutableDict("Γ₁"=>"A"),
    "-1"    => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ"),
    "2"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B"),
    "m"     => ImmutableDict("Γ₁"=>"A′", "Γ₂"=>"A′′"),
    "2/m"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"Bg", "Γ₂⁻"=>"Bᵤ"),
    "222"   => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B₁", "Γ₃"=>"B₃", "Γ₄"=>"B₂"),
    "mm2"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"B₂", "Γ₄"=>"B₁"),
    "mmm"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"B₁g", "Γ₂⁻"=>"B₁ᵤ", "Γ₃⁺"=>"B₃g", "Γ₃⁻"=>"B₃ᵤ", "Γ₄⁺"=>"B₂g", "Γ₄⁻"=>"B₂ᵤ"),
    "4"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃"=>"²E", "Γ₄"=>"¹E"),
    "-4"    => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃"=>"²E", "Γ₄"=>"¹E"),
    "4/m"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"Bg", "Γ₂⁻"=>"Bᵤ", "Γ₃⁺"=>"²Eg", "Γ₃⁻"=>"²Eᵤ", "Γ₄⁺"=>"¹Eg", "Γ₄⁻"=>"¹Eᵤ"),
    "422"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"B₁", "Γ₃"=>"A₂", "Γ₄"=>"B₂", "Γ₅"=>"E"),
    "4mm"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"B₁", "Γ₃"=>"B₂", "Γ₄"=>"A₂", "Γ₅"=>"E"),
    "-42m"  => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"B₁", "Γ₃"=>"B₂", "Γ₄"=>"A₂", "Γ₅"=>"E"), # setting 1
    "-4m2"  => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"B₁", "Γ₃"=>"A₂", "Γ₄"=>"B₂", "Γ₅"=>"E"), # setting 2 *** swapped B₂ and A₂; seems to be a typo in Bilbao? ***
    "4/mmm" => ImmutableDict("Γ₁⁺"=>"A₁g", "Γ₁⁻"=>"A₁ᵤ", "Γ₂⁺"=>"B₁g", "Γ₂⁻"=>"B₁ᵤ", "Γ₃⁺"=>"A₂g", "Γ₃⁻"=>"A₂ᵤ", "Γ₄⁺"=>"B₂g", "Γ₄⁻"=>"B₂ᵤ", "Γ₅⁺"=>"Eg", "Γ₅⁻"=>"Eᵤ"),
    "3"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"²E", "Γ₃"=>"¹E"),
    "-3"    => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"²Eg", "Γ₂⁻"=>"²Eᵤ", "Γ₃⁺"=>"¹Eg", "Γ₃⁻"=>"¹Eᵤ"),
    "312"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E"), # setting 1
    "321"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E"), # setting 2
    "3m1"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E"), # setting 1
    "31m"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E"), # setting 2
    "-31m"  => ImmutableDict("Γ₁⁺"=>"A₁g", "Γ₁⁻"=>"A₁ᵤ", "Γ₂⁺"=>"A₂g", "Γ₂⁻"=>"A₂ᵤ", "Γ₃⁺"=>"Eg", "Γ₃⁻"=>"Eᵤ"), # setting 1
    "-3m1"  => ImmutableDict("Γ₁⁺"=>"A₁g", "Γ₁⁻"=>"A₁ᵤ", "Γ₂⁺"=>"A₂g", "Γ₂⁻"=>"A₂ᵤ", "Γ₃⁺"=>"Eg", "Γ₃⁻"=>"Eᵤ"), # setting 2
    "6"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃"=>"²E₁", "Γ₄"=>"²E₂", "Γ₅"=>"¹E₁", "Γ₆"=>"¹E₂"),
    "-6"    => ImmutableDict("Γ₁"=>"A′", "Γ₂"=>"A′′", "Γ₃"=>"²E′", "Γ₄"=>"²E′′", "Γ₅"=>"¹E′", "Γ₆"=>"¹E′′"),
    "6/m"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"Bg", "Γ₂⁻"=>"Bᵤ", "Γ₃⁺"=>"²E₁g", "Γ₃⁻"=>"²E₁ᵤ", "Γ₄⁺"=>"²E₂g", "Γ₄⁻"=>"²E₂ᵤ", "Γ₅⁺"=>"¹E₁g", "Γ₅⁻"=>"¹E₁ᵤ", "Γ₆⁺"=>"¹E₂g", "Γ₆⁻"=>"¹E₂ᵤ"),
    "622"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"B₂", "Γ₄"=>"B₁", "Γ₅"=>"E₂", "Γ₆"=>"E₁"),
    "6mm"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"B₂", "Γ₄"=>"B₁", "Γ₅"=>"E₂", "Γ₆"=>"E₁"),
    "-62m"  => ImmutableDict("Γ₁"=>"A₁′", "Γ₂"=>"A₁′′", "Γ₃"=>"A₂′′", "Γ₄"=>"A₂′", "Γ₅"=>"E′", "Γ₆"=>"E′′"),  # setting 1
    "-6m2"  => ImmutableDict("Γ₁"=>"A₁′", "Γ₂"=>"A₁′′", "Γ₃"=>"A₂′′", "Γ₄"=>"A₂′", "Γ₅"=>"E′", "Γ₆"=>"E′′"),  # setting 2
    "6/mmm" => ImmutableDict("Γ₁⁺"=>"A₁g", "Γ₁⁻"=>"A₁ᵤ", "Γ₂⁺"=>"A₂g", "Γ₂⁻"=>"A₂ᵤ", "Γ₃⁺"=>"B₂g", "Γ₃⁻"=>"B₂ᵤ", "Γ₄⁺"=>"B₁g", "Γ₄⁻"=>"B₁ᵤ", "Γ₅⁺"=>"E₂g", "Γ₅⁻"=>"E₂ᵤ", "Γ₆⁺"=>"E₁g", "Γ₆⁻"=>"E₁ᵤ"),
    "23"    => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"¹E", "Γ₃"=>"²E", "Γ₄"=>"T"),
    "m-3"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"¹Eg", "Γ₂⁻"=>"¹Eᵤ", "Γ₃⁺"=>"²Eg", "Γ₃⁻"=>"²Eᵤ", "Γ₄⁺"=>"Tg", "Γ₄⁻"=>"Tᵤ"),
    "432"   => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E", "Γ₄"=>"T₁", "Γ₅"=>"T₂"),
    "-43m"  => ImmutableDict("Γ₁"=>"A₁", "Γ₂"=>"A₂", "Γ₃"=>"E", "Γ₄"=>"T₂", "Γ₅"=>"T₁"),
    "m-3m"  => ImmutableDict("Γ₁⁺"=>"A₁g", "Γ₁⁻"=>"A₁ᵤ", "Γ₂⁺"=>"A₂g", "Γ₂⁻"=>"A₂ᵤ", "Γ₃⁺"=>"Eg", "Γ₃⁻"=>"Eᵤ", "Γ₄⁺"=>"T₁g", "Γ₄⁻"=>"T₁ᵤ", "Γ₅⁺"=>"T₂g", "Γ₅⁻"=>"T₂ᵤ")
)

const PGIRLABS_CDML2MULLIKEN_3D_COREP = ImmutableDict(
    # Same as `PGIRLABS_CDML2MULLIKEN_3D` but with labels for physically real irreps 
    # (coreps); the label for real irreps are unchanged, but the labels for complex irreps
    # differ (e.g. ¹E and ²E becomes E). Point groups 1, -1, 2, m, 2/m, 222, mm2, mmm, 422,
    # 4mm, -42m, -4m2, 4/mmm, 312, 321, 3m1, 31m, -31m, -3m1, 622, 6mm, -62m, -6m2, 6/mmm, 
    # 432, -43m, and m-3m have only real irreps, so we don't include them here.
    "4"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃Γ₄"=>"E"),
    "-4"    => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃Γ₄"=>"E"),
    "4/m"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"Bg", "Γ₂⁻"=>"Bᵤ", "Γ₃⁺Γ₄⁺"=>"Eg", "Γ₃⁻Γ₄⁻"=>"Eᵤ"),
    "3"     => ImmutableDict("Γ₁"=>"A", "Γ₂Γ₃"=>"E"),
    "-3"    => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺Γ₃⁺"=>"Eg", "Γ₂⁻Γ₃⁻"=>"Eᵤ"),
    "6"     => ImmutableDict("Γ₁"=>"A", "Γ₂"=>"B", "Γ₃Γ₅"=>"E₁", "Γ₄Γ₆"=>"E₂"),
    "-6"    => ImmutableDict("Γ₁"=>"A′", "Γ₂"=>"A′′", "Γ₃Γ₅"=>"E′", "Γ₄Γ₆"=>"E′′"),
    "6/m"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺"=>"Bg", "Γ₂⁻"=>"Bᵤ", "Γ₃⁺Γ₅⁺"=>"E₁g", "Γ₃⁻Γ₅⁻"=>"E₁ᵤ", "Γ₄⁺Γ₆⁺"=>"E₂g", "Γ₄⁻Γ₆⁻"=>"E₂ᵤ"),
    "23"    => ImmutableDict("Γ₁"=>"A", "Γ₂Γ₃"=>"E", "Γ₄"=>"T"),
    "m-3"   => ImmutableDict("Γ₁⁺"=>"Ag", "Γ₁⁻"=>"Aᵤ", "Γ₂⁺Γ₃⁺"=>"Eg", "Γ₂⁻Γ₃⁻"=>"Eᵤ", "Γ₄⁺"=>"Tg", "Γ₄⁻"=>"Tᵤ"),
)

"""
$(TYPEDSIGNATURES)

Return the Mulliken label of a point group irrep `pgir`.

## Notes
This functionality is a simple mapping between the tabulated CDML point group irrep labels
and associated Mulliken labels, using the listings at the Bilbao Crystallographic
Database [^1].

Ignoring subscript, the rough rules associated with assignment of Mulliken labels are:

1. **Irrep dimensionality**: 
    - **1D irreps**: if a real irrep, assign A or B (B if antisymmetric under a principal 
      rotation); if a complex irrep, assigned label ¹E or ²E.
    - **2D irreps**: assign label E
    - **3D irreps**: assign label T
2. **_u_ and _g_ subscripts**: if the group contains inversion, indicate whether irrep is
   symmetric (g ~ gerade) or antisymmetric (ᵤ ~ ungerade) under inversion.
3. **Prime superscripts**: if the group contains a mirror *m* aligned with a principal 
   rotation axis, but does *not* contain inversion, indicate whether irrep is symmetric (′) 
   or antisymmetric (′′) under this mirror.
4. **Numeral subscripts**: the rules for assignment of numeral subscripts are too
   too complicated in general - and indeed, we are unaware of a general coherent rule -- to
   describe here.

## References
[^1]: Bilbao Crystallographic Database's
      [Representations PG program](https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=spg).
"""
function mulliken(pgir::PGIrrep{D}) where D
    pglab   = label(group(pgir))
    pgirlab = label(pgir)
    return _mulliken(pglab, pgirlab, iscorep(pgir))
end
function _mulliken(pglab, pgirlab, iscorep) # split up to let `SiteIrrep` overload `mulliken`
    if iscorep
        return PGIRLABS_CDML2MULLIKEN_3D_COREP[pglab][pgirlab]
    else
        return PGIRLABS_CDML2MULLIKEN_3D[pglab][pgirlab]
    end
end

#=
# ATTEMPT AT ACTUALLY IMPLEMENTING THE MULLIKEN NOTATION OURSELVES, USING THE "RULES"
# UNFORTUNATELY, THERE IS A LACK OF SPECIFICATION OF THESE RULES; THIS E.G. IMPACTS:
#       - ¹ and ² superscripts to E labels: no rules, whatsoever. Cannot reverse-engineer
#         whatever the "rule" is (if there is any) that fits all point groups; e.g., rule
#         seems different between e.g. {4, -4, 4/m} and {6, -6, 23, m-3}. the rule we went
#         with below works for {4, -4, 4/m}, but not {6, -6, 23, m-3}
#       - how to infer that a ₁₂₃ subscript to a label is not needed (i.e. that the label
#         is already unambiguous? there just doesn't seem to be any way to infer this
#         generally, without looking at all the different irreps at the same time.
#       - straaange corner cases for A/B label assignment; e.g. -6m2 is all A labels,
#         but the principal rotation (-6) has characters with both positive and negative 
#         sign; so the only way strictly A-type labels could be assigned if we pick some
#         other principal rotation, e.g. 3₀₀₁... that doesn't make sense.
#       - sometimes, subscript assignment differs, e.g. for 622: B₁ vs. B₂. the rule used
#         to pick subscripts are ambiguous here, because there are two sets of two-fold
#         rotation operations perpendicular to the principal axis - and they have opposite
#         signs; so the assignment of ₁₂ subscripts depends on arbitrarily picking one of
#         these sets.
#       - more issues probably exist...
# BECAUSE OF THIS, WE JUST OPT TO ULTIMATELY JUST NOT IMPLEMENT THIS OURSELVES, AND INSTEAD
# RESTRICT OURSELVES TO GETTING THE MULLIKEN NOTATION ONLY FOR TABULATED POINT GROUPS BY
# COMPARING WITH THE ASSOCIATED LABELS
#
# TENTATIVE IMPLEMENTATION:
# Rough guidelines are e.g. in http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Mulliken.html &
# xuv.scs.illinois.edu/516/handouts/Mulliken%20Symbols%20for%20Irreducible%20Representations.pdf
# The "canonical" standard/most explicit resource is probably Pure & Appl. Chem., 69, 1641
# (1997), see e.g. https://core.ac.uk/download/pdf/9423.pdf or Mulliken's own 'Report on
# Notation for the Spectra of Polyatomic Molecules' (https://doi.org/10.1063/1.1740655)
#
# Note that the convention's choices not terribly specific when it comes to cornercases,
# such as for B-type labels with 3 indices or E-type labels with 1/2 indices
function mulliken(ir::Union{PGIrrep{D}, LGIrrep{D}}) where D
    D ∈ (1,2,3) || _throw_invaliddim(D)
    if ir isa LGIrrep && !issymmorphic(num(ir), D)
        error("notation not defined for `LGIrrep`s of nonsymmorphic space groups")
    end
       
    g   = group(ir)
    χs  = characters(ir)
    irD = irdim(ir)

    # --- determine "main" label of the irrep (A, B, E, or T) ---
    # find all "maximal" (=principal) rotations, including improper ones
    idxs_principal = find_principal_rotations(g, include_improper=true)
    op_principal   = g[first(idxs_principal)]    # a representative principal rotation
    axis_principal = D == 1              ? nothing : # associcated rotation axis
                     isone(op_principal) ? nothing :
                     D == 2              ? SVector(0,0,1) :
                  #= D == 3 =#             rotation_axis_3d(op_principal)
    lab = if irD == 1     # A or B
        χ′s = @view χs[idxs_principal]
        # if character of rotation around _any_ principal rotation axis (i.e. maximum 
        # rotation order) is negative, then label is B; otherwise A. If complex character
        # label is ¹E or ²E
        if all(x->isapprox(abs(real(x)), 1, atol=DEFAULT_ATOL), χ′s) # => real irrep 
            any(isapprox(-1, atol=DEFAULT_ATOL), χ′s) ? "B" : "A"

        else                                                         # => complex irrep
            # must then be a complex rep; label is ¹E or ²E. the convention at e.g. Bilbao
            # seems to be to (A) pick a principal rotation with positive sense ("⁺"),  
            # then (B) look at its character, χ, then (C) if Imχ < 0 => assign ¹
            # superscript, if Imχ > 0 => assign ² superscript. this is obviously a very
            # arbitrary convention, but I guess it's as good as any
            idx⁺ = findfirst(idx->seitz(g[idx])[end] == '⁺', idxs_principal)
            idx⁺ === nothing && (idx⁺ = 1) # in case of 2-fold rotations where 2⁺ = 2⁻

            imag(χ′s[idx⁺]) < 0 ? "¹E" : "²E"
        end

    elseif irD == 2 # E (reality is always real for PGIrreps w/ irD = 3)
        "E"

    elseif irD == 3 # T (reality is always real for PGIrreps w/ irD = 3)
        "T"

    else
        throw(DomainError(irD, "the dimensions of a crystallographic point group irrep "*
                               "is expected to be between 1 and 3"))
        # in principle, 4 => "G" and 5 => "H", but irreps of dimension larger than 3 never
        # arise for crystallographic point groups; not even when considering time-reversal
    end

    # --> number subscript ₁, ₂, ₃
    # NOTE: these rules are very messy and also not well-documented; it is especially
    #       bad for E and T labels; there, I've just inferred a rule that works for the
    #       crystallographic point groups, simply by comparing to published tables (e.g.,
    #       Bilbao and http://symmetry.jacobs-university.de/)
    if axis_principal !== nothing
        if lab == "A" || lab == "B"
            # ordinarily, we can follow a "simple" scheme - but there is an arbitrary 
            # choice involved when we have point groups 222 or mmm - where we need to assign
            # 3 different labels to B-type irreps

            # special rules needed for "222" (C₂) and "mmm" (D₂ₕ) groups
            # check B-type label and point group 222 [identity + 3×(2-fold rotation)] or
            # mmm [identity + 3×(2-fold rotation + mirror) + inversion]
            istricky_B_case = if lab == "B"
                (length(g)==4 && count(op->rotation_order(op)==2, g) == 3)      || # 222
                (length(g)==8 && count(op->abs(rotation_order(op))==2, g) == 6)    # mmm
            else
                false
            end

            if !istricky_B_case
                # simple scheme: 
                # find a 2-fold rotation or mirror (h) whose axis is ⟂ to principal axis:
                #   χ(h) = +1  =>  ₁
                #   χ(h) = -1  =>  ₂
                idxᴬᴮ = if D == 3 
                    findfirst(g) do op
                        abs(rotation_order(op)) == 2 && 
                        dot(rotation_axis_3d(op), axis_principal) == 0
                    end

                elseif D == 2
                    findfirst(op -> rotation_order(op) == -2, g)

                end
                if idxᴬᴮ !== nothing
                    lab *= real(χs[idxᴬᴮ]) > 0 ? '₁' : '₂'
                end
            
            else # tricky B case (222 or mmm)
                # there's no way to do this independently of a choice of setting; we just 
                # follow what seems to be the norm and _assume_ the presence of 2₀₀₁, 2₀₁₀,
                # and 2₁₀₀; if they don't exist, things will go bad! ... no way around it
                idxˣ = findfirst(op->seitz(op)=="2₁₀₀", g)
                idxʸ = findfirst(op->seitz(op)=="2₀₁₀", g)
                idxᶻ = findfirst(op->seitz(op)=="2₀₀₁", g)
                if idxˣ === nothing || idxʸ === nothing || idxᶻ === nothing
                    error("cannot assign tricky B-type labels for nonconventional axis settings")
                end
                χsˣʸᶻ = (χs[idxˣ], χs[idxʸ], χs[idxᶻ])
                lab *= χsˣʸᶻ == (-1,-1,+1) ? '₁' :
                       χsˣʸᶻ == (-1,+1,-1) ? '₂' :
                       χsˣʸᶻ == (+1,-1,-1) ? '₃' : error("unexpected character combination")
            end


        elseif lab == "E" || lab == "¹E" || lab == "²E"
            # TODO
            # disambiguation needed for pgs 6 (C₆), 6/m (C₆ₕ), 622 (D₆), 6mm (C₆ᵥ), 
            # 6/mmm (D₆ₕ)
            # rule seems to be that we look at the 2-fold rotation operation aligned with
            # the principal axis; denoting this operation by 2, the rule is:
            #   χ(2) = -1  =>  E₁
            #   χ(2) = +1  =>  E₂ 
            idxᴱ = findfirst(g) do op
                rotation_order(op) == 2 && (D==2 || rotation_axis_3d(op) == axis_principal)
            end
            if idxᴱ !== nothing
                lab *= real(χs[idxᴱ]) < 0 ? '₁' : '₂'
            end

        elseif lab == "T"
            # disambiguation of T symbols is only needed for 432 (O), -43m (Td) and m-3m
            # (Oₕ) pgs; the disambiguation can be done by checking the sign of the character
            # of a principal rotation; which in this case is always a 4-fold proper or 
            # improper rotation ±4, the rule thus is:
            #   χ(±4) = -1  =>  T₂
            #   χ(±4) = +1  =>  T₁
            # to avoid also adding unnecessary subscripts to pgs 23 and m-3 (whose T-type
            # irreps are already disambiguated), we check if the principal rotation has
            # order 4 - if not, this automatically excludes 23 and m-3.

            if abs(rotation_order(op_principal)) == 4
                # we only need to check a single character; all ±4 rotations have the same
                # character in this cornercase
                lab *= real(χs[first(idxs_principal)]) > 0 ? '₁' : '₂'
            end
        end
    end

    # --> letter subscript g, ᵤ and prime superscript ′, ′′
    idx_inversion = findfirst(op -> isone(-rotation(op)), g)
    if idx_inversion !== nothing
        # --> letter subscript g, ᵤ
        # rule applies for groups with inversion -1:
        #   χ(-1) = +1  =>  _g [gerade ~ even]
        #   χ(-1) = -1  =>  ᵤ  [ungerade ~ odd]
        lab *= real(χs[idx_inversion]) > 0 ? 'g' : 'ᵤ'

    elseif length(g) > 1
        # --> prime superscript ′, ′′
        # rule applies for groups without inversion and with a mirror aligned with a 
        # principal rotation axis; most often m₀₀₁ in 3D _or_ for the 2-element group
        # consisting only of {1, "mᵢⱼₖ"}. Denoting the mirror by "m", the rule is:
        #   χ(m) = +1  =>  ′
        #   χ(m) = -1  =>  ′′
        # rule is relevant only for point groups m, -6, -62m in 3D and m in 2D
        idxᵐ = if length(g) == 2 
            findfirst(op -> occursin("m", seitz(op)), g) # check if {1, "mᵢⱼₖ"} case
            
        elseif D == 3
            # mirror aligned with principal axis casecan now assume that D == 3
            # (only occurs if D = 3; so if D = 2, we just move on)
            findfirst(g) do op # find aligned mirror
                (rotation_order(op) == -2) && (rotation_axis_3d(op) == axis_principal)
            end
        end

        if idxᵐ !== nothing
            lab *= real(χs[idxᵐ]) > 0 ? "′" : "′′"
        end
    end

    return lab
end
=#

"""
$(TYPEDSIGNATURES)

Return the the indices of the "maximal" rotations among a set of operations `ops`, i.e. 
those of maximal order (the "principal rotations").

## Keyword arguments
- `include_improper` (`=false`): if `true`, improper rotations are included in the
  consideration. If the order of improper and proper rotations is identical, only 
  the indices of the proper rotations are returned. If the maximal (signed) rotation
  order is -2 (a mirror), it is ignored, and the index of the identity operation is
  returned.
"""
function find_principal_rotations(ops::AbstractVector{SymOperation{D}}; 
                                  include_improper::Bool=false) where D
    # this doesn't distinguish between rotation direction; i.e., picks either ⁺ or ⁻ 
    # depending on ordering of `ops`. if there are no rotations, then we expect that there
    # will at least always be an identity operation - otherwise, we throw an error
    # may pick either a proper or improper rotation; picks a proper rotation if max order
    # of proper rotation exceeds order of max improper rotation
    rots = rotation_order.(ops)
    maxrot = maximum(rots)     # look for proper rotations
    rot = if include_improper
        minrot = minimum(rots) # look for improper rotations
        maxrot ≥ abs(minrot) ? maxrot : (minrot == -2 ? maxrot : minrot)
    else
        maxrot
    end
    idxs = findall(==(rot), rots)
    return idxs::Vector{Int}
end