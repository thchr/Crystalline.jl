"""
    schoenflies(sgnum::Integer) --> String

Returns the Schoenflies notation for a given space group number
`sgnum`. Schoenflies notation only applies to point groups and 
space groups, not plane groups, so this notation is only relevant
in three dimensions.
"""
schoenflies(sgnum::Integer) = SCHOENFLIES_TABLE[sgnum]

"""
    hermannmauguin(sgnum::Integer, dim::Integer=3) --> String

Returns the Hermann-Mauguin notation for a given space group number
`sgnum` and dimensionality `dim`, sometimes also called the IUC 
(International Union of Crystallography) or international notation 
(since it is used in the International Tables of Crystallography); 
accordingly, the functionality is aliased by `iuc(sgnum, dim)`. 
Hermann-Mauguin notation applies in two and three-dimensions.

For additional information see https://en.wikipedia.org/wiki/Hermann%E2%80%93Mauguin_notation.
"""
hermannmauguin(sgnum::Integer, dim::Integer=3) = HERMANNMAUGUIN_TABLE[dim][sgnum]
const iuc = hermannmauguin # alias

""" 
    centering(sgnum::Integer, dim::Integer=3) --> Char

Determines the conventional centering type of a given space/plane group number
`sgnum` by comparison with the Hermann-Mauguin notation's first letter. 

Possible output values, depending on dimensionality `dim`, are (see ITA Sec. 9.1.4):

    dim=2 ┌ 'p': no centring (primitive)
          └ 'c': face centered

    dim=3 ┌ 'P': no centring (primitive)
          ├ 'I': body centred (innenzentriert)
          ├ 'F': all-face centred
          ├ 'A', 'B', 'C': one-face centred, (b,c) or (c,a) or (a,b)
          └ 'R': hexagonal cell rhombohedrally centred
"""
centering(sgnum::Integer, dim::Integer=3) = first(hermannmauguin(sgnum, dim))

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

# Hermann-Mauguin notation, ordered relative to space/plane group number
const HERMANNMAUGUIN_TABLE = (
# ------------------------------------------------------------------------------------------
# line-group notation (one dimension) [see https://en.wikipedia.org/wiki/Line_group]
# ------------------------------------------------------------------------------------------    
("p1", "p1m"),
# ------------------------------------------------------------------------------------------
# plane-group notation (two dimensions) [see e.g. Table 19 of Cracknell, Adv. Phys. 1974]
# ------------------------------------------------------------------------------------------
(
# oblique
"p1",   "p211",
# rectangular ('p' or 'c' centering; c-centered lattices are rhombic in their primitive cell)
"p1m1", "p1g1", "c1m1", "p2mm", "p2mg", "p2gg", "c2mm",   
# square
"p4",   "p4mm", "p4gm",
# hexagonal
"p3",   "p3m1", "p31m", "p6",   "p6mm"
),
# ------------------------------------------------------------------------------------------
# space-group notation (three dimensions) [adapted from https://bruceravel.github.io/demeter/artug/atoms/space.html]
# ------------------------------------------------------------------------------------------
(
# triclinic
"P1",      "P-1",
# monoclinic
"P2",      "P21",     "C2",      "PM",      "PC",      "CM",       
"CC",      "P2/M",    "P21/M",   "C2/M",    "P2/C",    "P21/C",    
"C2/C",
# orthorhombic
"P222",    "P2221",   "P21212",  "P212121", "C2221",   "C222",
"F222",    "I222",    "I212121", "PMM2",    "PMC21",   "PCC2",
"PMA2",    "PCA21",   "PNC2",    "PMN21",   "PBA2",    "PNA21",
"PNN2",    "CMM2",    "CMC21",   "CCC2",    "AMM2",    "ABM2",
"AMA2",    "ABA2",    "FMM2",    "FDD2",    "IMM2",    "IBA2",
"IMA2",    "PMMM",    "PNNN",    "PCCM",    "PBAN",    "PMMA",
"PNNA",    "PMNA",    "PCCA",    "PBAM",    "PCCN",    "PBCM",
"PNNM",    "PMMN",    "PBCN",    "PBCA",    "PNMA",    "CMCM",
"CMCA",    "CMMM",    "CCCM",    "CMMA",    "CCCA",    "FMMM",
"FDDD",    "IMMM",    "IBAM",    "IBCA",    "IMMA",
# tetragonal
"P4",      "P41",     "P42",     "P43",     "I4",      "I41",
"P-4",     "I-4",     "P4/M",    "P42/M",   "P4/N",    "P42/N",
"I4/M",    "I41/A",   "P422",    "P4212",   "P4122",   "P41212",
"P4222",   "P42212",  "P4322",   "P43212",  "I422",    "I4122",
"P4MM",    "P4BM",    "P42CM",   "P42NM",   "P4CC",    "P4NC",
"P42MC",   "P42BC",   "I4MM",    "I4CM",    "I41MD",   "I41CD",
"P-42M",   "P-42C",   "P-421M",  "P-421C",  "P-4M2",   "P-4C2",
"P-4B2",   "P-4N2",   "I-4M2",   "I-4C2",   "I-42M",   "I-42D",
"P4/MMM",  "P4/MCC",  "P4/NBM",  "P4/NNC",  "P4/MBM",  "P4/MNC",
"P4/NMM",  "P4/NCC",  "P42/MMC", "P42/MCM", "P42/NBC", "P42/NNM",
"P42/MBC", "P42/MNM", "P42/NMC", "P42/NCM", "I4/MMM",  "I4/MCM",
"I41/AMD", "I41/ACD",
# trigonal
"P3",      "P31",     "P32",     "R3",      "P-3",     "R-3",
"P312",    "P321",    "P3112",   "P3121",   "P3212",   "P3221",
"R32",     "P3M1",    "P31M",    "P3C1",    "P31C",    "R3M",
"R3C",     "P-31M",   "P-31C",   "P-3M1",   "P-3C1",   "R-3M",
"R-3C",
# hexagonal
"P6",      "P61",     "P65",     "P62",     "P64",     "P63",
"P-6",     "P6/M",    "P63/M",   "P622",    "P6122",   "P6522",
"P6222",   "P6422",   "P6322",   "P6MM",    "P6CC",    "P63CM",
"P63MC",   "P-6M2",   "P-6C2",   "P-62M",   "P-62C",   "P6/MMM",
"P6/MCC",  "P63/MCM", "P63/MMC",
# cubic
"P23",     "F23",     "I23",     "P213",    "I213",    "PM3",
"PN3",     "FM3",     "FD3",     "IM3",     "PA3",     "IA3",
"P432",    "P4232",   "F432",    "F4132",   "I432",    "P4332",
"P4132",   "I4132",   "P-43M",   "F-43M",   "I-43M",   "P-43N",
"F-43C",   "I-43D",   "PM3M",    "PN3N",    "PM3N",    "PN3M",
"FM3M",    "FM3C",    "FD3M",    "FD3C",    "IM3M",    "IA3D"
)
)



""" 
    seitz(op::SymOperation) --> String

Computes the correponding Seitz notation {β|τ} for a symmetry operation in 
triplet form.

Implementation based on ITA5 Table 11.2.1.1 (for 3D)\n
        ________________________________________________
        |_detW_|_trW_|_-3_|_-2 |_-1 |__0_|__1_|__2_|__3_|
        |  1         |    |    |  2 |  3 |  4 |  6 |  1 |
        |__1_________|_-1_|_-6_|_-4_|_-3_|__m_|____|____|
with the elements of the table giving the type of symmetry operation in
in Hermann-Mauguin notation. The rotation axis and the rotation sense are 
computed following the rules in ITA6 Sec. 1.2.2.4(1)(b-c).
The implementation has been checked against the Tables 1.4.2.1-5 of ITA6.

Note that the orientation of axis (i.e. its sign) is not necessarily equal
to the orientation picked in those tables; it is a matter of convention,
and the conventions have not been explicated in ITA6.

For 2D operations, we elevate the operation to one in 3D that leaves the 
3rd coordinate invariant, and then compute results using the 3D procedure.
"""
function seitz(op::SymOperation)
    W = rotation(op); w = translation(op); dim = size(W,1)
    if dim == 2 # we just augment the 2D case by leaving z invariant
        W = [W zeros(2); 0.0 0.0 1.0]; 
        w = [w; 0]
    elseif dim == 1
        w_str = !iszero(w[1]) ? unicode_frac(w[1]) : "0"
        if isone(W[1])
            return "{1|"*w_str*"}"
        elseif isone(-W[1])
            return "{-1|"*w_str*"}"
        else
            throw(DomainError((W,w), "not a valid 1D symmetry operation"))
        end
    end

    detW = det(W); detW′, detW = detW, round(Int64, detW) # det, then round & flip
    isapprox(detW′, detW, atol=DEFAULT_ATOL) || throw(ArgumentError("det W must be an integer for a SymOperation {W|w}; got $(detW′)"))
    trW  = tr(W);  trW′,  trW  = trW, round(Int64, trW)   # tr, then round & flip
    isapprox(trW′, trW, atol=DEFAULT_ATOL) || throw(ArgumentError("tr W must be an integer for a SymOperation {W|w}; got $(trW′)"))

    # --- rotation order (and proper/improper determination) ---
    if detW == 1 # proper rotations
        if -1 ≤ trW ≤ 1 # 2-, 3-, or 4-fold rotation
            rot = trW + 3
        elseif trW == 2 # 6-fold rotation
            rot = 6
        elseif trW == 3 # identity operation
            rot = 1
        else 
            _throw_seitzerror(trW, detW)
        end
    elseif detW == -1 # improper rotations (rotoinversions)
        if trW == -3    # inversion
            rot = -1
        elseif trW == -2 # 6-fold rotoinversion
            rot= -6
        elseif -1 ≤ trW ≤ 0 # 4- and 3-fold rotoinversion
            rot = trW - 3
        elseif trW == 1  # mirror, note that "m" == "-2" conceptually
            rot = -2
        else
            _throw_seitzerror(trW, detW)
        end
    else
        _throw_seitzerror(trW, detW)
    end
    order = abs(rot)
    rot_str = rot == -2 ? "m" : string(rot)
    
    # --- rotation axis (for order ≠ 1)---
    # the rotation axis 𝐮 is determined from the product of
    # 𝐘ₖ(𝐖) ≡ (d𝐖)ᵏ⁻¹+(d𝐖)ᵏ⁻² + ... + (d𝐖) + 𝐈 where d ≡ det(𝐖) 
    # with an arbitrary vector 𝐯 that is not perpendicular to 𝐮
    # [cf. ITA6  Vol. A, p. 16, Sec. 1.2.2.4(1)(b)]
    if dim == 3 && order == 1 || dim == 2 && rot ≠ -2 # only need orientation in 2D for mirrors 
        axis_str = ""                                 # (w/ in plane normals; otherwise along [001])
        u = dim == 2 ? [0, 0, 1] : [0, 0, 0]
    else
        Yₖ = Matrix{Float64}(I, 3, 3) # calculate Yₖ by iteration
        for j=1:order-1
            term = W^j
            if detW^j == -1;
                Yₖ .-= term 
            else
                Yₖ .+= term
            end
        end
        u = zeros(Float64, 3)
        while iszero(u)
            v = rand(3); 
            u = Yₖ*v # there is near-infinitesimal chance that u is zero for random v, but we check anyway.
        end
        norm = minimum(Base.Filter(x->x>DEFAULT_ATOL,abs.(u))) # minimum nonzero element
        u ./= norm # normalize
        u′,  u  = u, round.(Int64, u) # convert from float to integer and check validity of conversion
        isapprox(u′, u, atol=DEFAULT_ATOL) || throw(ArgumentError("the rotation axis must be equivalent to an integer vector by appropriate normalization; got $(u′)"))
        # the sign of u is arbitrary: we adopt the convention of '-' elements
        # coming "before" '+' elements; e.g. [-1 -1 1] is picked over [1 1 -1]
        # and [-1 1 -1] is picked over [1 -1 1]; note that this impacts the 
        # sense of rotation which depends on the sign of the rotation axis
        negidx = findfirst(x->x<0, u)
        firstnonzero = findfirst(x-> abs(x) ≠ 0, u)
        if negidx ≠ nothing && (negidx ≠ firstnonzero || negidx === firstnonzero === 3); u .*= -1; end

        axis_str = subscriptify(join(string(u[i]) for i in 1:dim)) # for 2D, ignore z-component
    end
    
    # --- rotation sense (for order > 2}) ---
    # ±-rotation sense is determined from sign of det(𝐙) where
    # 𝐙 ≡ [𝐮|𝐱|det(𝐖)𝐖𝐱] where 𝐱 is an arbitrary vector that 
    # is not parallel to 𝐮. [ITA6  Vol. A, p. 16, Sec. 1.2.2.4(1)(c)]
    if order > 2
        while true
            global x = rand(Int64, 3)
            iszero(x×u) || break # check that generated 𝐱 is not parallel to 𝐮 (if it is, 𝐱×𝐮 = 0)
        end
        Z = [u x (detW*W*x)]
        sense_str = signbit(det(Z)) ? "⁻" : "⁺"
    else
        sense_str = ""
    end

    # --- nonsymmorphic part ---
    w_str = !iszero(w) ? join((unicode_frac(w[i]) for i in 1:dim), ',') : "0"
        

    # --- combine labels ---
    return '{' * rot_str * sense_str * axis_str * '|' * w_str * '}'
end
seitz(str::String) = seitz(SymOperation(str))
_throw_seitzerror(trW, detW) = throw(ArgumentError("trW = $(trW) for detW = $(detW) is not a valid symmetry operation; see ITA5 Vol A, Table 11.2.1.1"))