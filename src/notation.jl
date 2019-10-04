"""
    schoenflies(sgnum::Integer) --> String

Returns the Schoenflies notation for a given space group number
`sgnum`. Schoenflies notation only applies to point groups and 
space groups, not plane groups, so this notation is only relevant
in three dimensions.
"""
function schoenflies(sgnum::Integer) 
    return schoenflies_table[sgnum]
end

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
function hermannmauguin(sgnum::Integer, dim::Integer=3) 
    return hermannmauguin_table[dim][sgnum]
end
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
const schoenflies_table = [
# triclinic
"C_1^1",    "C_i^1", 
# monoclinic
"C_2^1",    "C_2^2",    "C_2^3",    "C_S^1",    "C_S^2",    "C_S^3", 
"C_S^4",    "C_2h^1",   "C_2h^2",   "C_2h^3",   "C_2h^4",   "C_2h^5", 
"C_2h^6",
# orthorhombic
"D_2^1",    "D_2^2",    "D_2^3",    "D_2^4",    "D_2^5",    "D_2^6", 
"D_2^7",    "D_2^8",    "D_2^9",    "C_2v^1",   "C_2v^2",   "C_2v^3", 
"C_2v^4",   "C_2v^5",   "C_2v^6",   "C_2v^7",   "C_2v^8",   "C_2v^9", 
"C_2v^10",  "C_2v^11",  "C_2v^12",  "C_2v^13",  "C_2v^14",  "C_2v^15", 
"C_2v^16",  "C_2v^17",  "C_2v^18",  "C_2v^19",  "C_2v^20",  "C_2v^21", 
"C_2v^22",  "D_2h^1",   "D_2h^2",   "D_2h^3",   "D_2h^4",   "D_2h^5", 
"D_2h^6",   "D_2h^7",   "D_2h^8",   "D_2h^9",   "D_2h^10",  "D_2h^11", 
"D_2h^12",  "D_2h^13",  "D_2h^14",  "D_2h^15",  "D_2h^16",  "D_2h^17", 
"D_2h^18",  "D_2h^19",  "D_2h^20",  "D_2h^21",  "D_2h^22",  "D_2h^23", 
"D_2h^24",  "D_2h^25",  "D_2h^26",  "D_2h^27",  "D_2h^28",
# tetragonal 
"C_4^1",    "C_4^2",    "C_4^3",    "C_4^4",    "C_4^5",    "C_4^6", 
"S_4^1",    "S_4^2",    "C_4h^1",   "C_4h^2",   "C_4h^3",   "C_4h^4", 
"C_4h^5",   "C_4h^6",   "D_4^1",    "D_4^2",    "D_4^3",    "D_4^4", 
"D_4^5",    "D_4^6",    "D_4^7",    "D_4^8",    "D_4^9",    "D_4^10", 
"C_4v^1",   "C_4v^2",   "C_4v^3",   "C_4v^4",   "C_4v^5",   "C_4v^6", 
"C_4v^7",   "C_4v^8",   "C_4v^9",   "C_4v^10",  "C_4v^11",  "C_4v^12", 
"D_2d^1",   "D_2d^2",   "D_2d^3",   "D_2d^4",   "D_2d^5",   "D_2d^6", 
"D_2d^7",   "D_2d^8",   "D_2d^9",   "D_2d^10",  "D_2d^11",  "D_2d^12", 
"D_4h^1",   "D_4h^2",   "D_4h^3",   "D_4h^4",   "D_4h^5",   "D_4h^6", 
"D_4h^7",   "D_4h^8",   "D_4h^9",   "D_4h^10",  "D_4h^11",  "D_4h^12", 
"D_4h^13",  "D_4h^14",  "D_4h^15",  "D_4h^16",  "D_4h^17",  "D_4h^18", 
"D_4h^19",  "D_4h^20",
# trigonal
"C_3^1",    "C_3^2",    "C_3^3",    "C_3^4",    "C_3i^1",   "C_3i^2", 
"D_3^1",    "D_3^2",    "D_3^3",    "D_3^4",    "D_3^5",    "D_3^6", 
"D_3^7",    "C_3v^1",   "C_3v^2",   "C_3v^3",   "C_3v^4",   "C_3v^5", 
"C_3v^6",   "D_3d^1",   "D_3d^2",   "D_3d^3",   "D_3d^4",   "D_3d^5", 
"D_3d^6",
# hexagonal
"C_6^1",    "C_6^2",    "C_6^3",    "C_6^4",    "C_6^5",    "C_6^6", 
"C_3h^1",   "C_6h^1",   "C_6h^2",   "D_6^1",    "D_6^2",    "D_6^3", 
"D_6^4",    "D_6^5",    "D_6^6",    "C_6v^1",   "C_6v^2",   "C_6v^3", 
"C_6v^4",   "D_3h^1",   "D_3h^2",   "D_3h^3",   "D_3h^4",   "D_6h^1", 
"D_6h^2",   "D_6h^3",   "D_6h^4",
# cubic
"T^1",      "T^2",      "T^3",      "T^4",      "T^5",      "T_h^1", 
"T_h^2",    "T_h^3",    "T_h^4",    "T_h^5",    "T_h^6",    "T_h^7", 
"O^1",      "O^2",      "O^3",      "O^4",      "O^5",      "O^6", 
"O^7",      "O^8",      "T_d^1",    "T_d^2",    "T_d^3",    "T_d^4", 
"T_d^5",    "T_d^6",    "O_h^1",    "O_h^2",    "O_h^3",    "O_h^4", 
"O_h^5",    "O_h^6",    "O_h^7",    "O_h^8",    "O_h^9",    "O_h^10"
]

# Hermann-Mauguin notation, ordered relative to space/plane group number
const hermannmauguin_table = Dict{Int64, Vector{String}}(
3 => # space-group notation (three dimensions) [from https://bruceravel.github.io/demeter/artug/atoms/space.html]
[
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
]
, 
2 => # plane group notation (two dimensions) [see e.g. Table 19 of Cracknell, Adv. Phys. 1974]
[
# oblique
"p1",   "p211",
# rectangular ('p' or 'c' centering; c-centered lattices are rhombic in their primitive cell)
"p1m1", "p1g1", "c1m1", "p2mm", "p2mg", "p2gg", "c2mm",   
# square
"p4",   "p4mm", "p4gm",
# hexagonal
"p3",   "p3m1", "p31m", "p6",   "p6mm"
]
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
    end

    detW = det(W); detW′, detW = detW, round(Int64, detW) # det, then round & flip
    detW′ ≈ detW || throw(ArgumentError("det W must be an integer for a SymOperation {W|w}"))
    trW  = tr(W);  trW′,  trW  = trW, round(Int64, trW)   # tr, then round & flip
    trW′ ≈ trW || throw(ArgumentError("tr W must be an integer for a SymOperation {W|w}"))

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
        norm = minimum(Base.Filter(!iszero,u)) # minimum nonzero element
        u ./= norm # normalize
        u′,  u  = u, round.(Int64, u) # convert from float to integer and check validity of conversion
        u′ ≈ u || throw(ArgumentError("the rotation axis must be equivalent to an integer vector by appropriate normalization; got $(u′)"))
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