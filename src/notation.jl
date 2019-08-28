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
"C_1^1",      "C_i^1", 
# monoclinic
"C_2^1",      "C_2^2",      "C_2^3",      "C_S^1",      "C_S^2",       "C_S^3", 
"C_S^4",      "C_2h^1",     "C_2h^2",     "C_2h^3",     "C_2h^4",      "C_2h^5", 
"C_2h^6",
# orthorhombic
"D_2^1",      "D_2^2",      "D_2^3",      "D_2^4",      "D_2^5",       "D_2^6", 
"D_2^7",      "D_2^8",      "D_2^9",      "C_2v^1",     "C_2v^2",      "C_2v^3", 
"C_2v^4",     "C_2v^5",     "C_2v^6",     "C_2v^7",     "C_2v^8",      "C_2v^9", 
"C_2v^10",    "C_2v^11",    "C_2v^12",    "C_2v^13",    "C_2v^14",     "C_2v^15", 
"C_2v^16",    "C_2v^17",    "C_2v^18",    "C_2v^19",    "C_2v^20",     "C_2v^21", 
"C_2v^22",    "D_2h^1",     "D_2h^2",     "D_2h^3",     "D_2h^4",      "D_2h^5", 
"D_2h^6",     "D_2h^7",     "D_2h^8",     "D_2h^9",     "D_2h^10",     "D_2h^11", 
"D_2h^12",    "D_2h^13",    "D_2h^14",    "D_2h^15",    "D_2h^16",     "D_2h^17", 
"D_2h^18",    "D_2h^19",    "D_2h^20",    "D_2h^21",    "D_2h^22",     "D_2h^23", 
"D_2h^24",    "D_2h^25",    "D_2h^26",    "D_2h^27",    "D_2h^28",
# tetragonal 
"C_4^1",      "C_4^2",      "C_4^3",      "C_4^4",      "C_4^5",       "C_4^6", 
"S_4^1",      "S_4^2",      "C_4h^1",     "C_4h^2",     "C_4h^3",      "C_4h^4", 
"C_4h^5",     "C_4h^6",     "D_4^1",      "D_4^2",      "D_4^3",       "D_4^4", 
"D_4^5",      "D_4^6",      "D_4^7",      "D_4^8",      "D_4^9",       "D_4^10", 
"C_4v^1",     "C_4v^2",     "C_4v^3",     "C_4v^4",     "C_4v^5",      "C_4v^6", 
"C_4v^7",     "C_4v^8",     "C_4v^9",     "C_4v^10",    "C_4v^11",     "C_4v^12", 
"D_2d^1",     "D_2d^2",     "D_2d^3",     "D_2d^4",     "D_2d^5",      "D_2d^6", 
"D_2d^7",     "D_2d^8",     "D_2d^9",     "D_2d^10",    "D_2d^11",     "D_2d^12", 
"D_4h^1",     "D_4h^2",     "D_4h^3",     "D_4h^4",     "D_4h^5",      "D_4h^6", 
"D_4h^7",     "D_4h^8",     "D_4h^9",     "D_4h^10",    "D_4h^11",     "D_4h^12", 
"D_4h^13",    "D_4h^14",    "D_4h^15",    "D_4h^16",    "D_4h^17",     "D_4h^18", 
"D_4h^19",    "D_4h^20",
# trigonal
"C_3^1",       "C_3^2",     "C_3^3",      "C_3^4",      "C_3i^1",      "C_3i^2", 
"D_3^1",       "D_3^2",     "D_3^3",      "D_3^4",      "D_3^5",       "D_3^6", 
"D_3^7",       "C_3v^1",    "C_3v^2",     "C_3v^3",     "C_3v^4",      "C_3v^5", 
"C_3v^6",      "D_3d^1",    "D_3d^2",     "D_3d^3",     "D_3d^4",      "D_3d^5", 
"D_3d^6",
# hexagonal
"C_6^1",       "C_6^2",     "C_6^3",       "C_6^4",     "C_6^5",       "C_6^6", 
"C_3h^1",      "C_6h^1",    "C_6h^2",      "D_6^1",     "D_6^2",       "D_6^3", 
"D_6^4",       "D_6^5",     "D_6^6",       "C_6v^1",    "C_6v^2",      "C_6v^3", 
"C_6v^4",      "D_3h^1",    "D_3h^2",      "D_3h^3",    "D_3h^4",      "D_6h^1", 
"D_6h^2",      "D_6h^3",    "D_6h^4",
# cubic
"T^1",         "T^2",       "T^3",         "T^4",       "T^5",         "T_h^1", 
"T_h^2",       "T_h^3",     "T_h^4",       "T_h^5",     "T_h^6",       "T_h^7", 
"O^1",         "O^2",       "O^3",         "O^4",       "O^5",         "O^6", 
"O^7",         "O^8",       "T_d^1",       "T_d^2",     "T_d^3",       "T_d^4", 
"T_d^5",       "T_d^6",     "O_h^1",       "O_h^2",     "O_h^3",       "O_h^4", 
"O_h^5",       "O_h^6",     "O_h^7",       "O_h^8",     "O_h^9",       "O_h^10"
]

# Hermann-Mauguin notation, ordered relative to space/plane group number
const hermannmauguin_table = Dict{Int64, Vector{String}}(
3 => # space-group notation (three dimensions) [from https://bruceravel.github.io/demeter/artug/atoms/space.html]
[
# triclinic
"P 1",        "P -1",
# monoclinic
"P 2",        "P 21",       "C 2",        "P M",        "P C",        "C M",        
"C C",        "P 2/M",      "P 21/M",     "C 2/M",      "P 2/C",      "P 21/C",     
"C 2/C",
# orthorhombic
"P 2 2 2",    "P 2 2 21",   "P 21 21 2",  "P 21 21 21", "C 2 2 21",   "C 2 2 2",
"F 2 2 2",    "I 2 2 2",    "I 21 21 21", "P M M 2",    "P M C 21",   "P C C 2",
"P M A 2",    "P C A 21",   "P N C 2",    "P M N 21",   "P B A 2",    "P N A 21",
"P N N 2",    "C M M 2",    "C M C 21",   "C C C 2",    "A M M 2",    "A B M 2",
"A M A 2",    "A B A 2",    "F M M 2",    "F D D 2",    "I M M 2",    "I B A 2",
"I M A 2",    "P M M M",    "P N N N",    "P C C M",    "P B A N",    "P M M A",
"P N N A",    "P M N A",    "P C C A",    "P B A M",    "P C C N",    "P B C M",
"P N N M",    "P M M N",    "P B C N",    "P B C A",    "P N M A",    "C M C M",
"C M C A",    "C M M M",    "C C C M",    "C M M A",    "C C C A",    "F M M M",
"F D D D",    "I M M M",    "I B A M",    "I B C A",    "I M M A",
# tetragonal 
"P 4",        "P 41",       "P 42",       "P 43",       "I 4",        "I 41",
"P -4",       "I -4",       "P 4/M",      "P 42/M",     "P 4/N",      "P 42/N",
"I 4/M",      "I 41/A",     "P 4 2 2",    "P 4 21 2",   "P 41 2 2",   "P 41 21 2",
"P 42 2 2",   "P 42 21 2",  "P 43 2 2",   "P 43 21 2",  "I 4 2 2",    "I 41 2 2",
"P 4 M M",    "P 4 B M",    "P 42 C M",   "P 42 N M",   "P 4 C C",    "P 4 N C",
"P 42 M C",   "P 42 B C",   "I 4 M M",    "I 4 C M",    "I 41 M D",   "I 41 C D",
"P -4 2 M",   "P -4 2 C",   "P -4 21 M",  "P -4 21 C",  "P -4 M 2",   "P -4 C 2",
"P -4 B 2",   "P -4N2",     "I -4 M 2",   "I -4 C 2",   "I -42 M",    "I -42 D",
"P 4/M M M",  "P 4/M C C",  "P 4/N B M",  "P 4/N N C",  "P 4/M B M",  "P 4/M N C",
"P 4/N M M",  "P 4/N C C",  "P 42/M M C", "P 42/M C M", "P 42/N B C", "P 42/N N M",
"P 42/M B C", "P 42/M N M", "P 42/N M C", "P 42/N C M", "I 4/M M M",  "I 4/M C M",
"I 41/A M D", "I 41/A C D",
# trigonal
"P 3",        "P 3 1",      "P 32",      "R 3",         "P -3",       "R -3",
"P 3 1 2",    "P 3 2 1",    "P 31 1 2",   "P 31 2 1",   "P 32 1 2",   "P 32 2 1",
"R 32",       "P 3 M 1",    "P 3 1 M",    "P 3 C 1",    "P 3 1 C",    "R 3 M",
"R 3C",       "P -3 1 M",   "P -3 1 C",   "P -3 M 1",   "P -3 C 1",   "R -3 M",
"R -3 C",
# hexagonal
"P 6",        "P 61",       "P 65",       "P 62",       "P 64",       "P 63",
"P -6",       "P 6/M",      "P 63/M",     "P 62 2",     "P 61 2 2",   "P 65 2 2",
"P 62 2 2",   "P 64 2 2",   "P 63 2 2",   "P 6 M M",    "P 6 C C",    "P 63 C M",
"P 63 M C",   "P -6 M 2",   "P -6 C 2",   "P -6 2 M",   "P -62 C",    "P 6/M M M",
"P 6/M C C",  "P 63/M C M", "P 63/M M C",
# cubic
"P 2 3",      "F 2 3",      "I 2 3",      "P 21 3",     "I 21 3",     "P M 3",
"P N 3",      "F M 3",      "F D 3",      "I M 3",      "P A 3",      "I A 3",
"P 4 3 2",    "P 42 3 2",   "F 4 3 2",    "F 41 3 2",   "I 4 3 2",    "P 43 3 2",
"P 41 3 2",   "I 41 3 2",   "P -4 3 M",   "F -4 3 M",   "I -4 3 M",   "P -4 3 N",
"F -4 3 C",   "I -4 3 D",   "P M 3 M",    "P N 3 N",    "P M 3 N",    "P N 3 M",
"F M 3 M",    "F M 3 C",    "F D 3 M",    "F D 3 C",    "I M 3 M",    "I A 3 D"
]
, 
2 => # plane group notation (two dimensions) [see e.g. Table 19 of Cracknell, Adv. Phys. 1974]
[
# oblique
"p 1",        "p 2 1 1",
# rectangular ('p' or 'c' centering; c-centered lattices are rhombic in their primitive cell)
"p 1 m 1",    "p 1 g 1",
"c 1 m 1",    "p 2 m m",
"p 2 m g",    "p 2 g g",
"c 2 m m",    
# square
"p 4",        "p 4 m m",
"p 4 g m",
# hexagonal
"p 3",        "p 3 m 1",
"p 3 1 m",    "p 6",
"p 6 m m"
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
    W = pg(op); w = translation(op); dim = size(W,1)
    if dim == 2 # we just augment the 2D case by leaving z invariant
        W = [W zeros(2); 0.0 0.0 1.0]; 
        w = [w, 0]
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
        if negidx ≠ nothing && negidx ≠ 1; u .*= -1; end

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