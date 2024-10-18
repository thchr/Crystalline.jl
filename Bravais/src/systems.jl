"""
    crystal(a, b, c, α, β, γ)  -->  DirectBasis{3}

Calculate basis vectors ``\\mathbf{R}_1``, ``\\mathbf{R}_2``, ``\\mathbf{R}_3`` in a 3D
Cartesian basis for a right-handed coordinate system with specified basis vector lengths
`a`, `b`, `c` (associated with ``\\mathbf{R}_1``, ``\\mathbf{R}_2``, & ``\\mathbf{R}_3``,
respectively) and specified interaxial angles 
`α` ``= ∠(\\mathbf{R}_2,\\mathbf{R}_3)``, `β` ``= ∠(\\mathbf{R}_3,\\mathbf{R}_1)``, 
`γ` ``= ∠(\\mathbf{R}_1,\\mathbf{R}_2)``, with ``∠`` denoting the angle between two vectors.

For definiteness, the ``\\mathbf{R}_1`` basis vector is oriented along the ``x``-axis of the
Cartesian coordinate system, and the ``\\mathbf{R}_2`` axis is placed in the ``xy``-plane.
"""
function crystal(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    # consistency checks on interaxial angles (equivalently, sides of the corresponding 
    # unit-spherical triangle)
    if !isvalid_sphericaltriangle(α,β,γ)
        throw(DomainError((α,β,γ), "The provided angles α,β,γ cannot be mapped to a spherical triangle, and thus do not form a valid axis system"))
    end
    # for many of the "special" input values (e.g., π/2, π/3), we can get better floating
    # point accuracy by dividing the angles by π and then using `cospi` instead of `cos`
    # etc.; e.g., this makes sure that if γ = π/2, then cosγ = 0.0, not ~6e-17.
    α′, β′, γ′ = α/π, β/π, γ/π
    sinγ, cosγ = sincospi(γ′)
    # R₁ and R₂ are easy
    R₁ = SVector{3,Float64}(a, 0.0, 0.0)
    R₂ = SVector{3,Float64}(b*cosγ, b*sinγ, 0.0)
    # R3 is harder
    cosα = cospi(α′)
    cosβ = cospi(β′)
    ϕ = atan(cosα - cosγ*cosβ, sinγ*cosβ)
    θ = asin(sign(β)*sqrt(cosα^2 + cosβ^2 - 2*cosα*cosγ*cosβ)/abs(sinγ)) # more stable than asin(cosβ/cosϕ) when β or γ ≈ π/2
    sinθ, cosθ = sincospi(θ/π)
    sinϕ, cosϕ = sincospi(ϕ/π)
    R₃ = SVector{3,Float64}(c.*(sinθ*cosϕ, sinθ*sinϕ, cosθ))

    Rs = DirectBasis(R₁, R₂, R₃)
    return Rs
end

"""
    crystal(a, b, γ)  -->  DirectBasis{2}

Calculate basis vectors ``\\mathbf{R}_1``, ``\\mathbf{R}_2`` in a 2D Cartesian basis for a 
right-handed coordinate system with specified basis vector lengths `a`, `b` (associated with
``\\mathbf{R}_1`` & ``\\mathbf{R}_2``, respectively) and specified interaxial angle
`γ` ``= ∠(\\mathbf{R}_1,\\mathbf{R}_2)``.

For definiteness, the ``\\mathbf{R}_1`` basis vector is oriented along the ``x``-axis of the
Cartesian coordinate system.
"""
function crystal(a::Real,b::Real,γ::Real) 
    R₁ = SVector{2,Float64}(a, 0.0)
    R₂ = SVector{2,Float64}(b.*(cos(γ), sin(γ)))

    return DirectBasis(R₁,R₂)
end

"""
    crystal(a)  -->  DirectBasis{1}
    
Return a one-dimensional crystal with lattice period `a`.
"""
crystal(a::Real) = DirectBasis(SVector{1,Float64}(float(a)))


# For a three-axis system, α, β, and γ are subject to constraints: specifically, 
# since they correspond to sides of a (unit-radius) spherical triangle, they 
# are subject to identical constraints. These constraints are
#     0 < α + β + γ < 2π,                           (1)
#     sin(s-α)*sin(s-β)*sin(s-γ)/sin(s) > 0,        (2)
# with s = (α + β + γ)/2. Constraint (2) can be identified from Eq. (38) of 
# http://mathworld.wolfram.com/SphericalTrigonometry.html; due to (1), it can 
# be simplified to sin(s-α)*sin(s-β)*sin(s-γ) > 0. This impacts generation 
# of triclinic and monoclinic crystals.
function isvalid_sphericaltriangle(α::Real, β::Real, γ::Real)
    s = (α+β+γ)/2
    check1 = zero(s) < s < π
    check2 = sin(s-α)*sin(s-β)*sin(s-γ) > zero(s)
    return check1 && check2
end

°(φ::Real) = deg2rad(φ)

""" 
    crystalsystem(Rs::DirectBasis{D})  -->  String
    crystalsystem(Gs::ReciprocalBasis{D})  -->  String

Determine the crystal system of a point lattice with `DirectBasis` `Rs`, assuming the
conventional setting choice defined in the International Tables of Crystallography [^ITA6].

If a `ReciprocalBasis` `Gs` is provided for the associated reciprocal point lattice, the
crystal system is determined by first transforming to the direct lattice.

There are 4 crystal systems in 2D and 7 in 3D (see Section 2.1.2(iii) of [^ITA5]):

| `D`    | System       | Conditions             | Free parameters      |
|:-------|--------------|------------------------|----------------------|
| **1D** | linear       | none                   | a                    |
| **2D** | square       | a=b & γ=90°            | a                    |
|        | rectangular  | γ=90°                  | a,b                  |
|        | hexagonal    | a=b & γ=120°           | a                    |
|        | oblique      | none                   | a,b,γ                |
| **3D** | cubic        | a=b=c & α=β=γ=90°      | a                    |
|        | hexagonal    | a=b & α=β=90° & γ=120° | a,c                  |
|        | trigonal     | a=b & α=β=90° & γ=120° | a,c (a,α for hR)     |
|        | tetragonal   | a=b & α=β=γ=90°        | a,c                  |
|        | orthorhombic | α=β=γ=90°              | a,b,c                |
|        | monoclinic   | α=γ=90°                | a,b,c,β≥90°          |
|        | triclinic    | none                   | a,b,c,α,β,γ          |

The `Rs` must specify a set of conventional basis vectors, i.e., not generally primitive.
For primitive basis vectors, the crystal system can be further reduced into 5 Bravais types
in 2D and 14 in 3D (see [`bravaistype`](@ref)).

[^ITA6]: M.I. Aroyo, International Tables of Crystallography, Vol. A, 6th ed. (2016): Tables
         3.1.2.1 and 3.1.2.2 (or Tables 2.1.2.1, 9.1.7.1, and 9.1.7.2 of [^ITA5]).

[^ITA5]: T. Hahn, International Tables of Crystallography, Vol. A, 5th ed. (2005).
"""
function crystalsystem(Rs::DirectBasis{D}) where D
    if D == 1
        # doesn't seem to exist a well-established convention for 1D? this is ours...
        system = "linear"
        
    elseif D == 2
        a,b = norm.(Rs)
        γ = angles(Rs)
        if a≈b && γ≈°(90)
            system = "square"
        elseif γ≈°(90)
            system = "rectangular"
        elseif a≈b && γ≈°(120)
            system = "hexagonal"
        else
            system = "oblique"
        end

            
    elseif D == 3 
        # TODO: Generalize this to work for non-standard orientations of the lattice/
        #       lattice-vector orderings
        a,b,c = norm.(Rs)
        α,β,γ = angles(Rs)
        if a≈b≈c && α≈β≈γ≈°(90)             # cubic        (cP, cI, cF)
            system = "cubic"
        elseif a≈b && γ≈°(120) && α≈β≈°(90) # hexagonal    (hR, hP)
            system = "hexagonal" 
        elseif a≈b≈c && α≈β≈γ               # trigonal     (? hP, hI ?)
            system = "trigonal"
                # rhombohedral axes                   (a = b = c, α=β=γ < 120° ≠ 90° ?)
                # hexagonal axes, triple obverse axes (a = b ≠ c, α=β=90°, γ=120° ?)
        elseif a≈b && α≈β≈γ≈°(90)           # tetragonal   (tP, tI) 
            system = "tetragonal"
        elseif α≈β≈γ≈°(90)                  # orthorhombic (oP, oI, oF, oC)
            system = "orthorhombic"
        elseif α≈γ≈°(90)                    # monoclinic   (mP, mC)
            system = "monoclinic"
        else                                # triclinic    (aP)
            system = "triclinic"
        end
    else
        throw(DomainError(D, "dimension must be 1, 2, or 3"))
    end
    return system
end
function crystalsystem(Gs::ReciprocalBasis{D}) where D
    Rs = DirectBasis{D}(reciprocalbasis(Gs).vs)
    return crystalsystem(Rs)
end


function crystalsystem(sgnum::Integer, D::Integer=3)
    if D == 1
        # doesn't seem to exist a well-established convention for 1D? this is ours...
        if      sgnum ∈ 1:2;   return "linear"       # lp
        else    _throw_invalid_sgnum(sgnum, D)
        end
    elseif D == 2
        if      sgnum ∈ 1:2;   return "oblique"      # mp
        elseif  sgnum ∈ 3:9;   return "rectangular"  # op, oc
        elseif  sgnum ∈ 10:12; return "square"       # tp
        elseif  sgnum ∈ 13:17; return "hexagonal"    # hp
        else    _throw_invalid_sgnum(sgnum, D)
        end
    
    elseif D == 3
        if      sgnum ∈ 1:2;     return "triclinic"     # aP
        elseif  sgnum ∈ 3:15;    return "monoclinic"    # mP, mC
        elseif  sgnum ∈ 16:74;   return "orthorhombic"  # oP, oI, oF, oC
        elseif  sgnum ∈ 75:142;  return "tetragonal"    # tP, tI
        elseif  sgnum ∈ 143:167; return "trigonal"      # hR, hP
        elseif  sgnum ∈ 168:194; return "hexagonal"     # hR, hP
        elseif  sgnum ∈ 195:230; return "cubic"         # cP, cI, cF
        else    _throw_invalid_sgnum(sgnum, D)
        end
    end
end


""" 
    directbasis(sgnum, D=3;    abclims, αβγlims)
    directbasis(sgnum, Val(D); abclims, αβγlims) --> DirectBasis{D}

Return a random (conventional) `DirectBasis` for a crystal compatible with the space group
number `sgnum` and dimensionality `D`.
Free parameters in the lattice vectors are chosen randomly, with limits optionally supplied
in `abclims` (lengths) and `αβγlims` (angles).
By convention, the length of the first lattice vector (`a`) is set to unity, such that the
second and third (`b` and `c`) lattice vectors' lengths are relative to the first.

Limits on the relative uniform distribution of lengths `b` and `c` can be specified as 
2-tuple kwarg `abclims`; similarly, limits on the angles `α`, `β`, `γ` can be set via
αβγlims (only affects oblique, monoclinic, & triclinic lattices).
"""
function directbasis(sgnum::Integer, Dᵛ::Val{D}=Val(3);
                     abclims::NTuple{2,Real}=(0.5,2.0), 
                     αβγlims::NTuple{2,Real}=(°(30),°(150))) where D
    system = crystalsystem(sgnum, D)
    if D == 1
        a = 1.0
        return crystal(a)
    elseif D == 2
        if     system == "square"      # a=b & γ=90° (free: a)
            a = b = 1.0
            γ = °(90)
        elseif system == "rectangular" # γ=90° (free: a,b)
            a = 1.0;    b = relrand(abclims)
            γ = °(90)           
        elseif system == "hexagonal"   # a=b & γ=120° (free: a)
            a = b = 1.0;
            γ = °(120)
        elseif system == "oblique"     # no conditions (free: a,b,γ)
            a = 1.0;    b = relrand(abclims)
            γ = uniform_rand(αβγlims...)
        else 
            throw(DomainError(system))
        end
        return crystal(a,b,γ)

    elseif D == 3
        if     system == "cubic"        # a=b=c & α=β=γ=90° (free: a)
            a = b = c = 1.0
            α = β = γ = °(90)
        elseif system == "hexagonal" || # a=b & α=β=90° & γ=120° (free: a,c)
               system == "trigonal"    
            a = b = 1.0;        c = relrand(abclims)
            α = β = °(90);      γ = °(120)
            # For the trigonal crystal system, a convention is adopted where 
            # the crystal basis matches the a hexagonal one, even when the 
            # Bravais type is rhombohedral (of course, the primitive basis 
            # differs). The conventional cell is always chosen to have (triple
            # obverse) hexagonal axes (see ITA6 Sec. 3.1.1.4, Tables 2.1.1.1, 
            # & 3.1.2.2).
            # For rhombohedral systems the primitive cell has a=b=c, α=β=γ<120°≠90°.
            # Note that the hexagonal and trigonal crystal systems also share 
            # the same crystal system abbreviation 'h' (see CRYSTALSYSTEM_ABBREV),
            # which already suggests this choice.
            # TODO: in principle, this means it would be more meaningful to 
            # branch on `CRYSTALSYSTEM_ABBREV[D][system]` than on `system`.
        elseif system == "tetragonal"   # a=b & α=β=γ=90° (free: a,c)
            a = b = 1.0;        c = relrand(abclims)
            α = β = γ = °(90)
        elseif system == "orthorhombic" # α=β=γ=90° (free: a,b,c)
            a = 1.0;            b, c = relrand(abclims), relrand(abclims)
            α = β = γ = °(90)
        elseif system == "monoclinic"   # α=γ=90° (free: a,b,c,β≥90°)
            a = 1.0;            b, c = relrand(abclims), relrand(abclims)
            α = γ = °(90);      β = uniform_rand(°(90), αβγlims[2])
            while !isvalid_sphericaltriangle(α,β,γ)
                # arbitrary combinations of α,β,γ need not correspond to a valid
                # axis-system; reroll until they do
                β = uniform_rand(°(90), αβγlims[2])
            end
        elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
            a = 1.0;            b, c = relrand(abclims), relrand(abclims)
            low, high = αβγlims
            α, β, γ = uniform_rand(low, high), uniform_rand(low, high), uniform_rand(low, high)
            while !isvalid_sphericaltriangle(α,β,γ)
                # arbitrary combinations of α,β,γ need not correspond to a valid
                # axis-system; reroll until they do
                α, β, γ = uniform_rand(low, high), uniform_rand(low, high), uniform_rand(low, high)
            end
        else 
            throw(DomainError(system))
        end        
        return crystal(a,b,c,α,β,γ)

    else 
        _throw_invalid_dim(D)
    end
end
function directbasis(sgnum::Integer, D::Integer;
            abclims::NTuple{2,Real}=(0.5,2.0), αβγlims::NTuple{2,Real}=(°(30),°(150)))
    directbasis(sgnum, Val(D); abclims=abclims, αβγlims=αβγlims)
end

const CRYSTALSYSTEM_ABBREV = (
    ImmutableDict("linear"=>'l'),                                                       # 1D
    ImmutableDict("oblique"=>'m', "rectangular"=>'o', "square"=>'t', "hexagonal"=>'h'), # 2D
    ImmutableDict("triclinic"=>'a', "monoclinic"=>'m', "orthorhombic"=>'o',             # 3D
        "tetragonal"=>'t', "trigonal"=>'h', "hexagonal"=>'h', "cubic"=>'c')
    )



# cached results of combining `crystalsystem(sgnum, D)` w/ a centering calc from `iuc`
const BRAVAISTYPE_2D = (
    "mp", "mp", "op", "op", "oc", "op", "op", "op", "oc", "tp", "tp", "tp", "hp", "hp",
    "hp", "hp", "hp")
const BRAVAISTYPE_3D = (
    "aP", "aP", "mP", "mP", "mC", "mP", "mP", "mC", "mC", "mP", "mP", "mC", "mP", "mP",
    "mC", "oP", "oP", "oP", "oP", "oC", "oC", "oF", "oI", "oI", "oP", "oP", "oP", "oP",
    "oP", "oP", "oP", "oP", "oP", "oP", "oC", "oC", "oC", "oA", "oA", "oA", "oA", "oF",
    "oF", "oI", "oI", "oI", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP", "oP",
    "oP", "oP", "oP", "oP", "oP", "oP", "oC", "oC", "oC", "oC", "oC", "oC", "oF", "oF",
    "oI", "oI", "oI", "oI", "tP", "tP", "tP", "tP", "tI", "tI", "tP", "tI", "tP", "tP",
    "tP", "tP", "tI", "tI", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI",
    "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI", "tI", "tI", "tP", "tP",
    "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI", "tI", "tI", "tP", "tP", "tP", "tP",
    "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tP", "tI", "tI",
    "tI", "tI", "hP", "hP", "hP", "hR", "hP", "hR", "hP", "hP", "hP", "hP", "hP", "hP",
    "hR", "hP", "hP", "hP", "hP", "hR", "hR", "hP", "hP", "hP", "hP", "hR", "hR", "hP",
    "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP",
    "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "hP", "cP", "cF",
    "cI", "cP", "cI", "cP", "cP", "cF", "cF", "cI", "cP", "cI", "cP", "cP", "cF", "cF",
    "cI", "cP", "cP", "cI", "cP", "cF", "cI", "cP", "cF", "cI", "cP", "cP", "cP", "cP",
    "cF", "cF", "cF", "cF", "cI", "cI")
    
"""
    bravaistype(sgnum::Integer, D::Integer=3; normalize::Bool=false)  -->  String

Return the Bravais type of `sgnum` in dimension `D` as a string (as the concatenation
of the single-character crystal abbreviation and the centering type).

## Keyword arguments

- **`normalize`:** If the centering type associated with `sgnum` is `'A'`, we can choose
  (depending on the keyword argument `normalize`, defaulting to `false`) to "normalize" to
  the centering type `'C'`, since the difference between `'A'` and `'C'` centering only
  amounts to a basis change.
  With `normalize=true` we then have only the canonical 14 Bravais type, i.e. 
  `unique(bravaistype.(1:230, 3), normalize=true)` returns only 14 distinct types, rather
  than 15.

  This only affects space groups 38-41 (normalizing their conventional Bravais types from
  `"oA"` to `"oC"`).
"""
function bravaistype(sgnum::Integer, D::Integer=3; normalize::Bool=false)
    @boundscheck boundscheck_sgnum(sgnum, D)
    if D == 3
        # If the centering type is 'A', then we could in fact always pick the basis
        # differently such that the centering would be 'C'; in other words, base-centered
        # lattices at 'A' and 'C' in fact describe the same Bravais lattice; there is no
        # significance in trying to differentiate them - if we do, we end up with 15
        # Bravais lattices in 3D rather than 14: so we manually fix that here:
        bt = @inbounds BRAVAISTYPE_3D[sgnum]
        if normalize
            return bt != "oA" ? bt : "oC"
        else
            return bt
        end
    elseif D == 2
        return @inbounds BRAVAISTYPE_2D[sgnum]
    else
        return "lp"
    end
end

"""
    centering(sgnum::Integer, D::Integer=3)  -->  Char

Return the conventional centering type `cntr` of the space group with number `sgnum` and
dimension `D`.

The centering type is equal to the first letter of the Hermann-Mauguin notation's label,
i.e., `centering(sgnum, D) == first(Crystalline.iuc(sgnum, D))`. Equivalently, the
centering type is the second and last letter of the Bravais type ([`bravaistype`](@ref)),
i.e., `centering(sgnum, D) == bravaistype(sgnum, D)`.

Possible values of `cntr`, depending on dimensionality `D`, are (see ITA Sec. 9.1.4):

- `D = 1`:
    - `cntr = 'p'`: no centering (primitive)
- `D = 2`:
    - `cntr = 'p'`: no centring (primitive)
    - `cntr = 'c'`: face centered
- `D = 3`: 
    - `cntr = 'P'`: no centring (primitive)
    - `cntr = 'I'`: body centred (innenzentriert)
    - `cntr = 'F'`: all-face centred
    - `cntr = 'A'` or `'C'`: one-face centred, (b,c) or (a,b)
    - `cntr = 'R'`: hexagonal cell rhombohedrally centred
"""
centering(sgnum::Integer, D::Integer=3) = last(bravaistype(sgnum, D))