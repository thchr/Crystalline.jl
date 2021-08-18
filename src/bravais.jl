"""
    crystal(a,b,c,α,β,γ) --> Rs::DirectBasis{3}

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
    # consistency checks on interaxial angles (equivalently, sides of the corresponding unit-spherical triangle)
    if !isvalid_sphericaltriangle(α,β,γ)
        throw(DomainError((α,β,γ), "The provided angles α,β,γ cannot be mapped to a spherical triangle, and thus do not form a valid axis system"))
    end
    sinγ, cosγ = sincos(γ)
    # R₁ and R₂ are easy
    R₁ = SVector{3,Float64}(a, 0.0, 0.0)
    R₂ = SVector{3,Float64}(b*cosγ, b*sinγ, 0.0)
    # R3 is harder
    cosα = cos(α)
    cosβ = cos(β)
    ϕ = atan(cosα - cosγ*cosβ, sinγ*cosβ)
    θ = asin(sign(β)*sqrt(cosα^2 + cosβ^2 - 2*cosα*cosγ*cosβ)/abs(sin(γ))) # more stable than asin(cosβ/cosϕ) when β or γ ≈ π/2
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    R₃ = SVector{3,Float64}(c.*(sinθ*cosϕ, sinθ*sinϕ, cosθ))

    Rs = DirectBasis(R₁, R₂, R₃)
    return Rs
end

"""
    crystal(a,b,γ) --> DirectBasis{2}

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
    crystal(a)  --> DirectBasis{1}
    
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
    crystalsystem(R::DirectBasis{D})

Determine the crystal system of a point lattice specified in a 
*conventional* `DirectBasis{D}` of dimension `D` using
Tables 2.1.2.1, 9.1.7.1, & 9.1.7.2 of the International Tables of 
Crystallography, Volume 1 (ITA). 
There are 4 crystal systems in 2D and 7 in 3D (see ITA 2.1.2(iii)):

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

The DirectBasis input is assumed to use *conventional* basis vectors; 
i.e. not necessarily primitive. For primitive basis vectors, the 
crystal system can be further reduced into 5 Bravais types in 2D and
14 in 3D.
"""
function crystalsystem(Rs::DirectBasis{D}) where D
    if D == 1
        # doesn't seem to exist a well-established convention for 1D? this is ours...
        system = "linear"
        
    elseif D == 2
        a,b = norms(Rs)
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
        # TODO: Generalize this to work for non-standard orientations of the lattice/lattice-vector orderings.
        #       If that is done, this can be safely applied to ReciprocalBasis as well (only a problem in 3D).
        a,b,c = norms(Rs)
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
        throw(DomainError("Support for higher-than-three dimensional lattices not implemented."))
    end
    return system
end


function crystalsystem(sgnum::Integer, D::Integer=3)
    if D == 1
        # doesn't seem to exist a well-established convention for 1D? this is ours...
        if      sgnum ∈ 1:2;   return "linear"       # lp
        else    throw(DomainError(sgnum, "There are only 2 one-dimensional line groups."))
        end
    elseif D == 2
        if      sgnum ∈ 1:2;   return "oblique"      # mp
        elseif  sgnum ∈ 3:9;   return "rectangular"  # op, oc
        elseif  sgnum ∈ 10:12; return "square"       # tp
        elseif  sgnum ∈ 13:17; return "hexagonal"    # hp
        else    throw(DomainError(sgnum, "There are only 17 two-dimensional plane groups."))
        end
    
    elseif D == 3
        if      sgnum ∈ 1:2;     return "triclinic"     # aP
        elseif  sgnum ∈ 3:15;    return "monoclinic"    # mP, mC
        elseif  sgnum ∈ 16:74;   return "orthorhombic"  # oP, oI, oF, oC
        elseif  sgnum ∈ 75:142;  return "tetragonal"    # tP, tI
        elseif  sgnum ∈ 143:167; return "trigonal"      # hR, hP
        elseif  sgnum ∈ 168:194; return "hexagonal"     # hR, hP
        elseif  sgnum ∈ 195:230; return "cubic"         # cP, cI, cF
        else    throw(DomainError(sgnum, "There are only 230 three-dimensional space groups."))
        end
    end
end

"""
    relrand(lims::NTuple{2,Real}) --> Float64

Computes a random number in the range specified by the two-element 
tuple `lims`. The random numbers are sampled from two uniform 
distributions, namely [`lims[1]`, 1] and [1, `lims[2]`], in such a
way as to ensure that the sampling is uniform over the joint  
interval [-1/`lims[1]`, -1] ∪ [1, `lims[2]`].

This is useful for ensuring an even sampling of numbers that are
either smaller or larger than unity. E.g. for `x = relrand((0.2,5.0))`,
`x` is equally probable to fall in inv(`x`)∈[1,5] or `x`∈[1,5].
"""
function relrand(lims::NTuple{2,<:Real})
    low, high = lims; invlow = inv(low)
    lowthres = (invlow - 1.0)/(invlow + high - 2.0)
    if rand() < lowthres && low < 1.0   # smaller than 1.0
        r = rand(_Uniform(low,1.0))
    elseif high > 1.0                   # bigger than 1.0
        r = rand(_Uniform(1.0,high))
    else                                # default
        return rand(_Uniform(low,high))
    end
end

""" 
    directbasis(sgnum, D=3;    abclims, αβγlims)
    directbasis(sgnum, Val(D); abclims, αβγlims) --> DirectBasis{D}

Return a (conventional) `DirectBasis` for a crystal compatible with 
the space group number `sgnum` and dimensionality `D`. Free parameters
in the lattice vectors are chosen randomly, with limits optionally
supplied in `abclims` (lengths) and `αβγlims` (angles).
By convention, the length of the first lattice vector (`a`) is set
to unity, such that the second and third (`b` and `c`) lattice 
vectors' lengths are relative to the first.

Limits on the relative uniform distribution of lengths `b` and `c`
can be specified as 2-tuple kwarg `abclims`; similarly, limits on 
the angles `α`, `β`, `γ` can be set via αβγlims (only affects 
oblique, monoclinic, & triclinic lattices).
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
            γ = rand(_Uniform(αβγlims...)) 
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
            # obverse) hexagonal axes (see ITA7 Sec. 3.1.1.4, Tables 2.1.1.1, 
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
            α = γ = °(90);      β = rand(_Uniform(°(90), αβγlims[2]))
            while !isvalid_sphericaltriangle(α,β,γ)  # arbitrary combinations of α,β,γ may not correspond 
                β = rand(_Uniform(°(90), αβγlims[2])) # to a valid axis-system; reroll until they do
            end
        elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
            a = 1.0;            b, c = relrand(abclims), relrand(abclims)
            U = _Uniform(αβγlims...)
            α, β, γ = rand(U), rand(U), rand(U)
            while !isvalid_sphericaltriangle(α,β,γ) # arbitrary combinations of α,β,γ may not correspond 
                α, β, γ = rand(U), rand(U), rand(U) # to a valid axis-system; reroll until they do
            end
        else 
            throw(DomainError(system))
        end        
        return crystal(a,b,c,α,β,γ)

    else 
        _throw_invaliddim(D)
    end
end
function directbasis(sgnum::Integer, D::Integer;
            abclims::NTuple{2,Real}=(0.5,2.0), αβγlims::NTuple{2,Real}=(°(30),°(150)))
    directbasis(sgnum, Val(D); abclims=abclims, αβγlims=αβγlims)
end


const CRYSTALSYSTEM_ABBREV = (ImmutableDict("linear"=>'l'),                                            # 1D
                              ImmutableDict("oblique"=>'m', "rectangular"=>'o', "square"=>'t',         # 2D
                                   "hexagonal"=>'h'),
                              ImmutableDict("triclinic"=>'a', "monoclinic"=>'m', "orthorhombic"=>'o',  # 3D
                                   "tetragonal"=>'t', "trigonal"=>'h', "hexagonal"=>'h', 
                                   "cubic"=>'c')
                             )

"""
    $(TYPEDSIGNATURES)

Return the Bravais type of `sgnum` in dimension `D` as a string (as the concatenation
of the single-character crystal abbreviation and the centering type).

## Keyword arguments

If the centering type associated with `sgnum` is `'A'`, we can choose (depending on the
keyword argument `normalize`, defaulting to `false`) to "normalize" to the centering type
`'C'`, since the difference between `A` and `C` centering only amounts to a basis change.
With `normalize=true` we then have only the canonical 14 Bravais type, i.e. 
`unique(bravaistype.(1:230, 3), normalize=true)` returns only 14 distinct types, rather
than 15.

This only affects space groups 38-41 (normalizing their conventional Bravais types from
`oA` to `oC`).
"""
@inline function bravaistype(sgnum::Integer, D::Integer=3; normalize::Bool=false)
    cntr = centering(sgnum, D)
    system = crystalsystem(sgnum, D)

    # If the centering type is 'A', then we could in fact always pick
    # the basis differently such that the centering would be 'C'; in 
    # other words, base-centered lattices at 'A' and 'C' in fact describe
    # the same Bravais lattice; there is no significance in trying to 
    # differentiate them - if we do, we end up with 15 Bravais lattices in 
    # 3D rather than 14: so we manually fix that here:
    if normalize
        cntr = cntr == 'A' ? 'C' : cntr
    end

    # pick the correct crystal system abbreviation from CRYSTALSYSTEM_ABBREV 
    # and return its concatenation with the (now-"normalized") centering type
    return CRYSTALSYSTEM_ABBREV[D][system]*cntr 
end



# Transformation matrices 𝐏 from a (direct-space) conventional basis 
# (𝐚 𝐛 𝐜) to a primitive basis (𝐚ₚ 𝐛ₚ 𝐜ₚ) via
#     (𝐚ₚ 𝐛ₚ 𝐜ₚ) = (𝐚 𝐛 𝐜)𝐏
# # with (𝐚 𝐛 𝐜) and (𝐚ₚ 𝐛ₚ 𝐜ₚ) interpreted as column matrices the 
# transformation matrix 𝐏 depends only on the centering type [note
# that centering type 'B' seems to not occur, by convention]
# The values of 𝐏 are taken from Table 2 of the Aroyo's Bilbao 
# publication (https://doi.org/:10.1107/S205327331303091X), which 
# give the coefficients of (𝐏ᵀ)⁻¹. See also Hinuma's 2016 paper
# (https://doi.org/10.1016/j.commatsci.2016.10.015) for details,
# though note that they use different matrices for 'A' and complicate
# the 'C' scenario (Table 3).
const PRIMITIVE_BASIS_MATRICES = (
    # 1D
    ImmutableDict('p'=>SMatrix{1,1,Float64}(1)),                # primitive
    # 2D
    ImmutableDict('p'=>SMatrix{2,2,Float64}([1 0; 0 1]),        # primitive/simple
                  'c'=>SMatrix{2,2,Float64}([1 1; -1 1]./2)),   # centered      
    # 3D
    ImmutableDict(
        'P'=>SMatrix{3,3,Float64}([1 0 0; 0 1 0; 0 0 1]),       # primitive/simple
        'F'=>SMatrix{3,3,Float64}([0 1 1; 1 0 1; 1 1 0]./2),    # face-centered
        'I'=>SMatrix{3,3,Float64}([-1 1 1; 1 -1 1; 1 1 -1]./2), # body-centered
        'R'=>SMatrix{3,3,Float64}([2 -1 -1; 1 1 -2; 1 1 1]./3), # rhombohedrally-centered
        'A'=>SMatrix{3,3,Float64}([2 0 0; 0 1 -1; 0 1 1]./2),   # base-centered (along x)
        'C'=>SMatrix{3,3,Float64}([1 1 0; -1 1 0; 0 0 2]./2))   # base-centered (along z)
    )

"""
    primitivebasismatrix(cntr::Char, ::Val{D}=Val(3)) -> SMatrix{D,D,Float64}

Return the transformation matrix `P` that transforms a conventional unit cell with centering
`cntr` to the corresponding primitive unit cell (in dimension `D`).
"""
@inline function primitivebasismatrix(cntr::Char, ::Val{D}=Val(3)) where D
    D∉1:3 && _throw_invaliddim(D)
    return PRIMITIVE_BASIS_MATRICES[D][cntr]
end

@inline function centeringtranslation(cntr::Char, ::Val{D}=Val(3)) where D
    if D == 3
        if cntr == 'P';     return zeros(SVector{3})
        elseif cntr == 'I'; return SVector((1,1,1)./2)
        elseif cntr == 'F'; return SVector((1,0,1)./2)
        elseif cntr == 'R'; return SVector((2,1,1)./3)
        elseif cntr == 'A'; return SVector((0,1,1)./2)
        elseif cntr == 'C'; return SVector((1,1,0)./2)
        else;               _throw_invalidcntr(cntr)
        end
    elseif D == 2
        if cntr == 'p';     return zeros(SVector{2})
        elseif cntr == 'c'; return SVector((1,1)./2)
        else;               _throw_invalidcntr(cntr)
        end
    elseif D == 1
        return zeros(SVector{1})
    else 
        _throw_invaliddim(D)
    end
end


"""
    reciprocalbasis(Rs::DirectBasis{D}) --> Gs::ReciprocalBasis{D}
    
Calculates the reciprocal basis associated with a `DirectBasis` `Rs`
(alternatively supplied as an `NTuple` of `Vector`s). Returns a `ReciprocalBasis`.
"""
function reciprocalbasis(Rs::Union{DirectBasis{D}, NTuple{D, Vector{<:Real}}}) where D
    if D == 3
        G₁′ = Rs[2]×Rs[3]
        pref = 2π/dot(Rs[1], G₁′)
        vecs = pref .* (G₁′, Rs[3]×Rs[1], Rs[1]×Rs[2])
    elseif D == 2
        G₁′ = (@SVector [-Rs[2][2], Rs[2][1]])
        pref = 2π/dot(Rs[1], G₁′)
        vecs = pref .* (G₁′, (@SVector [Rs[1][2], -Rs[1][1]]))
    elseif D == 1
        vecs = (SVector{1,Float64}(2π/first(Rs[1])),)
    else
        # The general definition of the reciprocal basis is [G₁ ... Gₙ]ᵀ = 2π[R₁ ... Rₙ]⁻¹; 
        # that form should generally be a bit slower than the above specific variants, cf. 
        # the inversion operation, so we only use it as a high-dimensional fallback. Since 
        # we use SVectors, however, either approach will probably have the same performance.
        Rm = basis2matrix(Rs)
        Gm = 2π.*inv(transpose(Rm))
        vecs = ntuple(i->Gm[:,i], Val(D))
    end

    return ReciprocalBasis{D}(vecs)
end


""" 
    primitivize(Vs::Basis, sgnum::Integer) --> Rs′::Basis

Transforms a conventional `Basis` (either `DirectBasis` or `ReciprocalBasis`) `Vs`
into its primitive equivalent `Vs′`, provided that its centering differs from
the conventional (P or p), by inferring the Bravais type from the space group number
`sgnum` and applying an applying an appropriate (`Basis`-type specific) transformation. 
"""
function primitivize(Vs::Basis{D}, sgnum::Integer) where D
    cntr = centering(sgnum, D)
    return primitivize(Vs, cntr)
end

function transform(Rs::DirectBasis{D}, P::AbstractMatrix{<:Real}) where D
    # Rm′ = Rm*P (w/ Rm a matrix w/ columns of untransformed direct basis vecs Rᵢ)
    Rm′ = basis2matrix(Rs)*P
    return DirectBasis{D}(ntuple(i->Rm′[:,i], Val(D)))
end

"""
    primitivize(Rs::DirectBasis, cntr::Char) --> Rs′::DirectBasis

Transforms a conventional DirectBasis `Rs` into its primitive 
equivalent `Rs′`, with the transformation dependent on the centering
type `cntr` (P, I, F, R, A, C, and p, c); for centering P and p, the 
conventional and primive bases coincide.
"""
function primitivize(Rs::DirectBasis{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return Rs
    else         
        P = primitivebasismatrix(cntr, Val(D))
        # Rm′ = Rm*P (w/ Rm a matrix w/ columns of conventional direct basis vecs Rᵢ)
        return transform(Rs, P)
    end  
end

"""
    conventionalize(Rs′::DirectBasis, cntr::Char) --> Rs::DirectBasis

Transforms a primitive DirectBasis `Rs′` into its conventional equivalent `Rs`, with the 
transformation dependent on the centering type `cntr` (P, I, F, R, A, C, and p, c); for
centering P and p, the conventional and primive bases coincide.
"""
function conventionalize(Rs′::DirectBasis{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return Rs′
    else         
        P = primitivebasismatrix(cntr, Val(D))
        # Rm = Rm′*P⁻¹ (w/ Rm′ a matrix w/ columns of primitive direct basis vecs Rᵢ′)
        return transform(Rs′, inv(P)) 
    end  
end

function transform(Gs::ReciprocalBasis{D}, P::AbstractMatrix{<:Real}) where D
    # While the direct basis (𝐚 𝐛 𝐜) transforms like 
    #       (𝐚′ 𝐛′ 𝐜′) = (𝐚 𝐛 𝐜)𝐏
    # under a basis change matrix 𝐏, the reciprocal basis (𝐚* 𝐛* 𝐜*) transforms like 
    #       (𝐚*′ 𝐛*′ 𝐜*′) = (𝐚* 𝐛* 𝐜*)(𝐏⁻¹)ᵀ
    # since (𝐚 𝐛 𝐜)(𝐚* 𝐛* 𝐜*)ᵀ = 2π𝐈 must be conserved after the basis change

    # Gm′ = Gm*(P⁻¹)ᵀ = Gm*(Pᵀ)⁻¹ (w/ Gm a matrix w/ columns of untransformed reciprocal
    # vecs Gᵢ)
    Gm′ = basis2matrix(Gs)/P'
    return ReciprocalBasis{D}(ntuple(i->Gm′[:,i], Val(D)))
end

"""
    primitivize(Gs::ReciprocalBasis, cntr::Char) --> Gs′::ReciprocalBasis
    
Calculates the **primitive** reciprocal basis associated with an assumed **conventional**
reciprocal basis `Gs::ReciprocalBasis` with centering type `cntr`.
"""
function primitivize(Gs::ReciprocalBasis{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return Gs
    else         
        P = primitivebasismatrix(cntr, Val(D))        
        return transform(Gs, P)
    end
end
function conventionalize(Gs′::ReciprocalBasis{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return Gs
    else         
        P = primitivebasismatrix(cntr, Val(D))        
        return transform(Gs′, inv(P))
    end
end

# Note that the _coefficients_ of a general 𝐤-vector transform
# differently than the reciprocal _basis_, which transforms
# from non-primed to primed variants via a basis matrix 𝐏
# according to (see also `primitivize(Gs::ReciprocalBasis)`):
# Specifically, a 𝐤-vector is specified by a product of a
# reciprocal basis (𝐚* 𝐛* 𝐜*) and a coefficient vector
# (k₁ k₂ k₃)ᵀ, ie. 𝐤 ≡ (𝐚* 𝐛* 𝐜*)(k₁ k₂ k₃)ᵀ [note that 
# (k₁ k₂ k₃)ᵀ is a column vector].
# As a result, (k₁ k₂ k₃)ᵀ transforms like 
#     (k₁′ k₂′ k₃′)ᵀ = Pᵀ (k₁ k₂ k₃)ᵀ
# since
#     𝐤 = (𝐚*′ 𝐛*′ 𝐜*′)(k₁′ k₂′ k₃′)ᵀ     (1)  [... by definition]
#       = (𝐚* 𝐛* 𝐜*)(𝐏⁻¹)ᵀ(k₁′ k₂′ k₃′)ᵀ       [... transformation of (𝐚* 𝐛* 𝐜*) under 𝐏]
#       = (𝐚* 𝐛* 𝐜*)(k₁ k₂ k₃)ᵀ           (2)  [... by definition]
# then, combining (1) and (2)
#     (𝐏⁻¹)ᵀ(k₁′ k₂′ k₃′)ᵀ = (k₁ k₂ k₃)ᵀ
#  ⇔ (k₁′ k₂′ k₃′)ᵀ = 𝐏ᵀ(k₁ k₂ k₃)ᵀ

@doc raw"""
    transform(v::T, P::AbstractMatrix) where T<:AbstractVec --> v′::T

Return a transformed coordinate vector `v′` from an original coordinate vector `v` using a
basis change matrix `P`.

Note that a basis change matrix `P` transforms direct coordinate vectors ([`RVec`](@ref)) as
``\mathbf{r}' = \mathbf{P}^{-1}\mathbf{r}`` but transforms reciprocal coordinates
([`KVec`](@ref)) as ``\mathbf{k}' = \mathbf{P}^{\mathrm{T}}\mathbf{k}`` (see e.g. ITA7 
Secs. 1.5.1.2 and 1.5.2.1).
"""
function transform(kv::KVec{D}, P::AbstractMatrix{<:Real}) where D
    k₀, kabc = parts(kv)
    return KVec{D}(P'*k₀, P'*kabc)
end

function transform(rv::RVec{D}, P::AbstractMatrix{<:Real}) where D
    rcnst, rfree = parts(rv)
    return RVec{D}(P\rcnst, P\rfree)
end

@doc raw"""
    primitivize(v::T, cntr::Char) where T<:AbstractVec{D} --> v′::T

Transforms a conventional coordinate vector `v` to a standard primitive basis (specified by
the centering type `cntr`), returning the primitive coordinate vector `v′`.

Note that a basis change matrix ``P`` (as returned e.g. by
[`Crystalline.primitivebasismatrix`](@ref)) transforms direct coordinate vectors
([`RVec`](@ref)) as ``\mathbf{r}' = \mathbf{P}^{-1}\mathbf{r}`` but transforms reciprocal
coordinates ([`KVec`](@ref)) as ``\mathbf{k}' = \mathbf{P}^{\mathrm{T}}\mathbf{k}`` (see
e.g. ITA7 Sec. 1.5.1.2 and 1.5.2.1).
Recall also the distinction between transforming a basis and the coordinates of a vector.
"""
function primitivize(v::AbstractVec{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p'
        return v
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(v, P)
    end
end

"""
    conventionalize(v′::T, cntr::Char) where T<:AbstractVec{D} --> v::T

Transforms a primitive coordinate vector `v′` back to a standard conventional basis
(specified by the centering type `cntr`), returning the conventional coordinate vector `v`.

See also [`primitivize(::AbstractVec, ::Char`](@ref) and `transform`.
"""
function conventionalize(v′::AbstractVec{D}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p'
        return v′
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(v′, inv(P))
    end
end