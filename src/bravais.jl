"""
    crystal(a,b,c,α,β,γ) -> Crystal{3}

Calculate basis vectors `R1`, `R2`, `R3` in a 3D Cartesian basis 
for a right-handed coordinate system with specified basis vector lengths 
`a`, `b`, `c` (associated with  `R1`, `R2`, `R3`, respectively)
and specified interaxial angles `α=angle(R2, R3)`, `β=angle(R3,R1)`, 
`γ=angle(R1,R2)`.

For definiteness, the `R1` basis vector is oriented along the
x-axis of the Cartesian coordinate system, and the `R2` axis is 
placed in the xy-plane.
"""
function crystal(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real)
    # consistency checks on interaxial angles (equivalently, sides of the corresponding unit-spherical triangle)
    if !isvalid_sphericaltriangle(α,β,γ)
        throw(DomainError((α,β,γ), "The provided angles α,β,γ cannot be mapped to a spherical triangle, and thus do not form a valid axis system"))
    end
    # R1 and R2 are easy
    R1 = Float64[a, 0, 0] 
    R2 = b.*[cos(γ), sin(γ), 0]
    # R3 is harder
    cosα = cos(α)
    cosβ = cos(β)
    sinγ,cosγ = sincos(γ)
    ϕ = atan(cosα - cosγ*cosβ, sinγ*cosβ)
    θ = asin(sign(β)*sqrt(cosα^2 + cosβ^2 -2*cosα*cosγ*cosβ)/abs(sin(γ))) # more stable than asin(cosβ/cosϕ) when β or γ ≈ π/2
    sinθ,cosθ = sincos(θ)
    sinϕ,cosϕ = sincos(ϕ)
    R3 = c.*[sinθ*cosϕ, sinθ*sinϕ, cosθ]

    C = Crystal((R1,R2,R3))
    return C
end

"""
    crystal(a,b,γ) -> Crystal{2}

Calculate basis vectors `R1`, `R2` in a 2D Cartesian basis for a 
right-handed coordinate system with specified basis vector lengths 
`a`, `b` (associated with  `R1`, `R2`, respectively) and specified 
interaxial angle `γ=angle(R1,R2)`.

For definiteness, the `R1` basis vector is oriented along the
x-axis of the Cartesian coordinate system.
"""
function crystal(a::Real,b::Real,γ::Real) 
    R1 = Float64[a, 0] 
    R2 = b.*[cos(γ), sin(γ)]

    return Crystal((R1,R2))
end

# For a three-axis system, α, β, and γ are subject to constraints: specifically, 
# since they correspond to sides of a (unit-radius) spherical triangle, they 
# are subject to identical constraints. These constraints are
#     0 < α + β + γ < 2π,                           (1)
#     sin(s-α)*sin(s-β)*sin(s-γ)/sin(s) > 0,        (2)
# with s = (α + β + γ)/2. Constraint (2) can be identified from Eq. (38) of 
# http://mathworld.wolfram.com/SphericalTrigonometry.html; due to (1), it can 
# be simplified to sin(s-α)*sin(s-β)*sin(s-γ) > 0. This impacts generation 
# of triclinic and monoclinic crystals.
function isvalid_sphericaltriangle(α,β,γ)
    s = (α+β+γ)/2
    check1 = 0 < s < π;                     
    check2 = sin(s-α)*sin(s-β)*sin(s-γ) > 0 
    return check1 && check2 
end


const origin_markeropts = (marker="o", markerfacecolor="white", markeredgecolor="black", markeredgewidth=1.5, markersize=4.5)

function plot(C::Crystal)
    R = basis(C)
    if dim(C) == 2
        corner = sum(R)
        for R′ in R
            plot([0, R′[1]], [0, R′[2]]; color="black") # basis vectors
            plot([R′[1], corner[1]], [R′[2], corner[2]]; color="grey") # remaining unit cell boundaries
        end
        plot([0,], [0,]; origin_markeropts...) # origin
    elseif dim(C) == 3
        corners = (R[1]+R[3], R[1]+R[2], R[2]+R[3])
        dirs = ((-1,1,-1), (-1,-1,1), (1,-1,-1))
        for (i,R) in enumerate(R)
            plot3D([0, R[1]], [0, R[2]], [0, R[3]]; color="black") # basis vectors
            for (corner,dir) in zip(corners,dirs) # remaining unit cell boundaries
                plot3D([corner[1], corner[1]+dir[i]*R[1]], 
                       [corner[2], corner[2]+dir[i]*R[2]], 
                       [corner[3], corner[3]+dir[i]*R[3]]; color="grey")
            end
        end
        plot3D([0,], [0,], [0,]; origin_markeropts...) # origin
        plt.gca().set_zlabel("z")
    end
    plt.gca().set_xlabel("x"); plt.gca().set_ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box") # seems broken in 3D (https://github.com/matplotlib/matplotlib/pull/13474)
    return nothing
end

°(φ::Real) = deg2rad(φ)

""" 
    crystalsystem(C::Crystal)

Determine the crystal system of a point lattice specified in a 
*conventional* basis using Table 2.1.2.1, 9.1.7.1, & 9.1.7.2 of 
the International Tables of Crystallography, Volume 1 (ITA). 
There are 4 crystal systems in 2D and 7 in 3D (see ITA 2.1.2(iii)):

      |_DIM_|_SYSTEM_______|_CONDITIONS_____________|_FREE PARAMS___|
      | 2D  | square       | a=b & γ=90°            | a             |
      |     | rectangular  | γ=90°                  | a,b           |
      |     | hexagonal    | a=b & γ=120°           | a             |
      |     | oblique      | none                   | a,b,γ         |
      |-----|--------------|------------------------|---------------|
      | 3D  | cubic        | a=b=c & α=β=γ=90°      | a             |
      |     | hexagonal    | a=b & α=β=90° & γ=120° | a,c           |
      |     | trigonal     | a=b=c & α=β=γ          | a,α or a,c    |
      |     | tetragonal   | a=b & α=β=γ=90°        | a,c           |
      |     | orthorhombic | α=β=γ=90°              | a,b,c         |
      |     | monoclinic   | α=γ=90°                | a,b,c,β≥90°   |
      |     | triclinic    | none                   | a,b,c,α,β,γ   |

The Crystal input is assumed to use *conventional* basis vectors; 
i.e. not necessarily primitive. For primitive basis vectors, the 
crystal system can be further reduced into 5 Bravais types in 2D and
14 in 3D.
"""
function crystalsystem(C::Crystal)
    if dim(C) == 2
        a,b = norms(C)
        γ = angles(C)
        if a≈b && γ≈°(90)
            system = "square"
        elseif γ≈°(90)
            system = "rectangular"
        elseif a≈b && γ≈°(120)
            system = "hexagonal"
        else
            system = "oblique"
        end

    elseif dim(C) == 3 
        a,b,c = norms(C)
        α,β,γ = angles(C)
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
    end
    return system
end


function crystalsystem(sgnum::Integer, dim::Integer=3)
    if dim == 2
        if      sgnum ∈ 1:2;   return "oblique"      # mp
        elseif  sgnum ∈ 3:9;   return "rectangular"  # op, oc
        elseif  sgnum ∈ 10:12; return "square"       # tp
        elseif  sgnum ∈ 13:17; return "hexagonal"    # hp
        else    throw(DomainError(sgnum, "There are only 17 two-dimensional plane groups."))
        end
    
    elseif dim == 3
        if      sgnum ∈ 1:2;     return "triclinic"     # aP
        elseif  sgnum ∈ 3:15;    return "monoclinic"    # mP, mC
        elseif  sgnum ∈ 16:74;   return "orthorhombic"  # oP, oI, oF, oC
        elseif  sgnum ∈ 75:142;  return "tetragonal"    # tP, tI
        elseif  sgnum ∈ 143:167; return "trigonal"      # ? hP, hI ?
        elseif  sgnum ∈ 168:194; return "hexagonal"     # hR, hP
        elseif  sgnum ∈ 195:230; return "cubic"         # cP, cI, cF
        else    throw(DomainError(sgnum, "There are only 230 three-dimensional space groups."))
        end
    end
end

"""
    relrand(lims::NTuple{2,Real}, N=1) --> Vector{Float64}

Computes a random number in the range specified by the two-element 
tuple `lims`. The random numbers are sampled from two uniform 
distributions, namely [lims[1], 1.0] and [1.0, lims[2]], in such a
way as to ensure that the sampling is uniform over the joint  
interval [-1/lims[1], -1.0] ∪ [1.0, lims[2]].

This is useful for ensuring an even sampling of numbers that are
either smaller or larger than unity. Eg. for `x=relrand((0.2,5.0))`,
`x` is equally probable to fall in inv(x)∈[1,5] or x∈[1,5].
"""
function relrand(lims::NTuple{2,Real})
    low, high = lims; invlow = inv(low)
    lowthres = (invlow - 1.0)/(invlow + high - 2.0)
    if rand() < lowthres && low < 1.0   # smaller than 1.0
        r = rand(Uniform(low,1.0))
    elseif high > 1.0                   # bigger than 1.0
        r = rand(Uniform(1.0,high))
    else                                # default
        return rand(Uniform(low,high))
    end
end
relrand(lims::NTuple{2,Real}, N) = [relrand(lims) for i=Base.OneTo(N)]

""" 
    gen_crystal(sgnum, dim=3; abclims, αβγlims) ---> Crystal{dim}

Generates a Crystal (in a conventional basis) compatible with the
space group number `sgnum`. By convention, the length of the first
lattice vector (= `a`) is set to unity, such that the second and
third (= `b` and `c`) lattice vectors' lengths are relative to the
first.

Limits on the relative uniform distribution of lengths `b` and `c`
can be specified as 2-tuple kwarg `abclims`; similarly, limits on 
the angles `α`, `β`, `γ` can be set via αβγlims (only affects 
oblique, monoclinic, & triclinic lattices).
"""
function gen_crystal(sgnum::Integer, dim=3;
                     abclims::NTuple{2,Real}=(0.5,2.0), 
                     αβγlims::NTuple{2,Real}=(°(30),°(150)))
    system = crystalsystem(sgnum, dim)
    if dim == 2
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
            γ = rand(Uniform(αβγlims...)) 
        else 
            throw(DomainError(system))
        end
        return crystal(a,b,γ)

    elseif dim == 3
        if     system == "cubic"       # a=b=c & α=β=γ=90° (free: a)
            a = b = c = 1.0
            α = β = γ = °(90)
        elseif system == "hexagonal"   # a=b & α=β=90° & γ=120° (free: a,c)
            a = b = 1.0;        c = relrand(abclims)
            α = β = °(90);      γ = °(120)
        elseif system == "trigonal"   # TODO 
            throw(DomainError(system, "The trigonal case (143-167) is not yet well thought-out"))
            # rhombohedral axes                   (a = b = c, α=β=γ < 120° ≠ 90° ?)
            # hexagonal axes, triple obverse axes (a = b ≠ c, α=β=90°, γ=120° ?)
            # maybe consult http://img.chem.ucl.ac.uk/sgp/large/sgp.htm and 
            # the setting choices in Bilbao & ISOTROPY
        elseif system == "tetragonal"  # a=b & α=β=γ=90° (free: a,c)
            a = b = 1.0;        c = relrand(abclims)
            α = β = γ = °(90)
        elseif system == "orthorhombic"# α=β=γ=90° (free: a,b,c)
            a = 1.0;            b, c = relrand(abclims, 2)
            α = β = γ = °(90)
        elseif system == "monoclinic"  # α=γ=90° (free: a,b,c,β≥90°)
            a = 1.0;            b, c = relrand(abclims, 2)
            α = γ = °(90);      β = rand(Uniform(°(90), αβγlims[2]))
            while !isvalid_sphericaltriangle(α,β,γ)  # arbitrary combinations of α,β,γ may not correspond 
                β = rand(Uniform(°(90), αβγlims[2])) # to a valid axis-system; reroll until they do
            end
        elseif system == "triclinic"   # no conditions (free: a,b,c,α,β,γ)
            a = 1.0;            b, c = relrand(abclims, 2)
            α, β, γ = rand(Uniform(αβγlims...),3)
            while !isvalid_sphericaltriangle(α,β,γ)   # arbitrary combinations of α,β,γ may not correspond 
                α, β, γ = rand(Uniform(αβγlims...),3) # to a valid axis-system; reroll until they do
            end
        else 
            throw(DomainError(system))
        end        
        return crystal(a,b,c,α,β,γ)

    else 
        _throw_invaliddim(dim)
    end
end

const CRYSTALSYSTEM_ABBREV_3D = Dict("triclinic"=>'a', "monoclinic"=>'m', "orthorhombic"=>'o', 
                                     "tetragonal"=>'t', "trigonal"=>'h', "hexagonal"=>'h', 
                                     "cubic"=>'c')
const CRYSTALSYSTEM_ABBREV_2D = Dict("oblique"=>'m', "rectangular"=>'o', "square"=>'t', 
                                     "hexagonal"=>'h')

function bravaistype(sgnum::Integer, dim::Integer=3)
    cntr = centering(sgnum, dim)
    system = crystalsystem(sgnum, dim)
    if dim == 3      # pick the correct abbreviation from a Dict
        return CRYSTALSYSTEM_ABBREV_3D[system]*cntr
    elseif dim == 2
        return CRYSTALSYSTEM_ABBREV_2D[system]*cntr
    end
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
const primitivematrix_3D = Dict(
         'P'=>[1 0 0; 0 1 0; 0 0 1],
         'F'=>[0 1 1; 1 0 1; 1 1 0]./2,
         'I'=>[-1 1 1; 1 -1 1; 1 1 -1]./2,
         'R'=>[2 -1 -1; 1 1 -2; 1 1 1]./3,
         'A'=>[2 0 0; 0 1 -1; 0 1 1]./2,
         'C'=>[1 1 0; -1 1 0; 0 0 2]./2
         )
const primitivematrix_2D = Dict(
         'c'=>[1 1; -1 1]./2, 
         'p'=>[1 0; 0 1])

"""
    primitivebasismatrix(cntr::Char, dim::Integer) -> ::matrix

Calculates a transformation matrix `P` from a conventional
to a primitive unit cell, using dictionary lookup.
"""
function primitivebasismatrix(cntr::Char, dim::Integer=3)
    if dim == 3
        return primitivematrix_3D[cntr];         
    elseif dim == 2
        return primitivematrix_2D[cntr];
    else
        _throw_invaliddim(dim)
    end
end

function centeringtranslation(cntr::Char, dim::Integer=3)
    if dim == 3      # pick the correct abbreviation from a Dict
        if cntr == 'P';     return zeros(Float64,3)
        elseif cntr == 'I'; return [1,1,1]/2
        elseif cntr == 'F'; return [1,0,1]/2
        elseif cntr == 'R'; return [2,1,1]/3
        elseif cntr == 'A'; return [0,1,1]/2
        elseif cntr == 'C'; return [1,1,0]/2
        else;               _throw_invalidcntr(cntr)
        end
    elseif dim == 2
        if cntr == 'p';     return zeros(Float64,2)
        elseif cntr == 'c'; return [1,1]/2
        else;               _throw_invalidcntr(cntr)
        end
    elseif dim == 1
        return zeros(Float64, 1)
    else 
        _throw_invaliddim(dim)
    end
end
@noinline _throw_invalidcntr(cntr::Char) = throw(DomainError(cntr, "input centering character must be {P,I,F,R,A,C} in 3D, {p,c} in 2D, or p in 1D"))
@noinline _throw_invaliddim(dim::Integer) = throw(DomainError(dim, "input dimension must be 1, 2, or 3"))
""" 
    primitivebasis(sgnum::Integer, C::Crystal) --> Cp::Crystal

Transforms the conventional basis of a Crystal `C` into its primitive 
equivalent `Cp`, provided that its centering differs from the conventional
(P or p), by inferring the Bravais type from the space group number
`sgnum` and applying an applying an appropriate transformation matrix. 
"""
function primitivebasis(sgnum::Integer, C::Crystal)
    cntr = centering(sgnum)
    return primitivebasis(C, cntr)
end

"""
    primitivebasis(C::Crystal, cntr::Char) --> Cp::Crystal

Transforms the conventional basis of a Crystal `C` into its primitive 
equivalent `Cp`, with the transformation dependent on the centering
type `cntr` (P, I, F, R, A, C, and p, c); for centering P and p, the 
conventional and primive bases coincide.
"""
function primitivebasis(C::Crystal, cntr::Char)
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return C
    else         
        P = primitivebasismatrix(cntr, dim(C))
        R_P = hcat(basis(C)...)*P # R_P = R_C*P (w/ R_C a matrix w/ columns of conventional direct basis vecs)
        newbasis = Tuple(collect(u) for u in eachcol(R_P)) # convert from matrix form back to tuple form
        return Crystal(newbasis)
    end  
end


"""
    reciprocalbasis(C::Crystal) --> G::NTuple{dim(C), Vector{Float64}}
    
Calculates the reciprocal basis vectors associated with a crystal `C`.
"""
reciprocalbasis(C::Crystal) = reciprocalbasis(basis(C))
function reciprocalbasis(R::NTuple{N, Vector{<:Real}}) where N
    if N == 3
        pref = 2π/dot(R[1], (R[2]×R[3]))
        return pref .* (R[2]×R[3], R[3]×R[1], R[1]×R[2])
    elseif N == 2
        pref = 2π/dot(R[1], [-R[2][2], R[2][1]])
        return pref .* ([-R[2][2], R[2][1]], [R[1][2], -R[1][1]])
    elseif N == 1
        return (2π/first(R[1]),)
    else
        # the general definition of the reciprocal basis is 
        # [G₁ G₂ ... Gₙ]ᵀ = 2π[R₁ R₂ ... Rₙ]⁻¹; that form is
        # a bit slower than the above specific variants, 
        # however, cf. the inversion operation, so we only 
        # use it as a hig-dimensional fallback (i.e. breadcrumbs)
        return tuple(eachrow((2π*I/hcat(R...)))...) 
    end
end

"""
    primitivereciprocalbasis(C::Crystal, cntr::Char) --> G_P::NTuple{dim(C), Vector{Float64}}
    
Calculates the **primitive** reciprocal basis vectors associated with a 
crystal `C` of centering type `cntr`.
"""
primitivereciprocalbasis(C::Crystal, cntr::Char) = primitivereciprocalbasis(basis(C), cntr::Char)
function primitivereciprocalbasis(R::NTuple{N, Vector{<:Real}}, cntr::Char) where N
    G_C = reciprocalbasis(R)
    
    # While the direct basis (𝐚 𝐛 𝐜) transforms like 
    #       (𝐚′ 𝐛′ 𝐜′) = (𝐚 𝐛 𝐜)𝐏
    # under a basis change matrix 𝐏, the direct basis
    # (𝐚* 𝐛* 𝐜*) transforms like 
    #       (𝐚*′ 𝐛*′ 𝐜*′) = (𝐚* 𝐛* 𝐜*)(𝐏⁻¹)ᵀ
    # since (𝐚 𝐛 𝐜)(𝐚* 𝐛* 𝐜*)ᵀ = 2π𝐈 must be conserved
    # after the basis change
    P = primitivebasismatrix(cntr, N)
    G_P = hcat(G_C...)/P' # G_P = G_C*(P⁻¹)ᵀ = G_C*(Pᵀ)⁻¹
                          # (w/ G_C a matrix w/ columns conventional reciprocal vecs)
    
    return Tuple(collect(v) for v in eachcol(G_P))    
end
# Note that the 𝑐𝑜𝑒𝑓𝑓𝑒𝑐𝑖𝑒𝑛𝑡𝑠 of a general 𝐤-vector transform
# differently than the reciprocal basis, which transforms
# from non-primed to primed variants via a basis matrix 𝐏
# according to (see also `primitivereciprocalbasis(...)`):
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

