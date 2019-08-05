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
    # R1 and R2 are easy
    R1 = Float64[a, 0, 0] 
    R2 = b.*[cos(γ), sin(γ), 0]
    # R3 is harder
    cosα = cos(α)
    cosβ = cos(β)
    sinγ,cosγ = sincos(γ)
    ϕ = atan(cosα - cosγ*cosβ, sinγ*cosβ)
    # TODO: There is some issue here in triclinic case: it can happen 
    # that cosβ < cosϕ. This logic needs more careful consideration
    # to ensure the appropriate domains of ϕ∈[0,2π] and θ∈[0,π]
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
        else    DomainError(sgnum, "There are only 17 two-dimensional plane groups.")
        end
    
    elseif dim == 3
        if      sgnum ∈ 1:2;     return "triclinic"     # aP
        elseif  sgnum ∈ 3:15;    return "monoclinic"    # mP, mC
        elseif  sgnum ∈ 16:74;   return "orthorhombic"  # oP, oI, oF, oC
        elseif  sgnum ∈ 75:142;  return "tetragonal"    # tP, tI
        elseif  sgnum ∈ 143:167; return "trigonal"      # ? hP, hI ?
        elseif  sgnum ∈ 168:194; return "hexagonal"     # hR, hP
        elseif  sgnum ∈ 195:230; return "cubic"         # cP, cI, cF
        else    DomainError(sgnum, "There are only 230 three-dimensional space groups.")
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
                     abclims::NTuple{2,Real}=(0.2,5.0), 
                     αβγlims::NTuple{2,Real}=(°(10),°(170)))
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
            error(DomainError(system))
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
            error("The trigonal case is not yet well thought-out")
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
        elseif system == "triclinic"   # no conditions (free: a,b,c,α,β,γ)
            a = 1.0;            b, c = relrand(abclims, 2)
            α, β, γ = rand(Uniform(αβγlims...),3)
            #if abs(α)+abs(β)<abs(γ) # cannot be satisfied; reroll
            #    return gen_crystal(sgnum, dim; abclims=abclims, αβγlims=αβγlims)
            #end
        else 
            error(DomainError(system))
        end        
        return crystal(a,b,c,α,β,γ)

    else 
        error(ArgumentError(dim, "dimension must be 2 or 3"))
    end
end

const crystalsystem_abbrev_3D = Dict("triclinic"=>"a", "monoclinic"=>"m", "orthorhombic"=>"o", 
                                     "tetragonal"=>"t", "trigonal"=>"h", "hexagonal"=>"h", 
                                     "cubic"=>"c")
const crystalsystem_abbrev_2D = Dict("oblique"=>"m", "rectangular"=>"o", "square"=>"t", 
                                     "hexagonal"=>"h")

function bravaistype(sgnum::Integer, dim::Integer=3)
    cntr = centering(sgnum, dim)
    system = crystalsystem(sgnum, dim)
    if cntr != 'R'
        if dim == 3      # pick the correct abbreviation from a Dict
            bravaistype = crystalsystem_abbrev_3D[system]*cntr
        elseif dim == 2
            bravaistype = crystalsystem_abbrev_2D[system]*cntr
        end
    else # TODO: special rules for trigonal systems that in the R (= rhombohedral) subset ?
        bravaistype = "rR"
    end
    return bravaistype
end




const primitivematrix_3D = Dict(
         "P"=>[1 0 0; 0 1 0; 0 0 1],        # transformation matrices P from conventional basis v_C 
         "I"=>[-1 1 1; 1 -1 1; 1 1 -1]./2,  # to primitive basis v_p, depending on centering types
         "F"=>[0 1 1; 1 0 1; 1 1 0]./2,     #     v_P = v_C*P
         "R"=>[-1 2 -1; -2 1 1; 1 1 1]./3,  # with v_P and v_C interpreted as matrices = [R1 R2 R3]
         "A"=>[2 0 0; 0 1 -1; 0 1 1]./2,
         "C"=>[1 1 0; -1 1 0; 0 0 2]./2)    # "B"=>[], # seems to not occur, by convention
const primitivematrix_2D = Dict("c"=>[1 1; -1 1]./2, "p"=>[1 0; 0 1])
"""
    primitivebasismatrix(cntr::Union{String, Char}, dim::Integer) -> ::matrix

    Calculates a transformation matrix `P` from a conventional
    to a primitive unit cell, using dictionary lookup.
"""
function primitivebasismatrix(cntr::Union{String, Char}, dim)
    if dim == 3
        return primitivematrix_3D[string(cntr)];         
    elseif dim == 2
        return primitivematrix_2D[string(cntr)];
    else
        error(ArgumentError("dim must be either 2 or 3"))
    end
end


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
    primitivebasis(C::Crystal, cntr::String) --> Cp::Crystal

    Transforms the conventional basis of a Crystal `C` into its primitive 
    equivalent `Cp`, with the transformation dependent on the centering
    type `cntr` (P, I, F, R, A, C, and p, c); for centering P and p, the 
    conventional and primive bases coincide.
"""
function primitivebasis(C::Crystal, cntr::String)
    if cntr == "P" || cntr == "p" # the conventional and primitive bases coincide
        return C
    else         
        P = primitivebasismatrix(cntr, dim(C))
        v_P = hcat(basis(C)...)*P # v_P = v_C*P (with v_C a matrix with columns = basis vecs)
        newbasis = Tuple(collect(u) for u in eachcol(v_P)) # convert from matrix form back to tuple form
        return Crystal(newbasis)
    end  
end


"""
    reciprocalbasis(C::Crystal) --> G::NTuple{dim(C), Vector{Float64}}
    
    Calculates the reciprocal basis vectors associated with a crystal C
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

# --- TESTING ---
C = crystal(1,1,1, π/3, π/3, π/3)
#display(B)
plt.close("all")
plot(C)
crystalsystem(C)

