# Transformation matrices P that map bases and coordinate vectors in either direct or
# reciprocal space from a conventional to a primitive setting:
#
# === Basis transformation ===
# --- Direct basis ---
# Under P, a direct-space conventional basis (𝐚 𝐛 𝐜) is mapped to a primitive basis
#   (𝐚′ 𝐛′ 𝐜′) = (𝐚 𝐛 𝐜)𝐏
# --- Reciprocal basis ---
# Under P, a reciprocal-space conventional basis (𝐚* 𝐛* 𝐜*) is mapped to a primitive basis
#   (𝐚*′ 𝐛*′ 𝐜*′) = (𝐚* 𝐛* 𝐜*)(𝐏⁻¹)ᵀ
# since (𝐚 𝐛 𝐜)(𝐚* 𝐛* 𝐜*)ᵀ = 2πI must be conserved after the basis change.
#
# === Coordinate vector transformation ===
# The _coefficients_ of a vector transform differently than its _bases_. Specifically:
# --- Direct coordinate vectors ---
# An 𝐫-vector specified in a conventional reciprocal basis (𝐚 𝐛 𝐜) with an associated
# coefficient vector (r₁ r₂ r₃)ᵀ, i.e. 𝐫 ≡ (𝐚 𝐛 𝐜)(r₁ r₂ r₃)ᵀ [w/ (r₁ r₂ r₃)ᵀ a column
# vector], is mapped to a primitive coefficient vector under P:
#     (r₁′ r₂′ r₃′)ᵀ = P⁻¹(r₁ r₂ r₃)ᵀ
# since
#     𝐤 = (𝐚′ 𝐛′ 𝐜′)(r₁′ r₂′ r₃′)ᵀ     (1)  [... by definition]
#       = (𝐚 𝐛 𝐜)P(r₁′ r₂′ r₃′)ᵀ            [... transformation of (𝐚 𝐛 𝐜) under P]
#       = (𝐚 𝐛 𝐜)(r₁ r₂ r₃)ᵀ           (2)  [... by definition]
# then, combining (1) and (2)
#     P(r₁′ r₂′ r₃′)ᵀ = (r₁ r₂ r₃)ᵀ
#  ⇔ (r₁′ r₂′ r₃′)ᵀ = P⁻¹(r₁ r₂ r₃)ᵀ
# --- Reciprocal coordinate vectors ---
# A 𝐤-vector specified in a conventional reciprocal basis (𝐚* 𝐛* 𝐜*) with an associated
# coefficient vector (k₁ k₂ k₃)ᵀ, i.e. 𝐤 ≡ (𝐚* 𝐛* 𝐜*)(k₁ k₂ k₃)ᵀ [w/ (k₁ k₂ k₃)ᵀ a column
# vector], is mapped to a primitive coefficient vector under P
#     (k₁′ k₂′ k₃′)ᵀ = Pᵀ(k₁ k₂ k₃)ᵀ
# since
#     𝐤 = (𝐚*′ 𝐛*′ 𝐜*′)(k₁′ k₂′ k₃′)ᵀ     (1)  [... by definition]
#       = (𝐚* 𝐛* 𝐜*)(P⁻¹)ᵀ(k₁′ k₂′ k₃′)ᵀ       [... transformation of (𝐚* 𝐛* 𝐜*) under P]
#       = (𝐚* 𝐛* 𝐜*)(k₁ k₂ k₃)ᵀ           (2)  [... by definition]
# then, combining (1) and (2)
#     (P⁻¹)ᵀ(k₁′ k₂′ k₃′)ᵀ = (k₁ k₂ k₃)ᵀ
#  ⇔ (k₁′ k₂′ k₃′)ᵀ = Pᵀ(k₁ k₂ k₃)ᵀ
#
# The values of P depend on convention. We adopt those of Table 2 of the Aroyo's Bilbao
# publication (https://doi.org/10.1107/S205327331303091X), which give the coefficients of
# (Pᵀ)⁻¹.
# See also the 2016 HPKOT/Hinuma paper (https://doi.org/10.1016/j.commatsci.2016.10.015)
# for additional details and context, though note that they use different matrices for 'A'
# and complicate the 'C' scenario (Table 3).
# Note that, by convention, the centering type 'B' never occurs among the space groups.

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

The choice of transformation matrix `P` (equivalently, the assumed setting choice) is
consistent with the choice in the International Table of Crystallography and the Bilbao
Crystallographic Server [^1].

[^1]: Aroyo et al., [Acta Cryst. A70, 126 (2014)](https://doi.org/10.1107/S205327331303091X):
      Table 2 gives (`P`ᵀ)⁻¹.
"""
@inline function primitivebasismatrix(cntr::Char, ::Val{D}=Val(3)) where D
    D∉1:3 && _throw_invaliddim(D)
    return PRIMITIVE_BASIS_MATRICES[D][cntr]
end

@inline function centeringtranslation(cntr::Char, ::Val{D}=Val(3)) where D
    if D == 3
        if     cntr == 'P'; return zeros(SVector{3})
        elseif cntr == 'I'; return SVector((1,1,1)./2)
        elseif cntr == 'F'; return SVector((1,0,1)./2)
        elseif cntr == 'R'; return SVector((2,1,1)./3)
        elseif cntr == 'A'; return SVector((0,1,1)./2)
        elseif cntr == 'C'; return SVector((1,1,0)./2)
        else;               _throw_invalidcntr(cntr)
        end
    elseif D == 2
        if     cntr == 'p'; return zeros(SVector{2})
        elseif cntr == 'c'; return SVector((1,1)./2)
        else;               _throw_invalidcntr(cntr)
        end
    elseif D == 1
        return zeros(SVector{1})
    else 
        _throw_invaliddim(D)
    end
end

function all_centeringtranslations(cntr::Char, Dᵛ::Val{D}=Val(3)) where D
    if D == 3 && cntr == 'F'
        # primitive cell has 1/4th the volume of conventional cell: 3 extra centers
        return [SVector((1,0,1)./2), SVector((0,1,1)./2), SVector((1,1,0)./2)]
    elseif D == 3 && cntr == 'R'
        # primitive cell has 1/3rd the volume of conventional cell: 2 extra centers
        return [SVector((2,1,1)./3), SVector((1,2,2)./3)]
    else
        # primitive cell has half the volume of conventional cell: 1 extra center
        return [centeringtranslation(cntr, Dᵛ)]
    end
end


"""
    reciprocalbasis(Rs)  -->  Gs::ReciprocalBasis{D}
    
Return the reciprocal basis `Gs` of a direct basis `Rs` in `D` dimensions, provided as a
`DirectBasis{D}` or a `D`-dimensional `NTuple` or `StaticVector` of `AbstractVector`s.
"""
function reciprocalbasis(Rs::Union{DirectBasis{D}, 
                                   NTuple{D, <:AbstractVector{<:Real}},
                                   StaticVector{D, <:AbstractVector{<:Real}}}) where D
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
        Rm = stack(Rs)
        Gm = 2π.*inv(transpose(Rm))
        vecs = ntuple(i->Gm[:,i], Val(D))
    end

    return ReciprocalBasis{D}(vecs)
end


""" 
    primitivize(Vs::AbstractBasis, sgnum::Integer) --> Vs′::AbstractBasis

Return the primtive `AbstractBasis` `Vs′` associated with the input (a conventional basis
`Vs` in dimension `D` in a crystal system consistent with space group number `sgnum`).
The space group number is used to infer the associated centering type which determines the
required transformation (see also [`centering`](@ref)).

For centering types `'P'` and `'p'`, the conventional and primitive bases coincide.
"""
function primitivize(Vs::AbstractBasis{D}, sgnum::Integer) where D
    cntr = centering(sgnum, D)
    return primitivize(Vs, cntr)
end

function transform(Rs::DirectBasis{D}, P::AbstractMatrix{<:Real}) where D
    # Rm′ = Rm*P (w/ Rm a matrix w/ columns of untransformed direct basis vecs Rᵢ)
    Rm′ = stack(Rs)*P
    return DirectBasis{D}(ntuple(i->Rm′[:,i], Val(D)))
end

"""
    primitivize(Rs::DirectBasis, cntr::Char) --> Rs′::DirectBasis

Return the primtive direct basis `Rs′` associated with the input (a conventional direct
basis `Rs` with centering type `cntr`).
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

Return the conventional direct basis `Rs` associated with the input (a primitive direct
basis `Rs′` with centering type `cntr`).
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
    # Gm′ = Gm*(P⁻¹)ᵀ = Gm*(Pᵀ)⁻¹ (w/ Gm a matrix w/ columns of untransformed reciprocal
    # vecs Gᵢ)
    Gm′ = stack(Gs)/P'
    return ReciprocalBasis{D}(ntuple(i->Gm′[:,i], Val(D)))
end

"""
    primitivize(Gs::ReciprocalBasis, cntr::Char) --> Gs′::ReciprocalBasis
    
Return the primitive reciprocal basis `Gs′` associated with the input (a conventional
reciprocal basis `Gs` with centering type `cntr`).
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