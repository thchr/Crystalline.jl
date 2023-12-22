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
# (Pᵀ)⁻¹. Equivalently, this is the "CDML" setting choices that can also be inferred by
# combining Tables 1.5.4.1 and 1.5.4.2 of the International Tables of Crystallography,
# Vol. B, Ed. 2, 2001 (ITB2). Of note, this is _not_ the standard ITA choice for the
# primitive cell for 'R' or 'A'-centered cells (see Tables 3.1.2.2 of ITA6; the CDML
# convention is more widespread, however, especially for k-vectors; hence our choice).
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

function canonicalize_centering(cntr, ::Val{D}, ::Val{P}) where {D,P}
    if D == P # space/plane/line groups
        return cntr
    elseif D == 3 && P == 2 # layer groups
        return cntr == '𝑝' ? 'p' : 
               cntr == '𝑐' ? 'c' : error(DomainError(cntr, "invalid layer group centering"))
    elseif D == 3 && P == 1 # rod groups
        return cntr == '𝓅' ? 'p' : error(DomainError(cntr, "invalid rod group centering"))
    elseif D == 2 && P == 1 # frieze groups
        return cntr == '𝓅' ? 'p' : error(DomainError(cntr, "invalid frieze group centering"))
    else
        throw(DomainError((D,P), "invalid combination of dimensionality D and periodicity P"))
    end
end

@doc raw"""
    primitivebasismatrix(cntr::Char, ::Val{D}=Val(3)) --> SMatrix{D,D,Float64}
    primitivebasismatrix(cntr::Char, ::Val{D}, ::Val{P}) --> SMatrix{D,D,Float64}

Return the transformation matrix ``\mathbf{P}`` that transforms a conventional unit cell
with centering `cntr` to the corresponding primitive unit cell (in dimension `D` and
periodicity `P`) in CDML setting.

If `P` is not provided, it default to `D` (as e.g., applicable to Crystalline.jl's 
`spacegroup`). If `D` and `P` differ, a subperiodic group setting is assumed (as e.g.,
applicable to Crystalline.jl's `subperiodicgroup`).

## Transformations in direct and reciprocal space

### Bases

The returned transformation matrix ``\mathbf{P}`` transforms a direct-space conventional
basis ``(\mathbf{a}\ \mathbf{b}\ \mathbf{c})`` to the direct-space primitive basis

```math
    (\mathbf{a}'\ \mathbf{b}'\ \mathbf{c}') =
    (\mathbf{a}\ \mathbf{b}\ \mathbf{c})\mathbf{P}.
```

Analogously, ``\mathbf{P}`` transforms a reciprocal-space conventional basis
``(\mathbf{a}^*\ \mathbf{b}^*\ \mathbf{c}^*)`` to 

```math
(\mathbf{a}^{*\prime}\ \mathbf{b}^{*\prime}\ \mathbf{c}^{*\prime}) =
(\mathbf{a}^*\ \mathbf{b}^*\ \mathbf{c}^*)(\mathbf{P}^{-1})^{\text{T}}.
```

see also [`transform(::DirectBasis, ::AbstractMatrix{<:Real})`](@ref) and
[`transform(::ReciprocalBasis, ::AbstractMatrix{<:Real})`](@ref)).

### Coordinates

The coordinates of a point in either direct or reciprocal space, each referred to a basis,
also transform under ``\mathbf{P}``. Concretely, direct- and reciprocal-space
conventional points ``\mathbf{r} = (r_1, r_2, r_3)^{\text{T}}`` and
``\mathbf{k} = (k_1, k_2, k_3)^{\text{T}}``, respectively, transform to a primitive setting
under ``\mathbf{P}`` according to:

```math
\mathbf{r}' = \mathbf{P}^{-1}\mathbf{r},\\
\mathbf{k}' = \mathbf{P}^{\text{T}}\mathbf{k}.
```

See also [`transform(::DirectPoint, ::AbstractMatrix{<:Real})`](@ref) and
[`transform(::ReciprocalPoint, ::AbstractMatrix{<:Real})`](@ref)).

## Setting conventions

The setting choice for the primitive cell implied by the returned ``\mathbf{P}`` follows the
widely adopted Cracknell-Davies-Miller-Love (CDML) convention.[^CDML]
This convention is explicated e.g. in Table 2 of [^Aroyo] (or, alternatively, can be
inferred from Tables 1.5.4.1 and 1.5.4.2 of [^ITB2]) and is followed e.g. on the Bilbao
Crystallographic Server[^BCS], in the CDML reference work on space group irreps[^CDML], and
in the C library `spglib`.[^spglib]

Note that this setting choice is _not_ what is frequently referred to as the "ITA primitive
setting", from which it differs for hP, hR, and oA Bravais types.

The setting choice is usually referred to as the CDML primitive setting, or, less
frequently and more ambiguously, as the crystallographic primitive setting.

[^CDML]: Cracknell, Davies, Miller, & Love, Kroenecker Product Tables, Vol. 1 (1979).

[^BCS]: Bilbao Crystallographic Server, [KVEC](https://www.cryst.ehu.es/cryst/get_kvec.html).

[^Aroyo]: Aroyo *et al.*, [Acta Cryst. A70, 126 (2014)]
          (https://doi.org/10.1107/S205327331303091X): Table 2 gives
          ``(\mathbf{P}^{-1})^{\text{T}}``.

[^ITB2]: Hahn, International Tables of Crystallography, Vol. B, 2nd edition (2001).

[^spglib]: Spglib documentation: [Transformation to the primitive setting]
           (https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell).
           Thus, Bravais.jl and [Spglib.jl](https://github.com/singularitti/Spglib.jl)
           transform to identical primitive settings and are hence mutually compatible.
"""
@inline function primitivebasismatrix(cntr::Char, ::Val{D}=Val(3), ::Val{D}=Val(D)) where D
    # space groups
    D ∉ 1:3 && _throw_invalid_dim(D)
    return PRIMITIVE_BASIS_MATRICES[D][cntr]
end
@inline function primitivebasismatrix(cntr::Char, Dᵛ::Val{3}, Pᵛ::Val{2})
    # layer and rod groups
    cntr = canonicalize_centering(cntr, Dᵛ, Pᵛ)
    P²ᴰ = PRIMITIVE_BASIS_MATRICES[2][cntr]
    return @SMatrix [P²ᴰ[1,1] P²ᴰ[2,1] 0.0; P²ᴰ[1,2] P²ᴰ[2,2] 0.0; 0.0 0.0 1.0]
end
@inline function primitivebasismatrix(cntr::Char, Dᵛ::Val{3}, Pᵛ::Val{1})
    # rod groups
    cntr = canonicalize_centering(cntr, Dᵛ, Pᵛ)
    P²ᴰ = PRIMITIVE_BASIS_MATRICES[1][cntr]
    return @SMatrix [P²ᴰ[1,1] 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]
end
@inline function primitivebasismatrix(cntr::Char, Dᵛ::Val{2}, Pᵛ::Val{1})
    # frieze groups
    cntr = canonicalize_centering(cntr, Dᵛ, Pᵛ)
    P¹ᴰ = PRIMITIVE_BASIS_MATRICES[1][cntr]
    return @SMatrix [P¹ᴰ[1,1] 0.0; 0.0 1.0]
end
function primitivebasismatrix(cntr::Char, ::Val{D}, ::Val{P}) where {D,P}
    # fall-back error
    throw(DomainError((D,P), "invalid combination of dimensionality D and periodicity P"))
end

@inline function centeringtranslation(cntr::Char,
                                      Dᵛ::Val{D}=Val(3), Pᵛ::Val{P}=Val(D)) where {D,P}
    D ∉ 1:3 && _throw_invalid_dim(D)
    P ∉ 1:D && _throw_invalid_dim(P)
    cntr = canonicalize_centering(cntr, Dᵛ, Pᵛ)
    if D == 3
        if     cntr == 'P'; return zeros(SVector{3})
        elseif cntr == 'I'; return SVector((1,1,1)./2)
        elseif cntr == 'F'; return SVector((1,0,1)./2)
        elseif cntr == 'R'; return SVector((2,1,1)./3)
        elseif cntr == 'A'; return SVector((0,1,1)./2)
        elseif cntr == 'C'; return SVector((1,1,0)./2)
        else;               _throw_invalid_cntr(cntr, 3)
        end
    elseif D == 2
        if     cntr == 'p'; return zeros(SVector{2})
        elseif cntr == 'c'; return SVector((1,1)./2)
        else;               _throw_invalid_cntr(cntr, 2)
        end
    elseif D == 1
        if     cntr == 'p'; return zeros(SVector{1})
        else;               _throw_invalid_cntr(cntr, 1)
        end
    else 
        _throw_invalid_dim(D)
    end
end

function all_centeringtranslations(cntr::Char,
                                   Dᵛ::Val{D}=Val(3), Pᵛ::Val{P}=Val(D)) where {D,P}
    D ∉ 1:3 && _throw_invalid_dim(D)
    P ∉ 1:D && _throw_invalid_dim(P)
    cntr = canonicalize_centering(cntr, Dᵛ, Pᵛ)
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

# ---------------------------------------------------------------------------------------- #

"""
    reciprocalbasis(Rs)  -->  ::ReciprocalBasis{D}
    
Return the reciprocal basis of a direct basis `Rs` in `D` dimensions, provided as a
`StaticVector` of `AbstractVector`s (e.g., a `DirectBasis{D}`) or a `D`-dimensional `NTuple`
of `AbstractVector`s, or a (or, type-unstably, as any iterable of `AbstractVector`s).
"""
function reciprocalbasis(Rs::Union{NTuple{D, <:AbstractVector{<:Real}},
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
reciprocalbasis(Rs) = reciprocalbasis(Tuple(Rs)) # type-unstable convenience accesor

# TODO: Provide a utility to go from ReciprocalBasis -> DirectBasis. Maybe deprecate
#       `reciprocalbasis` and have a more general function `dualbasis` instead?

# ---------------------------------------------------------------------------------------- #

@doc raw"""
    transform(Rs::DirectBasis, P::AbstractMatrix{<:Real})

Transform a direct basis `Rs` ``= (\mathbf{a}\ \mathbf{b}\ \mathbf{c})`` under the
transformation matrix `P` ``= \mathbf{P}``, returning 
`Rs′` ``= (\mathbf{a}'\ \mathbf{b}'\ \mathbf{c}') = (\mathbf{a}\ \mathbf{b}\ \mathbf{c})
\mathbf{P}``.
"""
function transform(Rs::DirectBasis{D}, P::AbstractMatrix{<:Real}) where D
    # Rm′ = Rm*P (w/ Rm a matrix w/ columns of untransformed direct basis vecs Rᵢ)
    Rm′ = stack(Rs)*P
    return DirectBasis{D}(ntuple(i->Rm′[:,i], Val(D)))
end

@doc raw"""
    transform(Gs::ReciprocalBasis, P::AbstractMatrix{<:Real})

Transform a reciprocal basis `Gs` ``= (\mathbf{a}^* \mathbf{b}^* \mathbf{c}^*)`` under the
transformation matrix `P` ``= \mathbf{P}``, returning 
`Gs′` ``= (\mathbf{a}^{*\prime}\ \mathbf{b}^{*\prime}\ \mathbf{c}^{*\prime}) =
(\mathbf{a}^*\ \mathbf{b}^*\ \mathbf{c}^*)(\mathbf{P}^{-1})^{\text{T}}``.
"""
function transform(Gs::ReciprocalBasis{D}, P::AbstractMatrix{<:Real}) where D
    # Gm′ = Gm*(P⁻¹)ᵀ = Gm*(Pᵀ)⁻¹ (w/ Gm a matrix w/ columns of reciprocal vecs Gᵢ)
    Gm′ = stack(Gs)/P'
    return ReciprocalBasis{D}(ntuple(i->Gm′[:,i], Val(D)))
end

@doc raw"""
    transform(r::DirectPoint, P::AbstractMatrix{<:Real})  -->  r′::typeof(r)

Transform a point in direct space `r` ``= (r_1, r_2, r_3)^{\text{T}}`` under the
transformation matrix `P` ``= \mathbf{P}``, returning `r′` ``= (r_1', r_2', r_3')^{\text{T}}
= \mathbf{P}^{-1}(r_1', r_2', r_3')^{\text{T}}``.
"""
function transform(r::DirectPoint{D}, P::AbstractMatrix{<:Real}) where D
    return DirectPoint{D}(P\parent(r))
end

@doc raw"""
    transform(k::ReciprocalPoint, P::AbstractMatrix{<:Real})  -->  k′::typeof(k)

Transform a point in reciprocal space `k` ``= (k_1, k_2, k_3)^{\text{T}}`` under the
transformation matrix `P` ``= \mathbf{P}``, returning `k′` ``= (k_1', k_2', k_3')^{\text{T}}
= \mathbf{P}^{\text{T}}(k_1', k_2', k_3')^{\text{T}}``.
"""
function transform(k::ReciprocalPoint{D}, P::AbstractMatrix{<:Real}) where D
    return ReciprocalPoint{D}(P'*parent(k))
end

# ---------------------------------------------------------------------------------------- #

function primitivize(V::Union{AbstractBasis{D}, AbstractPoint{D}}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return V
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(V, P)
    end
end

function primitivize(V::Union{AbstractBasis{D}, AbstractPoint{D}}, sgnum::Integer) where D
    cntr = centering(sgnum, D)
    return primitivize(V, cntr)
end

@doc """ 
    primitivize(V::Union{AbstractBasis, AbstractPoint}, 
                cntr_or_sgnum::Union{Char, <:Integer})   -->  V′::typeof(V)

Return the primitive basis or point `V′` associated with the input conventional
`AbstractBasis` or `AbstractPoint` `V`.

The assumed centering type is specified by `cntr_or_sgnum`, given either as a centering
character (`::Char`) or inferred from a space group number (`::Integer`) and the
dimensionality of `V` (see also [`centering(::Integer, ::Integer)`](@ref)).
"""
primitivize(::Union{AbstractBasis, AbstractPoint}, ::Union{Char, <:Integer})

@doc """
    primitivize(Rs::DirectBasis, cntr_or_sgnum::Union{Char, <:Integer}) --> Rs′::typeof(Rs)

Return the primitive direct basis `Rs′` corresponding to the input conventional direct basis
`Rs`.
"""
primitivize(::DirectBasis, ::Union{Char, <:Integer})

@doc """
    primitivize(Gs::ReciprocalBasis, cntr_or_sgnum::Union{Char, <:Integer}) --> Gs′::typeof(Gs)
    
Return the primitive reciprocal basis `Gs′` corresponding to the input conventional
reciprocal basis `Gs`.
"""
primitivize(::ReciprocalBasis, ::Union{Char, <:Integer})

@doc """
    primitivize(r::DirectPoint, cntr_or_sgnum::Union{Char, <:Integer}) --> r′::typeof(r)

Return the direct point `r′` with coordinates in a primitive basis, corresponding to
the input point `r` with coordinates in a conventional basis. 
"""
primitivize(::DirectPoint, ::Union{Char, <:Integer})

@doc """
    primitivize(k::ReciprocalPoint, cntr_or_sgnum::Union{Char, <:Integer}) --> k′::typeof(k)

Return the reciprocal point `k′` with coordinates in a primitive basis, corresponding to
the input point `k` with coordinates in a conventional basis. 
"""
primitivize(::ReciprocalPoint, ::Union{Char, <:Integer})

# ---------------------------------------------------------------------------------------- #

function conventionalize(V′::Union{AbstractBasis{D}, AbstractPoint{D}}, cntr::Char) where D
    if cntr == 'P' || cntr == 'p' # the conventional and primitive bases coincide
        return V′
    else
        P = primitivebasismatrix(cntr, Val(D))
        return transform(V′, inv(P))
    end
end

function conventionalize(V′::Union{AbstractBasis{D}, AbstractPoint{D}}, sgnum::Integer) where D
    cntr = centering(sgnum, D)
    return conventionalize(V′, cntr)
end

@doc """ 
    conventionalize(V′::Union{AbstractBasis, AbstractPoint}, 
                    cntr_or_sgnum::Union{Char, <:Integer})    -->  V::typeof(V′)

Return the conventional basis or point `V` associated with the input primitive
`AbstractBasis` or `AbstractPoint` `V′`.

The assumed centering type is specified by `cntr_or_sgnum`, given either as a centering
character (`::Char`) or inferred from a space group number (`::Integer`) and the
dimensionality of `V′` (see also [`centering(::Integer, ::Integer)`](@ref)).
"""
conventionalize(::Union{AbstractBasis, AbstractPoint}, ::Union{Char, <:Integer})

@doc """
    conventionalize(Rs′::DirectBasis, cntr_or_sgnum::Union{Char, <:Integer}) --> Rs::typeof(Rs′)

Return the conventional direct basis `Rs` corresponding to the input primitive direct basis
`Rs′`.
"""
conventionalize(::DirectBasis, ::Union{Char, <:Integer})

@doc """
    conventionalize(Gs′::ReciprocalBasis, cntr_or_sgnum::Union{Char, <:Integer}) --> Gs::typeof(Gs′)

Return the conventional reciprocal basis `Gs` corresponding to the input primitive
reciprocal basis `Gs′`.
"""
conventionalize(::ReciprocalBasis, ::Union{Char, <:Integer})

@doc """
    conventionalize(r′::DirectPoint, cntr_or_sgnum::Union{Char, <:Integer}) --> r::typeof(r′)

Return the direct point `r` with coordinates in a conventional basis, corresponding to the
input point `r′` with coordinates in a primitive basis.
"""
conventionalize(::DirectPoint, ::Union{Char, <:Integer})

@doc """
    conventionalize(k′::ReciprocalPoint, cntr_or_sgnum::Union{Char, <:Integer}) --> k::typeof(k′)

Return the reciprocal point `k` with coordinates in a conventional basis, corresponding to
the input point `k′` with coordinates in a primitive basis. 
"""
conventionalize(::ReciprocalPoint, ::Union{Char, <:Integer})

# ---------------------------------------------------------------------------------------- #

const _basis_explain_str = "Depending on the object, the basis may be inferrable " * 
            "directly from the object; if not, it must be supplied explicitly."

@doc """
    cartesianize!

In-place transform an object with coordinates in an lattice basis to an object with
coordinates in a Cartesian basis.

$_basis_explain_str
"""
function cartesianize! end

@doc """
    cartesianize

Transform an object with coordinates in an lattice basis to an object with coordinates in a
Cartesian basis.

$_basis_explain_str
@doc """
function cartesianize end

@doc """
    latticize!

In-place transform object with coordinates in a Cartesian basis to an object with
coordinates in a lattice basis.

$_basis_explain_str
"""
function latticize! end

@doc """
    latticize

Transform an object with coordinates in a Cartesian basis to an object with coordinates in
a lattice basis.

$_basis_explain_str
"""
function latticize end

@doc """
    cartesianize(v::AbstractVector{<:Real}, basis)

Transform a vector `v` with coordinates referred to a lattice basis to a vector with
coordinates referred to the Cartesian basis implied by the columns (or vectors) of `basis`.
"""
cartesianize(v::AbstractVector{<:Real}, basis::AbstractMatrix{<:Real}) = basis*v
function cartesianize(v::AbstractVector{<:Real},
                      basis::AbstractVector{<:AbstractVector{<:Real}})
    return v'basis
end

@doc """
    latticize(v::AbstractVector{<:Real}, basis)

Transform a vector `v` with coordinates referred to the Cartesian basis to a vector with
coordinates referred to the lattice basis implied by the columns (or vectors) of `basis`.
"""
latticize(v::AbstractVector{<:Real}, basis::AbstractMatrix{<:Real}) = basis\v
function latticize(v::AbstractVector{<:Real},
                   basis::AbstractVector{<:AbstractVector{<:Real}})
    return latticize(v, reduce(hcat, basis))
end
