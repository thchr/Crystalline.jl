struct ModdedLinearRelation{T}
    irlabs    :: Vector{String}
    coefs     :: Vector{T}
    ambiguity :: T
end

struct LinearRelation{T}
    irlabs :: Vector{String}
    coefs  :: Vector{T}
    equal  :: T
end
function Base.:+(x::LinearRelation, y::LinearRelation)
    x.irlabs == y.irlabs || error("irrep labels must match")
    return LinearRelation(x.irlabs, x.coefs + y.coefs, x.equal + y.equal)
end
Base.:-(x::LinearRelation) = LinearRelation(x.irlabs, -x.coefs, -x.equal)
Base.:-(x::LinearRelation, y::LinearRelation) = x + (-y)

function Base.show(io::IO, f::Union{ModdedLinearRelation, LinearRelation}) 
    # trivial case
    if all(iszero, f.coefs) && isone(f.ambiguity)
        print(io, "(mod 1)")
        return
    end

    # prefactor
    prefactor = gcd(f.coefs)
    printed_prefactor = false
    if prefactor ≠ 1
        if denominator(prefactor) == 1
            print(io, numerator(prefactor))
        else
            print(io, "(", numerator(prefactor), "/", denominator(prefactor))
        end
        print(io, ")×(")
        printed_prefactor = true
    end

    # components in formula
    first = true
    for (i, c) in enumerate(f.coefs)
        iszero(c) && continue
        if !first
            print(io, signbit(c) ? " - " : " + ")
        else
            print(io, signbit(c) ? "-" : "")
            first = false
        end
        c′ = c//prefactor
        n, d = numerator(c′), denominator(c′)
        if isone(d) 
            if !isone(abs(n))
                print(io, abs(n))
            end
        else
            print(io, "(", abs(n), "/", d, ")")
        end
        print(io, f.irlabs[i])
    end
    printed_prefactor && print(io, ")")

    # mod part
    if f isa ModdedLinearRelation
        print(io, " (mod ", numerator(f.ambiguity))
        if !isone(denominator(f.ambiguity))
            print(io, "/", denominator(f.ambiguity))
        end
        print(io, ")")
    elseif f isa LinearRelation
        print(io, " = ", f.equal)
    end
end


"""
    corner_anomaly(brs::BandRepSet, Qd::Dict{String, Rational{Int}})

Compute a formula for the corner anomaly ```Q``` in terms of the irrep multiplicities of a
 band representation set `brs`. 
The fractional corner anomaly assigned to a minimal set of bands occupying each individual
Wyckoff position, typically obtained from an explicit counting argument, must be provided
as a dictionary `Qd`.

## Extended description
The corner anomaly gives the fractional integrated probability density (or, for photonic
systems, normalized modal energy density) per symmetry-related sector of the unit cell
(colloquially, per "corner"); for electronic systems, this quantity is a related to charge
and is often called the ``corner charge''.

The construction of the corner anomaly formula using the algorithm described by
> K. Naito, R. Takahashi, H. Watanabe, and S. Murakami, *Fractional hinge and corner charges 
> in various crystal shapes with cubic symmetry*, 
> [Phys. Rev. B **105**, 045126 (2022)](https://doi.org/10.1103/PhysRevB.105.045126).

## Example
We may compute the corner anomaly formula in plane group 2:
```jl
julia> Qd = Dict("1a"=>0//1, "1b"=>0//1, "1c"=>0//1, "1d"=>1//2)
julia> corner_anomaly(bandreps(2, 2), Qd)
Q = (1/4)×(3Y₁ + 3B₁ + A₁ + Γ₁) (mod 1)
```
"""
function corner_anomaly(brs::BandRepSet, Qd::Dict{String, Rational{Int}})
    relation_d = wyckoff_occupation(brs)
    coefs = sum(Qd[w]*relation.coefs for (w, relation) in relation_d)
    ambiguity_coefs = [Qd[w]*relation.ambiguity for (w, relation) in relation_d]

    ambiguity = gcd(ambiguity_coefs)
    coefs .= evenmod.(coefs, ambiguity)

    return ModdedLinearRelation(first(values(relation_d)).irlabs, coefs, ambiguity)
end

function wyckoff_occupation(brs::BandRepSet)
    A = matrix(brs; includedim=true)
    Nⁱʳ, Nᴱᴮᴿ = size(A, 1), size(A, 2)
    irlabs = length(brs.irlabs) == Nⁱʳ ? brs.irlabs : vcat(brs.irlabs, "μ")

    F = smith(A) # A = SΛT
    S⁻¹, T⁻¹ = F.Sinv, F.Tinv
    Λᵍ = diagm(Nᴱᴮᴿ, Nⁱʳ, map(λ -> iszero(λ) ? 0//1 : 1//λ, F.SNF))
    Aᵍ = T⁻¹*Λᵍ*S⁻¹ # generalized inverse of A

    # map global EBR-indices (i) to local Wyckoff position indices (w,α), with w as Dict
    # keys and α = 1,…,n as local, linear indices per Wyckoff position
    wα_idxs = Dict{String, Vector{Int}}()
    w_mults = Dict{String, Int}()
    for (i, br) in enumerate(brs)
        w = br.wyckpos
        if haskey(wα_idxs, w)
            push!(wα_idxs[w], i)
        else
            wα_idxs[w] = [i]
            w_mults[w] = parse(Int, first(w))
        end
    end
    wps = keys(wα_idxs)

    siteir_dims_d = Dict{String, Dict{Int, Int}}()  # dim(w|α)/multiplicity(w) = dim(ρʷ_α)
    for (w, αs) in wα_idxs
        siteir_dims_d[w] = Dict(i => convert(Int, brs[i][end]//w_mults[w]) for i in αs)
    end
    
    # coefficients in the "occupation" formula
    δ(i,j) = i==j ? 1 : 0
    coefs_d           = Dict(wp => zeros(Rational{Int}, Nⁱʳ) for wp in wps) # Q = coefs⋅n mod p
    ambiguity_coefs_d = Dict(wp => zeros(Int, Nᴱᴮᴿ) for wp in wps)          # p = (..)⋅y
    
    AᵍA = convert.(Int, Aᵍ*A) # integer elements by definition
    for (w, αs) in wα_idxs
        for j in 1:Nⁱʳ
            coefs_d[w][j] = sum(siteir_dims_d[w][i]*Aᵍ[i,j] for i in αs)
        end
        for k in 1:Nᴱᴮᴿ
            ambiguity_coefs_d[w][k] = sum(siteir_dims_d[w][i]*(-δ(i,k) + AᵍA[i,k]) for i in αs)
        end
    end

    relation_d = Dict{String, ModdedLinearRelation}()
    for w in wps
        ambiguity = gcd(ambiguity_coefs_d[w])
        coefs′ = evenmod.(coefs_d[w], ambiguity) # reduce coefficient mod ambiguity
        relation_d[w] = ModdedLinearRelation(irlabs, coefs′, Rational(ambiguity))
    end
    return relation_d
end

# analogous to `mod(x, r)`, but returns an answer in an even interval around 0, i.e. in 
# `]-r/2, r/2]` rather than in `[0, r[`
function evenmod(x, r)
    y = mod(x, r)
    return y > r//2 ? y-r : y
end

function compatibility_relations(F::Smith)
    S⁻¹ = F.Sinv
    dᵇˢ = count(!iszero, F.SNF)
    return S⁻¹[dᵇˢ+1:end,:]
end

function pretty_compatibility_relations(brs)
    F = smith(matrix(brs; includedim=true))
    C = compatibility_relations(F)
    irlabs = vcat(brs.irlabs, "μ")
    return [LinearRelation(irlabs, Vector(r), 0) for r in eachrow(C)]
end


function compatibility_basis(F::Smith)
     # a basis of vectors that obey the compatibility relations; the non-empty null-space
     # of the compability relation matrix C
     C = compatibility_relations(F)
     Cᵍ = generalized_inv(C)
     nullC = convert.(Int, I-Cᵍ*C)
     return filter(!iszero, eachcol(nullC))
end

function generalized_inv(X)
    F = smith(X)
    Λᵍ = diagm(size(X,2), size(X,1), [iszero(x) ? 0//1 : 1//x for x in F.SNF])
    return F.Tinv*Λᵍ*F.Sinv
end