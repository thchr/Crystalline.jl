## --- TYPES ---

abstract type AbstractFourierLattice{D}; end
getcoefs(flat::AbstractFourierLattice) = flat.orbitcoefs
getorbits(flat::AbstractFourierLattice) = flat.orbits
dim(::AbstractFourierLattice{D}) where D = D
function (==)(flat1::AbstractFourierLattice, flat2::AbstractFourierLattice)
    return flat1.orbits == flat2.orbits && flat1.orbitcoefs == flat2.orbitcoefs
end
function Base.isapprox(flat1::AbstractFourierLattice, flat2::AbstractFourierLattice; kwargs...)
    return ( isapprox(flat1.orbits,     flat2.orbits;     kwargs...) && 
             isapprox(flat1.orbitcoefs, flat2.orbitcoefs; kwargs...) )
end

@doc """
    UnityFourierLattice{D} <: AbstractFourierLattice{D}

A general `D`-dimensional Fourier (plane wave) lattice specified by orbits of
reciprocal lattice vectors (`orbits`) and coefficient interrelations (`orbitcoefs`)).
The norm of all elements in `orbitcoefs` is unity. `orbits` (and associated
coefficients) are sorted in order of increasing norm (low to high).
"""
struct UnityFourierLattice{D} <: AbstractFourierLattice{D}
    orbits::Vector{Vector{SVector{D, Int}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}  # Vector of interrelations between coefficients of 𝐆-plane waves within an orbit; unit norm
end

@doc """
    ModulatedFourierLattice{D} <: AbstractFourierLattice{D}

A `D`-dimensional concrete Fourier (plane wave) lattice, derived from 
a [`UnityFourierLattice{D}`](@ref) by scaling and modulating its orbit coefficients 
by complex numbers; in general, the coefficients do not have unit norm.
"""
struct ModulatedFourierLattice{D} <: AbstractFourierLattice{D}
    orbits::Vector{Vector{SVector{D, Int}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}  # Vector of coefficients of 𝐆-plane waves within an orbit
end


## --- METHODS --- 

# Group orbits of plane waves G = (G)ᵀ under a symmetry operation Ô = {W|w}, 
# using that Ô acts as Ô⁻¹={W⁻¹|-W⁻¹w} when acting on functions, i.e.
#   Ôexp(iG⋅r) = Ôexp(iG⋅Ô⁻¹r) = exp[iG⋅(W⁻¹r-W⁻¹w)]
# and 
#   exp(iG⋅W⁻¹r) = exp(iGᵀW⁻¹r) = exp{i[(W⁻¹)ᵀG]ᵀ⋅r}
@doc """
    levelsetlattice(sgnum::Integer, D::Integer=2, idxmax::NTuple=ntuple(i->2,D))
        --> UnityFourierLattice{D}

Compute a "neutral"/uninitialized Fourier lattice basis, a [`UnityFourierLattice`](@ref),
consistent with the symmetries of the space group `sgnum` in dimension `D`. 
The resulting lattice `flat` is expanded in a Fourier basis split into symmetry-derived
orbits, with intra-orbit coefficients constrained by the symmetries of the space-group.
The inter-orbit coefficients are, however, free and unconstrained.

The Fourier resolution along each reciprocal lattice vector is controlled by `idxmax`:
e.g., if `D = 2` and `idxmax = (2, 3)`, the resulting Fourier lattice may contain 
reciprocal lattice vectors (k₁, k₂) with k₁∈[0,±1,±2] and k₂∈[0,±1,±2,±3], referred 
to a 𝐆-basis.

This "neutral" lattice can, and usually should, be subsequently modulated by
[`modulate`](@ref) (which modulates the inter-orbit coefficients, which may eliminate
"synthetic symmetries" that can exist in the "neutral" configuration, due to all 
inter-orbit coefficients being set to unity).

## Examples

Compute a `UnityFourierLattice`, modulate it with random inter-orbit coefficients 
via `modulate`, and finally plot it (via PyPlot.jl):

```julia-repl
julia> uflat = levelsetlattice(16, Val(2))
julia> flat  = modulate(uflat)
julia> Rs    = directbasis(16, Val(2)) 
julia> using PyPlot
julia> plot(flat, Rs)
```
"""
function levelsetlattice(sgnum::Integer, D::Integer=2,
                         idxmax::NTuple{D′,Int}=ntuple(i->2, D)) where D′
    D == D′ || throw(DomainError((D=D, idxmax=idxmax), "incompatible dimensions"))
    return levelsetlattice(sgnum, Val(D′), idxmax)
end
function levelsetlattice(sgnum::Integer, Dᵛ::Val{D}, idxmax::NTuple{D,Int}=ntuple(i->2, D)) where D
    # check validity of inputs
    (sgnum < 1)             && throw(DomainError(sgnum, "sgnum must be greater than 1"))
    D ∉ (1,2,3)             && _throw_invalid_dim(D)
    D ≠ length(idxmax)      && throw(DomainError((D, idxmax), "D must equal length(idxmax): got (D = $D) ≠ (length(idxmax) = $(length(idxmax)))"))
    (D == 2 && sgnum > 17)  || (D == 3 && sgnum > 230) && throw(DomainError(sgnum, "sgnum must be in range 1:17 in 2D and in 1:230 in 3D"))

    # prepare
    sg = spacegroup(sgnum, Dᵛ)
    sgops = operations(sg)
    Ws = rotation.(sgops) # operations W in R-basis (point group part)
    ws = translation.(sgops)

    # We define the "reciprocal orbit" associated with the action of W through (W⁻¹)ᵀ
    # calculating the operators (W⁻¹)ᵀ in the 𝐆-basis:
    # The action of a symmetry operator in an 𝐑-basis, i.e. W(𝐑), on a 𝐤 vector in a 
    # 𝐆-basis, i.e. 𝐤(𝐆), is 𝐤′(𝐆)ᵀ = 𝐤(𝐆)ᵀW(𝐑)⁻¹. To deal with column vectors, we 
    # transpose, obtaining 𝐤′(𝐆) = [W(𝐑)⁻¹]ᵀ𝐤(𝐆) [details in symops.jl, above littlegroup(...)].
    W⁻¹ᵀs = transpose.(inv.(Ws))

    # If idxmax is interpreted as (imax, jmax, ...), then this produces an iterator
    # over i = -imax:imax, j = -jmax:jmax, ..., where each call returns (..., j, i); 
    # note that the final order is anti-lexicographical; so we reverse it in the actual
    # loop for our own sanity's sake
    reviter = Iterators.product(reverse((:).(.-idxmax, idxmax))...)

    # --- compute orbits ---
    orbits = Vector{Vector{SVector{D,Int}}}() # vector to store orbits of G-vectors (in G-basis)
    for rG in reviter  
        G = SVector{D,Int}(reverse(rG)) # fix order and convert to SVector{D,Int} from Tuple

        skip = false # if G already contained in an orbit; go to next G
        for orb in orbits
            isapproxin(G, orb) && (skip=true; break) 
        end
        skip && continue
        
        neworb = _orbit(W⁻¹ᵀs, G) # compute orbit assoc with G-vector
        # the symmetry transformation may introduce round-off errors, but we know that 
        # the indices must be integers; fix that here, and check its validity as well
        neworb′ = [round.(Int,G′) for G′ in neworb] 
        if norm(neworb′ .- neworb) > DEFAULT_ATOL; 
            error("The G-combinations and their symmetry-transforms must be integers"); 
        end
        push!(orbits, neworb′) # add orbit to list of orbits
    end

    # --- restrictions on orbit coeffs. due to nonsymmorphic elements in space group ---
    orbitcoefs = Vector{Vector{ComplexF64}}()
    deleteidx = Vector{Int}()
    for (o,orb) in enumerate(orbits)
        start = true; prevspan = []
        for (W⁻¹ᵀ, w) in zip(W⁻¹ᵀs, ws)
            conds = zeros(ComplexF64, length(orb), length(orb))
            for (m, G) in enumerate(orb)
                G′ = W⁻¹ᵀ*G  # planewave G is transformed to by W⁻¹ᵀ
                diffs = norm.(Ref(G′) .- orb); 
                n = argmin(diffs) # find assoc linear index in orbit
                diffs[n] > DEFAULT_ATOL && error("Part of an orbit was miscalculated; diff = $(diffs[n])")
                # the inverse translation is -W⁻¹w; the phase is thus exp(-iG⋅W⁻¹w) which
                # is equivalent to exp[-i(W⁻¹ᵀG)w]. We use the latter, so we avoid an
                # unnecessary matrix-vector product [i.e. dot(G, W⁻¹w) = dot(G′, w)]
                conds[n,m] = cis(-2π*dot(G′, w)) # cis(x) = exp(ix)
            end

            nextspan = nullspace(conds-I, atol=NULL_ATOL)          
            if start
                prevspan = nextspan
                start = false
            elseif !isempty(prevspan) && !isempty(nextspan)
                spansect = nullspace([prevspan -nextspan], atol=NULL_ATOL)[size(prevspan, 2)+1:end,:]
                prevspan = nextspan*spansect
            else
                prevspan = nothing; break
            end
        end
                    
        if !isnothing(prevspan)
            if size(prevspan,2) != 1; error("Unexpected size of prevspan"); end
            coefbasis = vec(prevspan)
            coefbasis ./= coefbasis[argmax(norm(coefbasis, Inf))]
            push!(orbitcoefs, coefbasis)
        else 
            push!(deleteidx, o)
        end
    end
    deleteat!(orbits, deleteidx)

    # sort in order of descending wavelength (e.g., [0,0,...] term comes first; highest G-combinations come last)
    perm = sortperm(orbits, by=x->norm(first(x)))
    permute!(orbits, perm)
    permute!(orbitcoefs, perm)

    return UnityFourierLattice{D}(orbits, orbitcoefs)
end


@doc """
    _orbit(Ws, x)

Computes the orbit of a direct-space point `x` under a set of point-group operations `Ws`,
i.e. computes the set {gx | g∈G} where g denotes elements of the group G composed of all
operations in `Ws` (possibly iterated, to ensure full coverage).

It is important that `Ws` and `x` are given in the same basis. 

[``W' = PWP⁻¹`` if the basis change is from coordinates r to r' = Pr, corresponding 
to a new set of basis vectors (x̂')ᵀ=x̂ᵀP; e.g., when going from a direct basis
representation to a Cartesian one, the basis change matrix is P = [R₁ R₂ R₃],
with Rᵢ inserted as column vectors]
"""
function _orbit(Ws::AbstractVector{<:AbstractMatrix{<:Real}}, x::AbstractVector{<:Real})
    fx = float.(x)
    xorbit = [fx]
    for W in Ws
        x′ = fx
        while true
            x′ = W*x′
            if !isapproxin(x′, xorbit)
                push!(xorbit, x′)
            else 
                break
            end
        end
    end
    return sort!(xorbit) # convenient to sort it before returning, for future comparisons
end

function transform(flat::AbstractFourierLattice{D}, P::AbstractMatrix{<:Real}) where D
    # The orbits consist of G-vector specified as a coordinate vector 𝐤≡(k₁,k₂,k₃)ᵀ, referred
    # to the untransformed 𝐆-basis (𝐚* 𝐛* 𝐜*), and we want to instead express it as a coordinate
    # vector 𝐤′≡(k₁′,k₂′,k₃′)ᵀ referred to the transformed 𝐆-basis (𝐚*′ 𝐛*′ 𝐜*′)≡(𝐚* 𝐛* 𝐜*)(P⁻¹)ᵀ,
    # where P is the transformation matrix. This is achieved by transforming according to 𝐤′ = Pᵀ𝐤
    # or, equivalently, (k₁′ k₂′ k₃′)ᵀ = Pᵀ(k₁ k₂ k₃)ᵀ. See also `transform(::KVec, ...)` and 
    # `transform(::ReciprocalBasis, ...)`.

    # vec of vec of G-vectors (in a **untransformed** 𝐆-basis)
    orbits = getorbits(flat)
    # prealloc. a vec of vec of k-vecs (to be filled in the **transformed** 𝐆-basis)
    orbits′ = [Vector{SVector{D, Int}}(undef, length(orb)) for orb in orbits]
    # transform all k-vecs in the orbits
    for (i, orb) in enumerate(orbits)
        for (j, k) in enumerate(orb)
            k′ = P'*k
            int_k′ = round.(Int, k′)
            if !isapprox(k′, int_k′, atol=DEFAULT_ATOL)
                error("unexpectedly obtained non-integer k-vector in orbit")
            end
            orbits′[i][j] = int_k′ :: SVector{D,Int}
        end
    end
    # --- Comment regarding the `convert(SVector{D, Int}, ...)` call above: ---
    # Because primitive reciprocal basis Gs′≡(𝐚*′ 𝐛*′ 𝐜*′) consists of "larger" vectors
    # than the conventional basis Gs≡(𝐚* 𝐛* 𝐜*) (since the direct lattice shrinks when we
    # go to a primitive basis), not every conventional reciprocal lattice coordinate vector
    # 𝐤 has a primitive integer-coordinate vector 𝐤′=Pᵀ𝐤 (i.e. kᵢ∈ℕ does not imply kᵢ′∈ℕ).
    # However, since `flat` is derived consistent with the symmetries in a conventional
    # basis, the necessary restrictions will already have been imposed in the creation of
    # `flat` so that the primitivized version will have only integer coefficients (otherwise
    # the lattice would not be periodic in the primitive cell). I.e. we need not worry that
    # the conversion is impossible, so long that we transform to a meaningful basis.
    # The same issue of course isn't relevant for transforming in the reverse direction.
       
    # the coefficients of flat are unchanged; only the 𝐑- and 𝐆-basis change
    return typeof(flat)(orbits′, deepcopy(getcoefs(flat))) # return in the same type as `flat`
end

@doc raw"""
    primitivize(flat::AbstractFourierLattice, cntr::Char) --> ::typeof(flat)

Given `flat` referred to a conventional basis with centering `cntr`, compute the derived
(but physically equivalent) lattice `flat′` referred to the associated primitive basis. 

Specifically, if `flat` refers to a direct conventional basis `Rs`
``≡ (\mathbf{a} \mathbf{b} \mathbf{c})`` [with coordinate vectors
``\mathbf{r} ≡ (r_1, r_2, r_3)^{\mathrm{T}}``] then `flat′` refers to a direct primitive
basis `Rs′`
``≡ (\mathbf{a}' \mathbf{b}' \mathbf{c}') ≡ (\mathbf{a} \mathbf{b} \mathbf{c})\mathbf{P}``
[with coordinate vectors 
``\mathbf{r}' ≡ (r_1', r_2', r_3')^{\mathrm{T}} = \mathbf{P}^{-1}\mathbf{r}``], where
``\mathbf{P}`` denotes the basis-change matrix obtained from `primitivebasismatrix(...)`.

To compute the associated primitive basis vectors, see
[`primitivize(::DirectBasis, ::Char)`](@ref) [specifically, `Rs′ = primitivize(Rs, cntr)`].

## Examples

A centered ('c') lattice from plane group 5 in 2D, plotted in its 
conventional and primitive basis (requires `using PyPlot`):

```julia-repl
julia> sgnum = 5; D = 2; cntr = centering(sgnum, D)  # 'c' (body-centered)

julia> Rs   = directbasis(sgnum, Val(D))     # conventional basis (rectangular)
julia> flat = levelsetlattice(sgnum, Val(D)) # Fourier lattice in basis of Rs
julia> flat = modulate(flat)                 # modulate the lattice coefficients
julia> plot(flat, Rs)

julia> Rs′   = primitivize(Rs, cntr)    # primitive basis (oblique)
julia> flat′ = primitivize(flat, cntr)  # Fourier lattice in basis of Rs′

julia> using PyPlot
julia> plot(flat′, Rs′)
```
"""
function primitivize(flat::AbstractFourierLattice{D}, cntr::Char) where D
    # Short-circuit for lattices that have trivial transformation matrices
    (D == 3 && cntr == 'P') && return flat
    (D == 2 && cntr == 'p') && return flat
    D == 1 && return flat

    P = primitivebasismatrix(cntr, Val(D))
    return transform(flat, P)
end

"""
    conventionalize(flat::AbstractFourierLattice, cntr::Char) --> ::typeof(flat′)

Given `flat` referred to a primitive basis with centering `cntr`, compute the derived (but
physically equivalent) lattice `flat′` referred to the associated conventional basis. 

See also the complementary method
[`primitivize(::AbstractFourierLattice, ::Char)`](@ref) for additional details.
"""
function conventionalize(flat::AbstractFourierLattice{D}, cntr::Char) where D
    # Short-circuit for lattices that have trivial transformation matrices
    (D == 3 && cntr == 'P') && return flat
    (D == 2 && cntr == 'p') && return flat
    D == 1 && return flat

    P = primitivebasismatrix(cntr, Val(D))
    return transform(flat, inv(P))
end

@doc """
    modulate(flat::UnityFourierLattice{D},
    modulation::AbstractVector{ComplexF64}=rand(ComplexF64, length(getcoefs(flat))),
    expon::Union{Nothing, Real}=nothing, Gs::Union{ReciprocalBasis{D}, Nothing}=nothing)
                            --> ModulatedFourierLattice{D}

Derive a concrete, modulated Fourier lattice from a `UnityFourierLattice` `flat`
(containing the _interrelations_ between orbit coefficients), by 
multiplying the "normalized" orbit coefficients by a `modulation`, a _complex_
modulating vector (in general, should be complex; otherwise restores unintended
symmetry to the lattice). Distinct `modulation` vectors produce distinct 
realizations of the same lattice described by the original `flat`. By default,
a random complex vector is used.

An exponent `expon` can be provided, which introduces a penalty term to short-
wavelength features (i.e. high-|G| orbits) by dividing the orbit coefficients
by |G|^`expon`; producing a "simpler" and "smoother" lattice boundary
when `expon > 0` (reverse for `expon < 0`). This basically amounts to a 
continuous "simplifying" operation on the lattice (it is not necessarily a 
smoothing operation; it simply suppresses "high-frequency" components).
If `expon = nothing`, no rescaling is performed. If `Gs` is provided as `nothing`,
the orbit norm is computed in the reciprocal lattice basis (and, so, may not strictly
speaking be a norm if the lattice basis is not cartesian); to account for the basis
explicitly, `Gs` must be provided as a [`ReciprocalBasis`](@ref), see also
[`normscale`](@ref).
"""
function modulate(flat::AbstractFourierLattice{D},
                  modulation::Union{Nothing, AbstractVector{ComplexF64}}=nothing,
                  expon::Union{Nothing, Real}=nothing,
                  Gs::Union{ReciprocalBasis{D}, Nothing}=nothing) where D
    if isnothing(modulation)
        Ncoefs = length(getcoefs(flat))
        mod_r, mod_ϕ = rand(Float64, Ncoefs), 2π.*rand(Float64, Ncoefs)
        modulation = mod_r .* cis.(mod_ϕ) # ≡ reⁱᵠ (pick modulus and phase uniformly random)
    end
    orbits = getorbits(flat); orbitcoefs = getcoefs(flat) # unpacking ...
    
    # multiply the orbit coefficients by the overall `modulation` vector
    modulated_orbitcoefs = orbitcoefs.*modulation
    flat′ = ModulatedFourierLattice{D}(orbits, modulated_orbitcoefs)

    if !isnothing(expon) && !iszero(expon)
        # `expon ≠ 0`: interpret as a penalty term on short-wavelength orbits (high |𝐆|)
        # by dividing the orbit coefficients by |𝐆|^`expon`; producing smoother lattice
        # boundaries and simpler features for `expon > 0` (reverse for `expon < 0`)
        normscale!(flat′, expon, Gs)
    end

    return flat′
end

@doc """ 
    normscale(flat::ModulatedFourierLattice, expon::Real, 
              Gs::Union{ReciprocalBasis, Nothing} = nothing)  --> ModulatedFourierLattice

Applies inverse-orbit norm rescaling of expansion coefficients with a norm exponent `expon`.
If `Gs` is nothing, the orbit norm is computed in the lattice basis (and, so, is not
strictly a norm); by providing `Gs` as [`ReciprocalBasis`](@ref), the norm is evaluated
correctly in cartesian setting. See further discussion in [`modulate`](@ref).

An in-place equivalent is provided in [`normscale!`](@ref).
"""
function normscale(flat::ModulatedFourierLattice{D}, expon::Real,
                   Gs::Union{ReciprocalBasis{D}, Nothing} = nothing)  where D
    return normscale!(deepcopy(flat), expon, Gs)
end
@doc """
    normscale!(flat::ModulatedFourierLattice, expon::Real,
               Gs::Union{ReciprocalBasis, Nothing} = nothing) --> ModulatedFourierLattice

In-place equivalent of [`normscale`](@ref): mutates `flat`.
"""
function normscale!(flat::ModulatedFourierLattice{D}, expon::Real,
                    Gs::Union{ReciprocalBasis{D}, Nothing} = nothing) where D
    n₀ = isnothing(Gs) ? 1.0 : sum(norm, Gs) / D
    if !iszero(expon)
        orbits = getorbits(flat)
        @inbounds for i in eachindex(orbits)
            n = if isnothing(Gs)
                norm(first(orbits[i]))
            else
                norm(first(orbits[i])'*Gs) / n₀
            end
            rescale_factor = n^expon
            rescale_factor == zero(rescale_factor) && continue # for G = [0,0,0] case
            flat.orbitcoefs[i] ./= rescale_factor
        end
    end
    return flat
end

# -----------------------------------------------------------------------------------------
# The utilities and methods below are mostly used for plotting (see src/pyplotting.jl).
# We keep them here since they do not depend on PyPlot and have more general utility in 
# principle (e.g., exporting associated Meshes).


@doc raw"""
    (flat::AbstractFourierLattice)(xyz) --> Float64
    (flat::AbstractFourierLattice)(xyzs...) --> Float64

Evaluate an `AbstractFourierLattice` at the point `xyz` and return its real part, i.e.
    
```math
    \mathop{\mathrm{Re}}\sum_i c_i \exp(2\pi i\mathbf{G}_i\cdot\mathbf{r})
```

with ``\mathrm{G}_i`` denoting reciprocal lattice vectors in the allowed orbits of `flat`,
with ``c_i`` denoting the associated coefficients (and ``\mathbf{r} \equiv`` `xyz`).

`xyz` may be any iterable object with dimension matching `flat` consisting of real numbers
(e.g., a `Tuple`, `Vector`, or `SVector`). Alternatively, the coordinates can be supplied
individually (i.e., as `flat(x, y, z)`).
"""
function (flat::AbstractFourierLattice)(xyz)
    dim(flat) == length(xyz) || throw(DimensionMismatch("dimensions of flat and xyz are mismatched"))
    orbits = getorbits(flat)
    coefs  = getcoefs(flat)
    f = zero(Float64)
    for (orb, cs) in zip(orbits, coefs)
        for (G, c) in zip(orb, cs)
            # though one might naively think the phase would need a conversion between 
            # 𝐑- and 𝐆-bases, this is not necessary since P(𝐆)ᵀP(𝐑) = 2π𝐈 by definition
            exp_im, exp_re = sincos(2π*dot(G, xyz))
            f += real(c)*exp_re - imag(c)*exp_im    # ≡ real(c*exp(2π*1im*dot(G, xyz)))
        end
    end
    return f
end
(flat::AbstractFourierLattice{D})(xyzs::Vararg{Real, D}) where D = flat(SVector{D, Float64}(xyzs))

# note: 
# the "(x,y,z) ordering" depends on dimension D:
#       D = 2: x runs over cols (dim=2), y over rows (dim=1), i.e. "y-then-x"
#       D = 3: x runs over dim=1, y over dim=2, z over dim=3, i.e. "x-then-y-then-z"
# this is because plotting utilities usually "y-then-x", but e.g. Meshing.jl (for 3D
# isosurfaces) assumes the the more natural "x-then-y-then-z" sorting used here. This does
# require some care though, because if we export the output of a 3D calculation to Matlab,
# to use it for isosurface generation, it again requires a sorting like "y-then-x-then-z",
# so we need to permute dimensions 1 and 2 of the output of `calcfouriergridded` when used 
# with Matlab.
function calcfouriergridded!(vals, xyz, flat::AbstractFourierLattice{D}, 
                             N::Integer=length(xyz)) where D

    # evaluate `flat` over all gridpoints via broadcasting
    if D == 2
        # x along columns, y along rows: "y-then-x"
        broadcast!(flat, vals, reshape(xyz, (1,N)), reshape(xyz, (N,1)))
    elseif D == 3
        # x along dim 1, y along dim 2, z along dim 3: "x-then-y-then-z", equivalent to a
        # triple loop, ala `for x∈xyz, y∈xyz, z∈xyz; vals[ix,iy,iz] = f(x,y,z); end`
        broadcast!(flat, vals, reshape(xyz, (N,1,1)), reshape(xyz, (1,N,1)), reshape(xyz, (1,1,N)))
    end
    return vals
end
function calcfouriergridded(xyz, flat::AbstractFourierLattice{D},
                            N::Integer=length(xyz)) where D
    vals = Array{Float64, D}(undef, ntuple(i->N, D)...)
    return calcfouriergridded!(vals, xyz, flat, N)
end

# ---------------------------------------------------------------------------------------- #

# returns a lazy "vector" of the (real) value of `flat` over the entire unit cell, using
# `nsamples` per dimension; avoids double counting at unit cell edges
function _lazy_fourier_eval_over_uc_as_vector(
    flat::AbstractFourierLattice{D}, 
    nsamples::Integer
) where D
    step = 1.0/nsamples
    # last sample point is `0.5-step` to avoid double counting ±0.5 (cf. periodicity)
    samples = range(-0.5, 0.5-step, length=nsamples)
    if D == 1
        itr = (real(flat(x)) for x in samples)
    elseif D == 2
        itr = (real(flat(x,y)) for x in samples for y in samples)
    elseif D == 3
        itr = (real(flat(x,y,z)) for x in samples for y in samples for z in samples)
    end
end

"""
    $(TYPEDSIGNATURES)

Return the isovalue of `flat` such that the interior of the thusly defined levelset
isosurface encloses a fraction `filling` of the unit cell.

The keyword argument `nsamples` specifies the grid-resolution used in evaluating this
answer (via Statistics.jl's `quantile`).
"""
function filling2isoval(
    flat::AbstractFourierLattice{D}, 
    filling::Real=0.5,
    nsamples::Integer=51
) where D
    itr = _lazy_fourier_eval_over_uc_as_vector(flat, nsamples)
    return quantile(itr, filling)
end

"""
    $(TYPEDSIGNATURES)

Return the filling fraction of the interior of the isosurface defined by `flat` for an
isovalue `isoval`.

The keyword argument `nsamples` specifies the grid-resolution used in evaluating this
answer (using staircase integration).
"""
function isoval2filling(
    flat::AbstractFourierLattice{D}, 
    isoval::Real, 
    nsamples::Integer=51
) where D
    itr = _lazy_fourier_eval_over_uc_as_vector(flat, nsamples)
    return count(<(isoval), itr)/nsamples^D
end