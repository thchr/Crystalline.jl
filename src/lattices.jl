## --- TYPES ---

abstract type AbstractFourierLattice{D}; end
getcoefs(flat::AbstractFourierLattice) = flat.orbitcoefs
getorbits(flat::AbstractFourierLattice) = flat.orbits
dim(flat::AbstractFourierLattice{D}) where D = D

"""
UnityFourierLatticeFourierLattice{D} <: AbstractFourierLattice{D}

A general `D`-dimensional Fourier/plane wave lattice (specified 
by G-orbits and coefficient interrelations); specifies the allowable 
interrelations between coefficients within each orbit. The norm of 
all orbit coefficients is unity. The G-orbits `orbits` (& associated
coefficients) are sorted in order of increasing |G| (low to high).
"""
struct UnityFourierLattice{D} <: AbstractFourierLattice{D}
    orbits::Vector{Vector{SVector{D, Int}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}  # Vector of interrelations between coefficients of 𝐆-plane waves within an orbit; unit norm
end
function UnityFourierLattice(orbits, orbitcoefs)
    D = length(first(first(orbits)))
    UnityFourierLattice{D}(orbits, orbitcoefs)
end


"""
    ModulatedFourierLattice{D} <: AbstractFourierLattice{D}

A `D`-dimensional concrete Fourier/plane wave lattice, derived from 
a UnityFourierLattice by scaling/modulating its orbit coefficients 
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
"""
    levelsetlattice(sgnum::Integer, D::Integer=2, idxmax::NTuple=ntuple(i->2,dim))
        --> UnityFourierLattice{D}

Compute a "neutral"/uninitialized Fourier lattice basis, a UnityFourierLattice, consistent
with the symmetries of the space group `sgnum` in dimension `D`. The resulting lattice
`flat` is expanded in a Fourier basis split into symmetry-derived orbits, with intra-orbit 
coefficients constrained by the symmetries of the space-group. The inter-orbit coefficients
are, however, free and unconstrained.

The Fourier resolution along each reciprocal lattice vector is controlled by `idxmax`:
e.g., if `D = 2` and `idxmax = (2, 3)`, the resulting Fourier lattice may contain 
reciprocal lattice vectors (k₁, k₂) with k₁∈[0,±1,±2] and k₂∈[0,±1,±2,±3], referred 
to a 𝐆-basis.

This "neutral" lattice can, and usually should, be subsequently modulated by `modulate`
(modulates the inter-orbit coefficients, which will often eliminate symmetries that may
remain in the "neutral" configuration, where all inter-orbit coefficients are unity).

# Examples

Compute a UnityFourierLattice, modulate it with random inter-orbit coefficients via `modulate`,
and finally plot it:

```julia-repl
julia> uflat = levelsetlattice(16, 2)
julia> flat  = modulate(uflat)
julia> Rs    = directbasis(16, 2) 
julia> plot(flat, Rs)
```
"""
function levelsetlattice(sgnum::Integer, D::Integer=2, 
                         idxmax::NTuple=ntuple(i->2, D))
    # check validity of inputs
    (sgnum < 1)               && throw(DomainError(sgnum, "sgnum must be greater than 1"))
    !(D == 2 || D == 3)   && throw(DomainError(D, "D must be equal to 2 or 3"))
    D ≠ length(idxmax)      && throw(DomainError((D, idxmax), "D must equal length(idxmax): got (D = $D) ≠ (length(idxmax) = $(length(idxmax)))"))
    (D == 2 && sgnum > 17)  || (D == 3 && sgnum > 230) && throw(DomainError(sgnum, "sgnum must be in range 1:17 in 2D and in 1:230 in 3D"))

    # prepare
    sg = get_sgops(sgnum, D)
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
    orbits = Vector{Vector{SVector{D,Int64}}}() # vector to store orbits of G-vectors (in G-basis)
    for rG in reviter  
        G = SVector{D,Int64}(reverse(rG)) # fix order and convert to SVector{D,Int64} from Tuple

        skip = false # if G already contained in an orbit; go to next G
        for orb in orbits
            isapproxin(G, orb) && (skip=true; break) 
        end
        skip && continue
        
        neworb = orbit(W⁻¹ᵀs, G) # compute orbit assoc with G-vector
        # the symmetry transformation may introduce round-off errors, but we know that 
        # the indices must be integers; fix that here, and check its validity as well
        neworb′ = [round.(Int64,G′) for G′ in neworb] 
        if norm(neworb′ .- neworb) > DEFAULT_ATOL; 
            error("The G-combinations and their symmetry-transforms must be integers"); 
        end
        push!(orbits, neworb′) # add orbit to list of orbits
    end

    # --- restrictions on orbit coeffs. due to nonsymmorphic elements in space group ---
    orbitcoefs = Vector{Vector{ComplexF64}}()
    deleteidx = Vector{Int64}()
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

    return UnityFourierLattice(orbits, orbitcoefs)
end


"""
    orbit(Ws, x)

Computes the orbit of `x` under a set of point-group operations `Ws`,
i.e. computes the set `{gx | g∈G}` where `g` denotes elements of the group
`G` composed of all operations in `Ws` (possibly iterated, to ensure
full coverage).

At the moment, we only consider _point group_ operations; i.e. there are 
no nonsymmorphic `Ws` parts. 

It is important that `Ws` and `x` are given in the same basis. 

[W' = PWP⁻¹ if the basis change is from coordinates r to r' = Pr, corresponding 
to a new set of basis vectors (x̂')ᵀ=x̂ᵀP; e.g., when going from a direct basis
representation to a Cartesian one, the basis change matrix is P = [R₁ R₂ R₃],
with Rᵢ inserted as column vectors]
"""
function orbit(Ws::AbstractVector{<:AbstractMatrix{<:Real}}, x::AbstractVector{<:Real})
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

"""
    primitivize(flat::AbstractFourierLattice, cntr::Char) --> AbstractFourierLattice

Given `flat` referred to a conventional basis with centering `cntr`, compute the 
derived (but physically equivalent) lattice `flat′` referred to the associated 
primitive basis. 

Specifically, if `flat` refers to a direct conventional basis Rs≡(𝐚 𝐛 𝐜) [with 
coordinate-vectors 𝐫≡(r₁, r₂, r₃)ᵀ] then `flat′` refers to a direct primitive 
basis Rs′≡(𝐚′ 𝐛′ 𝐜′)≡(𝐚 𝐛 𝐜)P [with coordinate-vectors 𝐫′≡(r₁′, r₂′, r₃′)ᵀ=P⁻¹𝐫],
where P denotes the basis-change matrix obtained from `primitivebasismatrix(...)`.

To compute the associated primitive basis vectors, see `primitivize(::DirectBasis)`
[specifically, `Rs′ = primitivize(Rs, cntr)`].


# Examples

A centered ('c') lattice from plane group 5 in 2D, plotted in its 
conventional and primitive basis:

```julia-repl
julia> sgnum = 5; D = 2; cntr = centering(sgnum, D)  # 'c' (body-centered)

julia> Rs   = directbasis(sgnum, D)     # conventional basis (rectangular)
julia> flat = levelsetlattice(sgnum, D) # Fourier lattice in basis of Rs
julia> flat = modulate(flat)            # modulate the lattice coefficients
julia> plot(flat, Rs)

julia> Rs′   = primitivize(Rs, cntr)    # primitive basis (oblique)
julia> flat′ = primitivize(flat, cntr)  # Fourier lattice in basis of Rs′
julia> plot(flat′, Rs′)
```

"""
function primitivize(flat::AbstractFourierLattice{D}, cntr::Char) where D
    # Short-circuit for lattices that have trivial transformation matrices
    (D == 3 && cntr == 'P') && return flat
    (D == 2 && cntr == 'p') && return flat
    D == 1 && return flat

    # The orbits consist of G-vector specified as a coordinate vector 𝐤≡(k₁,k₂,k₃)ᵀ, referred
    # to the conventional 𝐆-basis (𝐚* 𝐛* 𝐜*), and we want to instead express it as a coordinate
    # vector 𝐤′≡(k₁′,k₂′,k₃′)ᵀ referred to the primitive 𝐆-basis (𝐚*′ 𝐛*′ 𝐜*′)≡(𝐚* 𝐛* 𝐜*)(P⁻¹)ᵀ,
    # where P is the transformation matrix. This is achieved by transforming according to 𝐤′ = Pᵀ𝐤
    # or, equivalently, (k₁′ k₂′ k₃′)ᵀ = Pᵀ(k₁ k₂ k₃)ᵀ. See also `primitivize(::KVec)` and 
    # `primitivize(::ReciprocalBasis)`.
    P = primitivebasismatrix(cntr, D)

    orbits = getorbits(flat) # vec of vec of G-vectors (in a **conventional** 𝐆-basis)
    orbits′ = [[SVector{D, Int}(ntuple(_->0,D)) for j in eachindex(orb)] for orb in orbits] # prealloc. a vec of vec of G-vecs (to be filled in the **primitive** 𝐆-basis)
    for (i, orb) in enumerate(orbits)
        for (j, k) in enumerate(orb)
            # Note that, because the primitive reciprocal basis Gs′≡(𝐚*′ 𝐛*′ 𝐜*′) is "larger"
            # vectors than the conventional basis Gs≡(𝐚* 𝐛* 𝐜*) (since the direct lattice shrinks
            # when we go to a primitive basis), not every conventional reciprocal lattice 
            # coordinate vector 𝐤 has a primitive integer-coordinate vector 𝐤′=Pᵀ𝐤 (i.e. kᵢ∈ℕ does 
            # not imply kᵢ′∈ℕ). However, since `flat` is derived consistent with the symmetries 
            # in a conventional basis, the necessary restrictions will already have been imposed
            # in the creation of `flat` so that the primivized version will have only integer
            # coefficients (otherwise the lattice would not be periodic in the primitive cell).
            orbits′[i][j] = convert(SVector{D, Int}, P'*k)
        end
    end

    # the coefficients of flat are unchanged; only the 𝐑- and 𝐆-basis change
    return typeof(flat)(orbits′, deepcopy(getcoefs(flat))) # return in the same type as `flat`
end

"""
    modulate(flat::UnityFourierLattice{dim},
    modulation::AbstractVector{ComplexF64}=rand(ComplexF64, length(getcoefs(flat))),
    expon::Union{Nothing, Real}=nothing)
                            --> ModulatedFourierLattice

Derive a concrete, modulated Fourier lattice from `flat`, a UnityFourierLattice 
struct (that contains the _interrelations_ between orbit coefficients), by 
multiplying the "normalized" orbit coefficients by a `modulation`, a _complex_
modulating vector (in general, should be complex; otherwise restores unintended
symmetry to the lattice). Distinct `modulation` vectors produce distinct 
realizations of the same lattice described by the original `flat`. By default,
a random complex vector is used.

An exponent `expon` can be provided, which introduces a penalty term to short-
wavelength features (i.e. high-|G| orbits) by dividing the orbit coefficients
by |G|^`expon`; producing a more "localized" and "smooth" lattice boundary
when `expon > 0` (reverse for `expon < 0`). This basically amounts to a 
continuous "simplifying" operation on the lattice (it is not necessarily a 
smoothing operation; it simply suppresses "high-frequency" components).
If `expon = nothing`, no rescaling is performed. 

The `normscale(!)` methods exists to perform subsequent `expon` norm-rescaling 
of a `ModulatedFourierLattice`.
"""
function modulate(flat::UnityFourierLattice,
                  modulation::AbstractVector{ComplexF64}=rand(ComplexF64, length(getcoefs(flat))),
                  expon::Union{Nothing, Real}=nothing)

    orbits = getorbits(flat); orbitcoefs = getcoefs(flat); # unpacking ...
    
    # `expon ≠ 0` is provided, we will interpret it as a penalty term on 
    # short-wavelength orbits (i.e., high |𝐆|) by dividing the orbit 
    # coefficients by |𝐆|ᵉˣᵖᴼⁿ; this produces more "localized" and "smooth"
    # lattice boundaries for `expon > 0` (reverse for `expon < 0`).
    if !isnothing(expon) && !iszero(expon) 
        @inbounds for i in 2:length(orbits) # leaves the constant term untouched 
                                            # (there will _always_ be a constant term)...
            modulation[i] /= (norm(first(orbits[i])))^expon
        end
    end

    # scale the orbit coefficients by the overall `modulation` vector
    modulated_orbitcoefs = orbitcoefs.*modulation

    return ModulatedFourierLattice(orbits, modulated_orbitcoefs)
end

""" 
    normscale(flat::ModulatedFourierLattice, expon::Real) --> ModulatedFourierLattice

Applies subsequent norm-rescaling via `expon`; see detailed description 
in `modulate`. An in-place variant is provided as `normscale!`.
"""
normscale(flat::ModulatedFourierLattice, expon::Real) = normscale!(deepcopy(flat), expon)
"""
    normscale!(flat::ModulatedFourierLattice, expon::Real) --> ModulatedFourierLattice

In-place equivalent of `normscale`: changes `flat`.
"""
function normscale!(flat::ModulatedFourierLattice, expon::Real)
    if !iszero(expon)
        @inbounds for i in 2:length(getorbits(flat))
            rescale_factor = norm(first(getorbits(flat)[i]))^expon
            flat.orbitcoefs[i] ./= rescale_factor
        end
    end
    return flat
end

""" 
    calcfourier(xyz, flat::AbstractFourierLattice) --> Float64

Compute the real part of the function evaluation of `flat` at a
point `xyz` (a tuple, SVector, or a vector), i.e. return
    Re[∑ᵢ cᵢexp(2πi𝐆ᵢ⋅𝐫)]
with 𝐆ᵢ denoting a 𝐆-vector in an allowed orbit in `flat`, and 
cᵢ an associated coefficient (and with 𝐫 ≡ `xyz`).
"""
calcfourier(xyz, flat::AbstractFourierLattice) = calcfourier(xyz, getorbits(flat), getcoefs(flat))
function calcfourier(xyz, orbits, orbitcoefs)
    f = zero(Float64)
    for (orb, coefs) in zip(orbits, orbitcoefs)
        for (G, c) in zip(orb, coefs)
            # though one might naively think the phase would need a conversion between 
            # 𝐑- and 𝐆-bases, this is not necessary since P(𝐆)ᵀP(𝐑) = 2π𝐈 by definition
            exp_im, exp_re = sincos(2π*dot(G, xyz))
            f += real(c)*exp_re - imag(c)*exp_im    # ≡ real(exp(2π*1im*dot(G, xyz)))
        end
    end
    return f
end

"""
    plot(flat::AbstractFourierLattice, Rs::DirectBasis)

Plots a lattice `flat::AbstractFourierLattice` with lattice vectors
given by `Rs::DirectBasis`. Possible kwargs are (default in brackets) 

- `N`: resolution [`100`]
- `filling`: determine isovalue from relative filling fraction [`0.5`]
- `isoval`: isovalue [nothing (inferred from `filling`)]
- `repeat`: if not `nothing`, repeats the unit cell an integer number of times [`nothing`]
- `fig`: figure handle to plot [`nothing`, i.e. opens a new figure]

If both `filling` and `isoval` kwargs simultaneously not equal 
to `nothing`, then `isoval` takes precedence.
"""
function plot(flat::AbstractFourierLattice, Rs::DirectBasis;
              N::Integer=100, 
              filling::Union{Real, Nothing}=0.5, 
              isoval::Union{Real, Nothing}=nothing,
              repeat::Union{Integer, Nothing}=nothing,
              fig=nothing)
 
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    if isnothing(isoval)
        isoval = !isnothing(filling) ? quantile(Iterators.flatten(vals), filling) : zero(Float64)
    end
    plotiso(xyz,vals,isoval,Rs,repeat,fig)

    return xyz,vals,isoval
end


function calcfouriergridded!(vals, xyz, flat::AbstractFourierLattice, 
                             N::Integer=length(xyz))
    f = (coords...)-> calcfourier(coords, flat)
    # evaluate f over all gridpoints via broadcasting
    if dim(flat) == 2
        broadcast!(f, vals, reshape(xyz, (1,N)), reshape(xyz, (N,1)))
    elseif dim(flat) == 3
        # VERIFY: unclear if this leads to the right ordering of vals wrt x,y,z and plotting packages
        broadcast!(f, vals, reshape(xyz, (N,1,1)), reshape(xyz, (1,N,1)), reshape(xyz, (1,1,N)))
    end
    return vals
end
function calcfouriergridded(xyz, flat::AbstractFourierLattice{D},
                            N::Integer=length(xyz)) where D
    vals = Array{Float64, D}(undef, ntuple(i->N, D)...)
    return calcfouriergridded!(vals, xyz, flat, N)
end


ivec(i,dim) = begin v=zeros(dim); v[i] = 1.0; return v end # helper function
# show isocontour of data
function plotiso(xyz, vals, isoval::Real=0.0, 
                 Rs=ntuple(i->ivec(i,length(ndims(vals))), length(ndims(vals))),
                 repeat::Union{Integer, Nothing}=nothing, 
                 fig=nothing)  
    dim = ndims(vals)
    if dim == 2
        # convert to a cartesian coordinate system rather than direct basis of Ri
        N = length(xyz) 
        X = broadcast((x,y) -> x*Rs[1][1] + y*Rs[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        Y = broadcast((x,y) -> x*Rs[1][2] + y*Rs[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        uc = [[0 0]; Rs[1]'; (Rs[1]+Rs[2])'; (Rs[2])'; [0 0]] .- (Rs[1]+Rs[2])'./2
        pad = abs((-)(extrema(uc)...))/25

        if isnothing(fig)
            fig = plt.figure()
        else
            fig.clf()
        end
        fig.gca().contourf(X,Y,vals; 
                           levels=[-1e12, isoval, 1e12],
                           cmap=plt.get_cmap("gray",2))# is also good
        fig.gca().contour(X,Y,vals,levels=[isoval], colors="w", linestyles="solid")
        fig.gca().plot(uc[:,1], uc[:,2], color="C4",linestyle="solid")
        fig.gca().scatter([0],[0],color="C4",s=30, marker="+")
        

        if !isnothing(repeat) # allow repetitions of unit cell in 2D
            for r1 in -repeat:repeat
                for r2 in -repeat:repeat
                    if r1 == r2 == 0; continue; end
                    offset = Rs[1].*r1 .+ Rs[2].*r2
                    fig.gca().contourf(X.+offset[1],Y.+offset[2],vals,levels=[minimum(vals), isoval, maximum(vals)]; cmap=plt.get_cmap("gray",256)) #get_cmap(coolwarm,3) is also good
                    fig.gca().contour(X.+offset[1],Y.+offset[2],vals,levels=[isoval], colors="w", linestyles="solid")
                end
            end
            xd = -(-)(extrema(uc[:,1])...); yd = -(-)(extrema(uc[:,2])...)
            plt.xlim([extrema(uc[:,1])...].+[-1,1].*repeat*xd.+[-1,1].*pad); 
            plt.ylim([extrema(uc[:,2])...].+[-1,1].*repeat*yd.+[-1,1].*pad);
        else
            plt.xlim([extrema(uc[:,1])...].+[-1,1].*pad); plt.ylim([extrema(uc[:,2])...].+[-1,1].*pad);
        end
        fig.gca().set_aspect("equal", adjustable="box")
        fig.gca().set_axis_off()

    elseif dim == 3
        scene=Scene()
        Makie.contour!(scene, xyz,xyz,xyz, vals,
                       levels=[isoval],colormap=:blues, linewidth=.1)
        Makie.display(scene)

        # marching cubes algorithm to find isosurfaces
        algo = MarchingCubes(iso=isoval, eps=1e-3)
        verts, faces = isosurface(vals, algo; 
                                  origin = SVector(-0.5,-0.5,-0.5), 
                                  widths = SVector(1.0,1.0,1.0))
        verts′ = [verts[i][j] for i = 1:length(verts), j = 1:3]
        faces′ = [faces[i][j] for i = 1:length(faces), j = 1:3]

        #println("Mesh: $(length(verts)) vertices\n", " "^6, "$(length(faces)) faces")
        isomesh = convert_arguments(Mesh, verts′, faces′)[1]

        # plot isosurface
        scene = Scene()
        mesh!(isomesh, color=:grey)
        display(scene)
    end
    return nothing
end