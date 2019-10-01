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
    orbits::Vector{Vector{SVector{D, Int64}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}    # Vector of interrelations between coefficients of 𝐆-plane waves within an orbit; unit norm
end
UnityFourierLattice(orbits, orbitcoefs) = begin
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
    orbits::Vector{Vector{SVector{D, Int64}}} # Vector of orbits of 𝐆-vectors (in 𝐆-basis)
    orbitcoefs::Vector{Vector{ComplexF64}}    # Vector of coefficients of 𝐆-plane waves within an orbit
end


## --- METHODS --- 

# Group orbits of plane waves G = (G)ᵀ under a symmetry operation Ô = {W|w}, 
# using that Ô acts as Ô⁻¹={W⁻¹|-W⁻¹w} when acting on functions, i.e.
#   Ôexp(iG⋅r) = Ôexp(iG⋅Ô⁻¹r) = exp[iG⋅(W⁻¹r-W⁻¹w)]
# and 
#   exp(iG⋅W⁻¹r) = exp(iGᵀW⁻¹r) = exp{i[(W⁻¹)ᵀG]ᵀ⋅r}
function levelsetlattice(sgnum::Int64, dim::Int64=2, 
                         idxmax::NTuple=ntuple(i->2,dim))
    sg = get_symops(sgnum, dim)
    symops = operations(sg)
    Ws = rotation.(symops) # operations W in R-basis (point group part)
    ws = translation.(symops)

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
    orbits = Vector{Vector{SVector{dim,Int64}}}() # vector to store orbits of G-vectors (in G-basis)
    for rG in reviter  
        G = SVector{dim,Int64}(reverse(rG)) # fix order and convert to SVector{dim,Int64} from Tuple

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
    @inbounds for i in 2:length(getorbits(flat))
        rescale_factor = norm(first(getorbits(flat)[i]))^expon
        flat.orbitcoefs[i] ./= rescale_factor
    end
    return flat
end

calcfourier(xyz, flat::AbstractFourierLattice) = calcfourier(xyz, getorbits(flat), getcoefs(flat))
function calcfourier(xyz, orbits, orbitcoefs)
    f = zero(ComplexF64)
    for (orb, coefs) in Iterators.zip(orbits, orbitcoefs)
        for (G, c) in Iterators.zip(orb, coefs)
            # though one might naively think the phase would need a conversion between 
            # 𝐑- and 𝐆-bases, this is not necessary since P(𝐆)ᵀP(𝐑) = 2π𝐈 by definition
            f += c*cis(2π*dot(G, xyz)) # cis(x) = exp(ix)
        end
    end
    return f
end

function plotfourier(flat::AbstractFourierLattice, 
                     C::Crystal, N::Integer=100, 
                     filling::Union{Real, Nothing}=0.5, 
                     repeat::Union{Integer, Nothing}=nothing)
 
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    isoval = !isnothing(filling) ? quantile(Iterators.flatten(vals), filling) : zero(Float64)
    
    plotiso(xyz,vals,isoval,basis(C),repeat)

    return xyz,vals,isoval
end


function calcfouriergridded!(vals, xyz, flat::AbstractFourierLattice, 
                             N::Integer=length(xyz))
    f = (coords...)-> real(calcfourier(coords, flat))
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
function plotiso(xyz, vals, isoval=0, 
                 R=ntuple(i->ivec(i,length(ndims(vals))), length(ndims(vals))),
                 repeat::Union{Integer, Nothing}=nothing)  
    dim = ndims(vals)
    if dim == 2
        # convert to a cartesian coordinate system rather than direct basis of Ri
        N = length(xyz) 
        X = broadcast((x,y) -> x*R[1][1] + y*R[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        Y = broadcast((x,y) -> x*R[1][2] + y*R[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        uc = [[0 0]; R[1]'; (R[1]+R[2])'; (R[2])'; [0 0]] .- (R[1]+R[2])'./2
        pad = abs((-)(extrema(uc)...))/25

        fig = plt.figure()
        fig.gca().contourf(X,Y,vals,levels=[minimum(vals), isoval, maximum(vals)]; cmap=plt.get_cmap("gray",256)) #get_cmap(coolwarm,3) is also good
        fig.gca().contour(X,Y,vals,levels=[isoval], colors="w", linestyles="solid")
        fig.gca().plot(uc[:,1], uc[:,2], color="C4",linestyle="solid")
        fig.gca().scatter([0],[0],color="C4",s=30, marker="+")
        

        if repeat !== nothing # allow repetitions of unit cell in 2D
            for r1 in -repeat:repeat
                for r2 in -repeat:repeat
                    if r1 == r2 == 0; continue; end
                    offset = R[1].*r1 .+ R[2].*r2
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
    end
    return nothing
end