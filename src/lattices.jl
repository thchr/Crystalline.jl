

_boolify(x, tol=0.5) = x > tol ? one(x) : zero(x)

function symmetrize!(sg::SpaceGroup, dat::Array{T} where T<:Number) #only for cubic lattices for now
    #if crystalsystem(C) != "cubic"; error("Not yet implemented"); end
    if !issymmorph(sg);       error("Not yet implemented"); end
    opdat = similar(dat)
    for op in operations(sg)
        transform!(opdat, op, dat)
        dat .+= opdat
        dat ./= 2
    end
    return dat
end

symmetrize!(sgnum::Integer, dim::Integer, dat::Array) = symmetrize!(get_symops(sgnum,dim,verbose=false), dat)

function transform(op::SymOperation, dat::Array{Number})
    return transform!(similar(op), op, dat)
end

# this should probably be done using something like IndirectArrays.jl to achieve appropriate generality
# the other option would be to use interpolation, which might be more appropriate for getting 
# general fractional translations right.
function transform!(opdat::Array{T}, op::SymOperation, dat::Array{T}) where T<:Number
    @assert axes(opdat) == axes(dat)

    perms = (findall(.!iszero.(oprow)) for oprow in eachrow(rotation(op)))
    display(op)
    display(collect(perms))
    perm = findfirst.(eachrow(.!iszero.(rotation(op))))

    display(perm)
    permutedims!(opdat, dat, perm)

    for i=1:dim(op)
        invert = op[i,perm[i]]<0
        if invert
            opdat = reverse(opdat, dims=i)
        end
    end

    return opdat
end




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

    # we define the "reciprocal orbit" associated with the action of W through (W⁻¹)ᵀ
    # Calculates the operators (W⁻¹)ᵀ in the G-basis:
    # The action of a symmetry operator in an 𝐑-basis, i.e. W(𝐑), on a 𝐤 vector in a 
    # 𝐆-basis, i.e. 𝐤(𝐆), is 𝐤′(𝐆)ᵀ = 𝐤(𝐆)ᵀW(𝐑)⁻¹. To deal with column vectors, we 
    # transpose, obtaining 𝐤′(𝐆) = [W(𝐑)⁻¹]ᵀ𝐤(𝐆) [details in symops.jl, above littlegroup(...)].
    W⁻¹ᵀs = transpose.(inv.(Ws))

    # if idxmax is interpreted as (imax, jmax, ...), then this produces an iterator
    # over i = -imax:imax, j = -jmax:jmax, ..., where each call returns (..., j, i); 
    # note that the final order is anti-lexicographical; so we reverse it in the actual
    # loop for our own sanity's sake
    reviter = Iterators.product(reverse((:).(.-idxmax, idxmax))...)

    orbits = Vector{Vector{SVector{dim,Int64}}}() # vector to store orbits of G-indices into G-basis
    for rG in reviter  # --- compute orbits ---
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
                conds[n,m] = exp(-1im*2π*dot(G′, w)) 
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
    perm = sortperm(orbits, by=x->norm(x[1]))
    permute!(orbits, perm)
    permute!(orbitcoefs, perm)

    return FourierLattice(orbits, orbitcoefs)
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

calcfourier(xyz, flat::FourierLattice) = calcfourier(xyz, flat.orbits, flat.orbitcoefs)
function calcfourier(xyz, orbits, orbitcoefs)
    f = zero(ComplexF64)
    for (orb, coefs) in Iterators.zip(orbits, orbitcoefs)
        for (G, c) in Iterators.zip(orb, coefs)
            # though one might naively think the phase would need a conversion between 
            # 𝐑- and 𝐆-bases, this is not necessary since P(𝐆)ᵀP(𝐑) = 2π𝐈 by definition
            f += c*exp(1im*2π*dot(G,xyz))
        end
    end
    return f
end

function plotfourier(flat::FourierLattice{dim}, C::Crystal, N=100, expon=2, filling=0.5, repeat=nothing) where dim
    orbits = flat.orbits; orbitcoefs = flat.orbitcoefs; R = basis(C) # unpacking ...
    xyz = range(-.5, .5, length=N)
    modulated_orbitcoefs = orbitcoefs.*rand(ComplexF64, length(orbitcoefs)) # in general, a complex value; otherwise we restore unintended symmetry
    if !isnothing(expon) && !iszero(expon) # divide shorter wavelength terms by their (multiplicity*norm)^expon (to get a more "localized" and smooth character)
        for i = 2:length(orbits) # leave the constant term untouched ...
        modulated_orbitcoefs[i] ./= (norm(orbits[i])).^expon
        end
    end
    
    vals = Array{Float64, dim}(undef, ntuple(i->N, dim)...)
    calcfouriergridded!(vals, xyz, orbits, modulated_orbitcoefs, dim, N)

    if !isnothing(filling)
        isoval = quantile(Iterators.flatten(vals), filling)
    else
        isoval = 0
    end
    
    plotiso(xyz,vals,isoval,R,repeat)

    return xyz,vals,isoval
end


function calcfouriergridded!(vals, xyz, orbits, orbitcoefs, dim, N)
    f = (coords...)-> real(calcfourier(coords, orbits, orbitcoefs))
    # evaluate f over all gridpoints via broadcasting
    if dim == 2
        broadcast!(f, vals, reshape(xyz, (1,N)), reshape(xyz, (N,1)))
    elseif dim == 3
        # unclear if this leads to the right ordering of vals wrt x,y,z and plotting packages
        broadcast!(f, vals, reshape(xyz, (N,1,1)), reshape(xyz, (1,N,1)), reshape(xyz, (1,1,N)))
    end
    return vals
end
function calcfouriergridded(xyz, orbits, orbitcoefs)
    dim = length(orbits[1][1]); N = length(xyz)
    vals = Array{Float64, dim}(undef, ntuple(i->N, dim)...)
    return calcfouriergridded!(vals, xyz, orbits, orbitcoefs, dim, N)
end


ivec(i,dim) = begin v=zeros(dim); v[i] = 1.0; return v end # helper function
# show isocontour of data
function plotiso(xyz, vals, isoval=0, 
                 R=ntuple(i->ivec(i,length(ndims(vals))), length(ndims(vals))),
                 repeat=nothing)  
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