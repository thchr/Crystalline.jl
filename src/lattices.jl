

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

    perms = (findall(.!iszero.(oprow)) for oprow in eachrow(pg(op)))
    display(op)
    display(collect(perms))
    perm = findfirst.(eachrow(.!iszero.(pg(op))))

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





function levelsetlattice(sgnum::Int64, dim::Int64=2, 
                         idxmax::NTuple=ntuple(i->2,dim))
    sg = get_symops(sgnum, dim)
    C = gen_crystal(sgnum, dim)
    R = basis(C)
    G = reciprocalbasis(R)

    symops = operations(sg)
    Wops = pg.(symops) # operations W in R-basis (point group part)
    PR = hcat(R...) # matrix with cols of R[i]; rCartesian = PR*rDirect
    PG = hcat(G...) # matrix with cols of G[i]; gCartesian = PG*gDirect
    
    #= # sym ops in a direct G basis (e.g., [1,1,0] means G₁-G₂); basis change 
       # via W' = PG⁻¹*PR*W*PR⁻¹*PG = iPGPR*W*iPGPR⁻¹
       iPGPR = PG\PR # PG⁻¹*PR
       Wops_Gbasis = [iPGPR*op/iPGPR for op in Wops]  =#
    PGᵀPR = PG'*PR  # (⋆)  TODO: this is actually just 2πI since PGᵀ = 2πPR⁻¹ by definition ...  
                    # (⋆⋆) TODO: Actually, this may suggest a deeper problem; see the TODO before littlegroup(..) in bravais.jl
    # Calculates the operators (W⁻¹)ᵀ in the G-basis; note that although
    # the transformations are orthogonal in the Cartesian basis (i.e. W⁻¹=Wᵀ),
    # this is not generally the case in other bases. Here, what we are really
    # manipulating is (this is what we need to be invariant essentially)
    #   exp[idot(G, W⁻¹r)] = exp(iGᵀW⁻¹r) = exp{i[(W⁻¹)ᵀG]ᵀr}
    # so we define the "reciprocal orbit" associated with the action of W through (W⁻¹)ᵀ
    W⁻¹ᵀops_Gbasis = [(PGᵀPR*op)/PGᵀPR for op in Wops] 
    # TODO: due to (⋆), we actually end up getting W⁻¹ᵀops_Gbasis = Wops, so this could 
    # all be removed and replaced by appropriate comments... Requires that (⋆⋆) isn't a concern though.
    
    # if idxmax is interpreted as (imax, jmax, ...), then this produces an iterator
    # over i = -imax:imax, j = -jmax:jmax, ..., where each call returns (..., j, i); 
    # note that the final order is anti-lexicographical; so we reverse it in the actual
    # loop for our own sanity's sake
    reviter = Iterators.product(reverse((:).(.-idxmax, idxmax))...)

    ijkorbits = Vector{Vector{SVector{dim,Int64}}}() # vector to store orbits of ijk-indices into G-basis
    for rijk in reviter  # --- compute orbits ---
        ijk = SVector{dim,Int64}(reverse(rijk)) # fix order and convert to SVector{dim,Int64} from Tuple

        skip = false # if ijk already contained in an orbit; go to next ijk
        for orb in ijkorbits
            isapproxin(ijk, orb) && (skip=true; break) 
        end
        skip && continue
        
        neworb = orbit(W⁻¹ᵀops_Gbasis, ijk) # compute orbit assoc with ijk-combination
        # the symmetry transformation may introduce round-off errors, but we know that 
        # the indices must be integers; fix that here, and check its validity as well
        neworb′ = [round.(Int64,ijk′) for ijk′ in neworb] 
        if norm(neworb′ .- neworb) > 1e-10; 
            error("The ijk-combinations and their symmetry-transforms must be integers"); 
        end
        push!(ijkorbits, neworb′) # add orbit to list of orbits
    end

    # compute restrictions on orbit coefficients due to any nonsymmorphic elements
    # in the space group
    #if !issymmorph(sg)
    wops = translation.(symops)
    # calculate inverse translation: this is equal to -W⁻¹w; we do it first 
    # in the direct R basis, and then transform to Cartesian basis
    W⁻¹w = [PR*(W\w) for (W,w) = zip(Wops, wops)] 

    orbcoefs = Vector{Vector{ComplexF64}}()
    deleteidx = Vector{Int64}()
    for (o,orb) in enumerate(ijkorbits)
        start = true; prevspan = []
        for (Wop, wop) in zip(W⁻¹ᵀops_Gbasis, W⁻¹w)
            conds = zeros(ComplexF64, length(orb), length(orb))
            for (m, ijk) in enumerate(orb)
                ijk′ = Wop*ijk  # where the ijk is transformed to by Wop
                diffs = norm.(Ref(ijk′) .- orb); 
                n = argmin(diffs) # find assoc linear index in orbit
                diffs[n] > 1e-10 && error("Part of an orbit was miscalculated")
                conds[n,m] = exp(-1im*dot(PG*ijk, wop)) 
            end

            nextspan = nullspace(conds-I, atol=1e-12)
            if start 
                prevspan = nextspan
                start = false
            elseif !isempty(prevspan) && !isempty(nextspan)
                spansect = nullspace([prevspan -nextspan], atol=1e-12)[size(prevspan, 2)+1:end,:]
                prevspan = nextspan*spansect
            else
                prevspan = nothing; break
            end
        end
                    
        if !isnothing(prevspan)
            if size(prevspan,2) != 1; error("Unexpected size of prevspan"); end
            coefbasis = vec(prevspan)
            coefbasis ./= coefbasis[argmax(norm(coefbasis, Inf))]
            push!(orbcoefs, coefbasis)
        else 
            push!(deleteidx, o)
        end
    end

    ijkorbits_forbidden = ijkorbits[deleteidx]
    deleteat!(ijkorbits, deleteidx)

    # sort in order of descending wavelength (e.g., [0,0,...] term comes first; highest ijk-combinations come last)
    perm = sortperm(ijkorbits, by=x->norm(x[1]))
    permute!(ijkorbits, perm)
    permute!(orbcoefs, perm)

    return ijkorbits, orbcoefs, R
end


"""
    orbit(Wops, x)

    Computes the orbit of `x` under a set of point-group operations `Wops`,
    i.e. computes the set `{gx | g∈G}` where `g` denotes elements of the group
    `G` composed of all operations in `Wops` (possibly iterated, to ensure
    full coverage).
    At the moment, we only consider _point group_ operations; i.e. there are 
    no nonsymmorphic `wops` parts. 
    It is important that `Wops` and `x` are given in the same basis. 
    [W' = PWP⁻¹ if the basis change is from coordinates r to r' = Pr, corresponding 
    to a new set of basis vectors (x̂')ᵀ=x̂ᵀP; e.g., when going from a direct basis
    representation to a Cartesian one, the basis change matrix is P = [R₁ R₂ R₃],
    with Rᵢ inserted as column vectors]
"""
function orbit(Wops, x)
    fx = float.(x)
    xorbit = [fx]
    for W in Wops
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

function isapproxin(x, itr)
    for y in itr
        if isapprox(x, y)
            return true
        end
    end
    return false
end


function calcfourier(xyz, ijkorbits, orbcoefs)
    f = zero(ComplexF64)
    for (orb, coefs) in Iterators.zip(ijkorbits, orbcoefs)
        for (ijk, c) in Iterators.zip(orb,coefs)
            # though one might naively think the phase would need a conversion between 
            # R and G bases, this is not necessary since PG'*PR = 2πI by definition
            f += c*exp(1im*2*pi*dot(ijk,xyz))
        end
    end
    return f
end

function plotfourier(ijkorbits, orbcoefs, R, N=100, expon=2, filling=0.5)
    dim = length(ijkorbits[end][end])
    xyz = range(-.5, .5, length=N)
    modulated_orbcoefs = orbcoefs.*rand(ComplexF64, length(orbcoefs)) # should be a general complex value; otherwise we restore unintended symmetry
    if !isnothing(expon) && !iszero(expon) # divide shorter wavelength terms by their (multiplicity*norm)^expon (to get a more "localized" and smooth character)
        for i = 2:length(ijkorbits) # leave the constant term untouched ...
        modulated_orbcoefs[i] ./= (norm(ijkorbits[i])).^expon
        end
    end
    
    vals = Array{Float64, dim}(undef, ntuple(i->N, dim)...)
    calcfouriergridded!(vals, xyz, ijkorbits, modulated_orbcoefs, dim, N)

    if !isnothing(filling)
        isoval = quantile(Iterators.flatten(vals), filling)
    else
        isoval = 0
    end
    
    plotiso(xyz,vals,isoval,R)

    return xyz,vals,isoval
end


function calcfouriergridded!(vals, xyz, ijkorbits, orbcoefs, dim, N)
    f = (coords...)-> real(calcfourier(coords, ijkorbits, orbcoefs))
    # evaluate f over all gridpoints via broadcasting
    if dim == 2
        broadcast!(f, vals, reshape(xyz, (1,N)), reshape(xyz, (N,1)))
    elseif dim == 3
        # unclear if this leads to the right ordering of vals wrt x,y,z and plotting packages
        broadcast!(f, vals, reshape(xyz, (N,1,1)), reshape(xyz, (1,N,1)), reshape(xyz, (1,1,N)))
    end
    return vals
end
function calcfouriergridded(xyz, ijkorbits, orbcoefs)
    dim = length(ijkorbits[1][1]); N = length(ijkorbits[1][1])
    vals = Array{Float64, dim}(undef, ntuple(i->N, dim)...)
    return calcfouriergridded!(vals, xyz, ijkorbits, orbcoefs, dim, N)
end


ivec(i,dim) = begin v=zeros(dim); v[i] = 1.0; return v end # helper function
# show isocontour of data
function plotiso(xyz, vals, isoval=0, R=ntuple(i->ivec(i,length(ndims)), length(ndims)))  
    dim = ndims(vals)
    if dim == 2
        # convert to a cartesian coordinate system rather than direct basis of Ri
        N = length(xyz) 
        X = broadcast((x,y) -> x*R[1][1] + y*R[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        Y = broadcast((x,y) -> x*R[1][2] + y*R[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
        uc = [[0 0]; R[1]'; (R[1]+R[2])'; (R[2])'; [0 0]] .- (R[1]+R[2])'./2
        pad = abs((-)(extrema(uc)...))/25

        plt.close("all")
        fig = plt.figure()
        fig.gca().contourf(X,Y,vals,levels=[minimum(vals), isoval, maximum(vals)]; cmap=plt.get_cmap("gray",256)) #get_cmap(coolwarm,3) is also good
        fig.gca().contour(X,Y,vals,levels=[isoval], colors="w", linestyles="solid")
        fig.gca().plot(uc[:,1], uc[:,2], color="C4",linestyle="solid")
        fig.gca().scatter([0],[0],color="C4",s=30, marker="+")
        plt.xlim([extrema(uc[:,1])...].+[-1,1].*pad); plt.ylim([extrema(uc[:,2])...].+[-1,1].*pad);
        fig.gca().set_aspect("equal", adjustable="box")
        fig.gca().set_axis_off()
    elseif dim == 3
        @show sum(vals.-isoval .> 0)/length(vals)
        scene=Scene()
        Makie.contour!(scene, xyz,xyz,xyz,vals,levels=[isoval],color=:blue,fillrange=true)
        Makie.display(scene)
    end
    return nothing
end
# fig = plt.figure()
# fig.gca().pcolor(real.(plotfourier(levelsetlattice(17,2,(3,3))[1])))