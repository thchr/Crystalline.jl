function gen_lattice(sgnum::Integer, dim::Integer=2, 
                     Ns::Tuple=ntuple(x->10, dim))
    sg = get_symops(sgnum, dim)
    C = gen_crystal(sgnum, dim)

    ϵ = rand(Float64, Ns...)
    ϵ .= boolify.(ϵ, 0.85)
    
    symmetrize!(sg, ϵ)
end

boolify(x, tol=0.5) = x > tol ? one(x) : zero(x)


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


function plotlattice(dat)
    plt.close("all")
    dim = ndims(dat)
    fig=plt.figure()
    if dim == 3
        ax = fig.gca(projection="3d")
        ax.voxels(dat)
    elseif dim == 2
        fig.gca().pcolor(boolify.(dat, 0.15))
    end
end