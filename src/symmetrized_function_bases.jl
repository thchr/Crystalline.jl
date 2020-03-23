using SGOps, StaticArrays, LinearAlgebra, SparseArrays, PyPlot

struct CrystalPW{D}
    Gs::Vector{SVector{D,Int}}
    cs::Vector{ComplexF64}
end

const DEFAULT_αβγ = [0,0,0]
# This essentially follows the discussion in Blokker, "Symmetry Projection of Crystal Wave 
# Functions by Means of a Computer", Journal of Computational Physics 12, 471 (1973), 
# [https://doi.org/10.1016/0021-9991(73)90099-5], especially §4 and Eq. (4.13) and below.
function planewave_projector(sgnum::Integer, klab::String, Dᵛ::Val{D}, idxmax::NTuple{D,Int}) where D
    lgirs = find_lgirreps(sgnum, klab, Dᵛ)
    lg = group(first(lgirs))
    kv = kvec(lg)(DEFAULT_αβγ)

    # Reciprocal lattice vectors 𝐊 in a basis of 𝐆ᵢ
    Ks = SVector{D, Int}.(Tuple.(CartesianIndices(map(i->-i:i,idxmax))))
    idxs = LinearIndices(Ks)
    s = length(Ks)

    Γs = Vector{SparseMatrixCSC{ComplexF64,Int}}(undef, length(lg)) # phases of SymOperation action on PWs
    M = Vector{Int64}(undef, s); V = Vector{ComplexF64}(undef, s)
    for (i, opᵢ) in enumerate(lg)
        for (n, Kₙ) in enumerate(Ks)
            Kₙⁱ = rotation(opᵢ)'\Kₙ # see also compose(::SymOperation, ::KVec)

            # TODO: This will not work for trigonal/hexagonal lattices where Kₙⁱ may not be
            #       contained in Ks (rotates "outside the square box"): impacts sgs 143-194
            #       More generally, it seems better to just start from some set of Ks and 
            #       then dynamically expand the set until we cover sufficient ground
            M[n] = idxs[findfirst(==(Kₙⁱ), Ks)]
            V[n] = cis(-2π*dot(kv+Kₙⁱ, translation(opᵢ)))
        end
        Γs[i] = sparse(M, eachindex(Ks), V, s, s)
    end

    # construct projection matrices
    Ps = [[spzeros(ComplexF64, s,s) for _ in Base.OneTo(SGOps.irdim(lgir))] for lgir in lgirs]
    for (j, lgir) in enumerate(lgirs)
        ʲΓ = irreps(lgir, DEFAULT_αβγ)
        dimʲΓ = SGOps.irdim(lgir)
        for d in Base.OneTo(dimʲΓ)
            for (i, opᵢ) in enumerate(lg)
                Ps[j][d] .+= conj(ʲΓ[i][d,d])*Γs[i]
            end
            Ps[j][d] .*= dimʲΓ/length(lg)
        end
    end

    return Ks, Ps
end


function planewave_symfuncs(sgnum::Integer, klab::String, Dᵛ::Val{D}, idxmax::NTuple{D,Int}) where D
    Ks, Ps = planewave_projector(sgnum, klab, Dᵛ, idxmax)
    s = prod(i->2i+1, idxmax)
    PWs = [[CrystalPW{D}[] for d in 1:size(Ps[j],1)] for j in eachindex(Ps)]
    for j in eachindex(Ps)
        for d in Base.OneTo(size(Ps[j],1))
            # This effectively calculates the _range_ of the projection matrix; i.e. its 
            # column space. We include only those elements with singular values equal to 1.
            # (see https://stackoverflow.com/a/43267856/9911781)
            U, σs, _ = svd(convert(Matrix{ComplexF64},Ps[j][d]))
            for n in 1:size(U,2)
                if isapprox(σs[n], 1.0, atol = 1e-11)
                    keep = filter(i->abs(U[i,n])>1e-11, 1:size(U,1))
                    push!(PWs[j][d], CrystalPW{D}(Ks[keep], U[keep,n]))
                else
                    break # use that σs are sorted
                end
            end
        end
    end
        

    # We ought to reorthogonalize PWs and create linear combinations that "sort" the distinct
    # PWs within a specific [j][d] combination after increasing mean(Ks) contributions
    return PWs

end


"""

Example:

    sgnum = 220
    PWs=planewave_symfuncs(sgnum,"Γ", Val(3), 2 .* (1,1,1));
    j = 3; d = 1; n = 1
    plot_xyslice(PWs[j][d][n], directbasis(sgnum))
"""
function plot_xyslice(PW::CrystalPW{D}, Rs::DirectBasis{D}, z::Real=0, N::Integer=100) where D
    f = (coords...)-> calcfourier(SVector{3,Float64}(coords..., z), PW)
    vals = Matrix{ComplexF64}(undef, N, N)
    xyz = range(-.5, .5, length=N)

    # TODO: Note that right now, we do not include the Bloch phase, i.e. the result is the
    #       periodic envelope u rather than ψ. Maybe would be good with a toggle.
    broadcast!(f, vals, reshape(xyz, (N,1)), reshape(xyz, (1,N)))

    fig = plt.figure()
    
    X = broadcast((x,y) -> x*Rs[1][1] + y*Rs[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
    Y = broadcast((x,y) -> x*Rs[1][2] + y*Rs[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
    
    fig.gca().contourf(X,Y,real.(vals); levels=25, cmap=plt.get_cmap("coolwarm",25))
end

function calcfourier(xyz, PW::CrystalPW{D}) where D
    f = zero(ComplexF64)
    for (G, c) in zip(PW.Gs, PW.cs)
        f += c*cis(2π*dot(G, xyz)) # ≡ exp(2πiGᵀr))
    end
    return f
end