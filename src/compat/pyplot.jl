using .PyPlot
import .PyPlot: plot, plot3D, plt
using Meshing

# Defaults --------------------------------------------------------------------------------

const ORIGIN_MARKER_OPTS = (marker="o", markerfacecolor="white", markeredgecolor="black", 
                            markeredgewidth=1.5, markersize=4.5)

# ::DirectBasis ---------------------------------------------------------------------------

function plot(Rs::DirectBasis{D}, 
              cntr::SVector{D, <:Real}=zeros(SVector{D, Float64}),
              ax=plt.figure().gca(projection = D==3 ? (using3D(); "3d") : "rectilinear")) where D

    Rs′ = Rs .+ Ref(cntr) # basis vectors translated by cntr
    
    if D == 1
        ax.plot((cntr[1], Rs′[1]), (0, 0))
        ax.plot((cntr[1],), (0,); ORIGIN_MARKER_OPTS...) # origin

    elseif D == 2
        corner = sum(Rs) + cntr
        for R′ in Rs′
            ax.plot((cntr[1], R′[1]), (cntr[2], R′[2]); color="black") # basis vectors
            ax.plot((R′[1], corner[1]), (R′[2], corner[2]); color="grey") # remaining unit cell boundaries
        end
        ax.plot((cntr[1],), (cntr[2],); ORIGIN_MARKER_OPTS...) # origin
    elseif D == 3
        corners = (Rs[1]+Rs[3], Rs[1]+Rs[2], Rs[2]+Rs[3]) .+ Ref(cntr)
        dirs = ((-1,1,-1), (-1,-1,1), (1,-1,-1))
        for R′ in Rs′
            ax.plot3D((cntr[1], R′[1]), (cntr[2], R′[2]), (cntr[3], R′[3]); color="black") # basis vectors
        end
        for (i,R) in enumerate(Rs)
            for (corner,dir) in zip(corners,dirs) # remaining unit cell boundaries
                ax.plot3D((corner[1], corner[1]+dir[i]*R[1]), 
                          (corner[2], corner[2]+dir[i]*R[2]), 
                          (corner[3], corner[3]+dir[i]*R[3]); color="grey")
            end
        end
        ax.plot3D((cntr[1],), (cntr[2],), (cntr[3],); ORIGIN_MARKER_OPTS...) # origin
        ax.set_zlabel("z")
    end
    ax.set_xlabel("x"); ax.set_ylabel("y")
    # seems broken in 3D (https://github.com/matplotlib/matplotlib/pull/13474); 
    # FIXME: may raise an error in later matplotlib releases
    ax.set_aspect("equal", adjustable="box")
    return nothing
end

# ::AbstractFourierLattice ----------------------------------------------------------------

"""
    plot(flat::AbstractFourierLattice, Rs::DirectBasis)

Plots a lattice `flat::AbstractFourierLattice` with lattice vectors
given by `Rs::DirectBasis`. Possible kwargs are (defaults in brackets) 

- `N`: resolution [`100` in 2D, `20` in 3D]
- `filling`: determine isovalue from relative filling fraction [`0.5`]
- `isoval`: isovalue [nothing (inferred from `filling`)]
- `repeat`: if not `nothing`, repeats the unit cell an integer number of times [`nothing`]
- `fig`: figure handle to plot [`nothing`, i.e. opens a new figure]

If both `filling` and `isoval` kwargs simultaneously not equal 
to `nothing`, then `isoval` takes precedence.
"""
function plot(flat::AbstractFourierLattice, Rs::DirectBasis{D};
              N::Integer=(D==2 ? 100 : 20), 
              filling::Union{Real, Nothing}=0.5, 
              isoval::Union{Real, Nothing}=nothing,
              repeat::Union{Integer, Nothing}=nothing,
              fig=nothing) where D
 
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    if isnothing(isoval)
        isoval = !isnothing(filling) ? quantile(Iterators.flatten(vals), filling) : zero(Float64)
    end
    plotiso(xyz,vals,isoval,Rs,repeat,fig)

    return xyz,vals,isoval
end

# plot isocontour of data
function plotiso(xyz, vals, isoval::Real, Rs::DirectBasis{D},
                repeat::Union{Integer, Nothing}=nothing, 
                fig=nothing) where D

    if isnothing(fig)
        fig = plt.figure()
    else
        fig.clf()
    end
    ax = fig.gca(projection= D==3 ? (using3D(); "3d") : "rectilinear")

    if D == 2
    # convert to a cartesian coordinate system rather than direct basis of Ri
    N = length(xyz) 
    X = broadcast((x,y) -> x*Rs[1][1] + y*Rs[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
    Y = broadcast((x,y) -> x*Rs[1][2] + y*Rs[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))

    ax.contourf(X,Y,vals; levels=(-1e12, isoval, 1e12), cmap=plt.get_cmap("gray",2))
    ax.contour(X,Y,vals,levels=(isoval,), colors="w", linestyles="solid")
    origo = sum(Rs)./2
    plot(Rs, -origo, ax) # plot unit cell
    ax.scatter([0],[0], color="C4", s=30, marker="+")

    uc = (zeros(eltype(Rs)), Rs[1], Rs[1]+Rs[2], Rs[2]) .- Ref(origo)

    pad = (maximum(maximum.(uc)) .- minimum(minimum.(uc)))/25
    xlims = extrema(getindex.(uc, 1)); ylims = extrema(getindex.(uc, 2))
    if !isnothing(repeat) # allow repetitions of unit cell in 2D
    minval, maxval = extrema(vals)
    for r1 in -repeat:repeat
    for r2 in -repeat:repeat
        if r1 == r2 == 0; continue; end
        offset = Rs[1].*r1 .+ Rs[2].*r2
        X′ = X .+ offset[1]; Y′ = Y .+ offset[2]
        ax.contourf(X′, Y′, vals; levels=(minval, isoval, maxval),
            cmap=plt.get_cmap("gray",256)) #get_cmap(coolwarm,3) is also good
        ax.contour(X′, Y′, vals; levels=(isoval,), colors="w", linestyles="solid")
    end
    end

    xd = -(-)(xlims...)*repeat; yd = -(-)(ylims...)*repeat
    plt.xlim(xlims .+ (-1,1).*xd .+ (-1,1).*pad) 
    plt.ylim(ylims .+ (-1,1).*yd .+ (-1,1).*pad)
    else
    plt.xlim(xlims .+ (-1,1).*pad); plt.ylim(ylims .+ (-1,1).*pad);
    end
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    elseif D == 3
    # TODO: Isocaps (when Meshing.jl supports it)

    # Calculate a triangular meshing of the Fourier lattice using Meshing.jl
    verts′, faces′ = mesh_3d_levelsetlattice(vals, isoval, Rs)

    # TODO: All plot utilities should probably be implemented as PlotRecipies or similar
    # In principle, it would be good to plot the isosurface using Makie here; but it 
    # just doesn't make sense to take on Makie as a dependency, solely for this purpose.
    # It seems more appropriate to let a user bother with that, and instead give some
    # way to extract the isosurface's verts and faces. Makie can then plot it with
    #     using Makie
    #     verts′, faces′ = Crystalline.mesh_3d_levelsetlattice(flat, isoval, Rs)
    #     isomesh = convert_arguments(Mesh, verts′, faces′)[1]
    #     scene = Scene()
    #     mesh!(isomesh, color=:grey)
    #     display(scene)
    # For now though, we just use matplotlib; the performance is awful though, since it
    # vector-renders each face (of which there are _a lot_).
    plot_trisurf(verts′[:,1], verts′[:,2], verts′[:,3], triangles = faces′ .- 1)
    plot(Rs, -sum(Rs)./2, ax) # plot the boundaries of the lattice's unit cell

    end
    return nothing
end

# Meshing utilities (for ::AbstractFourierLattice) ----------------------------------------
using .Meshing

# TODO: Refactor meshing part out of plots parts more cleanly for 3D case.
function mesh_3d_levelsetlattice(vals, isoval::Real, Rs::DirectBasis{3})
    # marching cubes algorithm to find isosurfaces (using Meshing.jl)
    algo = MarchingCubes(iso=isoval, eps=1e-3)
    verts, faces = isosurface(vals, algo; 
                              origin = SVector(-0.5,-0.5,-0.5), 
                              widths = SVector(1.0,1.0,1.0))    
    
    # transform to Cartesian basis & convert from N-vectors of 3-vectors to N×3 matrices
    verts′, faces′ = _mesh_to_cartesian(verts, faces, Rs)

    return verts′, faces′
end
function mesh_3d_levelsetlattice(flat::AbstractFourierLattice, isoval::Real, 
                                 Rs::DirectBasis{3})

    # marching cubes algorithm to find isosurfaces (using Meshing.jl)
    algo = MarchingCubes(iso=isoval, eps=1e-3)
    verts, faces = isosurface(xyz->calcfourier(xyz, flat), algo; 
                              origin = SVector(-0.5,-0.5,-0.5), 
                              widths = SVector(1.0,1.0,1.0))    
    
    # transform to Cartesian basis & convert from N-vectors of 3-vectors to N×3 matrices
    verts′, faces′ = _mesh_to_cartesian(verts, faces, Rs)
   
    return verts′, faces′
end

function _mesh_to_cartesian(verts::AbstractVector, faces::AbstractVector, Rs::DirectBasis{3})
    Nᵛᵉʳᵗˢ = length(verts); Nᶠᵃᶜᵉˢ = length(faces)
    verts′ = Matrix{Float64}(undef, Nᵛᵉʳᵗˢ, 3)
    @inbounds @simd for j in (1,2,3)  # Cartesian xyz-coordinates
        R₁ⱼ, R₂ⱼ, R₃ⱼ = Rs[1][j], Rs[2][j], Rs[3][j]
        for i in Base.OneTo(Nᵛᵉʳᵗˢ) # vertices
            verts′[i,j] = verts[i][1]*R₁ⱼ + verts[i][2]*R₂ⱼ + verts[i][3]*R₃ⱼ
        end
    end
    # convert `faces` from Nᶠᵃᶜᵉˢ-vector of 3-vectors to Nᶠᵃᶜᵉˢ×3 matrix
    faces′ = [faces[i][j] for i = Base.OneTo(Nᶠᵃᶜᵉˢ), j = Base.OneTo(3)]
    return verts′, faces′
end

# Plot from MPB-data (::AbstractFourierLattice) -------------------------------------------

function plot_lattice_from_mpbparams(filepath::String; kwargs...)
    Rs, flat, isoval, epsin, epsout = lattice_from_mpbparams(filepath)
    plot(flat, Rs; isoval=isoval, kwargs...)
    return nothing
end

# ::KVec ----------------------------------------------------------------------------------

# Plotting a single `KVec`
function plot(kv::KVec, 
              ax=plt.figure().gca(projection = dim(kv)==3 ? (using3D(); "3d") : "rectilinear"))   
    D = dim(kv)
    freeαβγ = freeparams(kv)
    nαβγ = count(freeαβγ)
    nαβγ == 3 && return ax # general point/volume (nothing to plot)

    _scatter = D == 3 ? ax.scatter3D : ax.scatter
    _plot    = D == 3 ? ax.plot3D : ax.plot
 
    if nαβγ == 0 # point
        k = kv()
        _scatter(k...)
    elseif nαβγ == 1 # line
        k⁰, k¹ = kv(zeros(D)), kv(freeαβγ.*0.5)
        ks = [[k⁰[i], k¹[i]] for i in 1:D]
        _plot(ks...)
    elseif nαβγ == 2 && D > 2 # plane
        k⁰⁰, k¹¹ = kv(zeros(D)), kv(freeαβγ.*0.5)
        αβγ, j = (zeros(3), zeros(3)), 1
        for i = 1:3
            if freeαβγ[i]
                αβγ[j][i] = 0.5
                j += 1
            end
        end
        k⁰¹, k¹⁰ = kv(αβγ[1]), kv(αβγ[2])
        # calling Poly3DCollection is not so straightforward: follow the advice
        # at https://discourse.julialang.org/t/3d-polygons-in-plots-jl/9761/3
        verts = ([tuple(k⁰⁰...); tuple(k¹⁰...); tuple(k¹¹...); tuple(k⁰¹...)],)
        plane = PyPlot.PyObject(art3D).Poly3DCollection(verts, alpha = 0.15)
        PyPlot.PyCall.pycall(plane.set_facecolor, PyPlot.PyCall.PyAny, [52, 152, 219]./255)
        PyPlot.PyCall.pycall(ax.add_collection3d, PyPlot.PyCall.PyAny, plane)
    end
    return ax
end

# Plotting a vector of `KVec`s
function plot(kvs::AbstractVector{KVec})
    D = dim(first(kvs))
    ax = plt.figure().gca(projection= D==3 ? (using3D(); "3d") : "rectilinear")
    for kv in kvs
        plot(kv, ax)
    end
    return ax
end
# Plotting a `LittleGroup`
plot(lgs::AbstractVector{<:LittleGroup}) = plot(kvec.(lgs))