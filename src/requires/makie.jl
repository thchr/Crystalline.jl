using Makie
using Meshing
using Statistics: quantile

# Defaults --------------------------------------------------------------------------------


# ::AbstractBasis ---------------------------------------------------------------------------

# @recipe(BasisPlot) do scene
#     Attributes(
#         marker = :circle,
#         markercolor= :white,
#         markeredgecolor = :black,
#         markerlinewidth = 1.5,
#         markersize = 4.5,
#     )
# end

function Makie.plot!(
            ax::Block,
            Vs::Union{Observable{<:AbstractBasis{D}}, <:AbstractBasis{D}};
            cntr::AbstractVector{<:Real} = zeros(SVector{D, Float64})
            kws...) where D

    Vs′ = Vs .+ Ref(cntr) # basis vectors translated by cntr
    
    if D == 1
        lines!(ax, @SVector [cntr[1], Vs′[1]], [0, 0])
        scatter!(ax, @SVector cntr[1], 0; #= TODO: marker opts =#) # origin

    elseif D == 2
        corner = sum(Vs) + cntr
        for V′ in Vs′
            lines!(ax, [cntr[1], V′[1]], [cntr[2], V′[2]]; color=:black) # basis vectors
            lines!(ax, [V′[1], corner[1]], [V′[2], corner[2]]; color=:grey) # remaining unit cell boundaries
        end
        scatter!(ax, cntr[1], cntr[2]; #= TODO: marker opts =#) # origin
    elseif D == 3
        corners = (Vs[1]+Vs[3], Vs[1]+Vs[2], Vs[2]+Vs[3]) .+ Ref(cntr)
        dirs = ((-1,1,-1), (-1,-1,1), (1,-1,-1))
        for V′ in Vs′
            lines!(ax, [cntr[1], V′[1]], [cntr[2], V′[2]], [cntr[3], V′[3]]; color=:black) # basis vectors
        end
        for (i,V) in enumerate(Vs)
            for (corner,dir) in zip(corners,dirs) # remaining unit cell boundaries
                lines!(ax, [corner[1], corner[1]+dir[i]*V[1]], 
                           [corner[2], corner[2]+dir[i]*V[2]], 
                           [corner[3], corner[3]+dir[i]*V[3]]; color=:grey)
            end
        end
        scatter!(ax, cntr[1], cntr[2], cntr[3]; #= TODO: marker opts =#) # origin
    end
    return ax
end

function Makie.plot!(
            Vs::Union{Observable{<:AbstractBasis{D}}, <:AbstractBasis{D}};
            cntr::AbstractVector{<:Real} = zeros(SVector{D, Float64})
            kws...) where D
    Makie.plot!(Makie.current_axis(), Vs; cntr, kws...)
end

# ::AbstractFourierLattice ----------------------------------------------------------------

"""
    plot(flat::AbstractFourierLattice, Rs::DirectBasis; kws...)

Plots a lattice `flat::AbstractFourierLattice` with lattice vectors
given by `Rs::DirectBasis`. 

## Keywords `kws`

- `N`: resolution (default:, `100` in 2D, `20` in 3D)
- `filling`: determine isovalue from relative filling fraction (default: `0.5`)
- `isovalue`: isovalue (default: `nothing`, i.e. inferred from `filling`)
- `repeat`: if not `nothing`, repeats the unit cell an integer number of times
  (default: `nothing`)

If both `filling` and `isoval` kwargs are simultaneously not equal to `nothing`, then
`isovalue` takes precedence.
"""
function Makie.plot!(
            ax::Makie.Block,
            flat::AbstractFourierLattice{D},
            Rs::DirectBasis{D};
            N::Integer = (D==2 ? 100 : 20),
            filling::Union{Real, Nothing} = 0.5,
            isovalue::Union{Real, Nothing} = nothing,
            repeat::Union{Integer, Nothing} = nothing) where D
 
    xyz, vals, isovalue = _create_isosurf_plot_data(flat; N, filling, isovalue)
    plotiso!(ax, xyz, vals, isovalue, Rs, repeat, fig)

    return xyz, vals, isovalue
end

function _create_isosurf_plot_data(
            flat::AbstractFourierLattice;
            N::Integer=(D==2 ? 100 : 20),
            filling::Union{Real, Nothing}=0.5,
            isoval::Union{Real, Nothing}=nothing) where D
 
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    if isnothing(isoval)
        isnothing(filling) && error(ArgumentError("`filling` and `isoval` cannot both be `nothing`"))
        # we don't want to "double count" the BZ edges - so to avoid that, exclude the last 
        # index of each dimension (same approach as in `MPBUtils.filling2isoval`)
        isoidxs = OneTo(N-1)
        vals′ = if D == 2;     (@view vals[isoidxs, isoidxs])
                elseif D == 3; (@view vals[isoidxs, isoidxs, isoidxs])
                end
        isoval = quantile(Iterators.flatten(vals′), filling)
    end

    return xyz, vals, isoval
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
    ax = fig.add_subplot(projection= D==3 ? (using3D(); "3d") : "rectilinear")

    if D == 2
    # convert to a cartesian coordinate system rather than direct basis of Ri
    N = length(xyz) 
    X = broadcast((x,y) -> x*Rs[1][1] + y*Rs[2][1], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
    Y = broadcast((x,y) -> x*Rs[1][2] + y*Rs[2][2], reshape(xyz,(1,N)), reshape(xyz, (N,1)))
    # note: the x-vs-y ordering has to be funky this way, because plotting routines expect 
    #       x to vary across columns and y to vary across rows - sad :(. See also the
    #       `calcfouriergridded` method.

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
    #caps_verts′, caps_faces′ = mesh_3d_levelsetisocaps(vals, isoval, Rs)

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
    verts, faces = isosurface(flat, algo; 
                              origin = SVector(-0.5,-0.5,-0.5), 
                              widths = SVector(1.0,1.0,1.0))    
    
    # transform to Cartesian basis & convert from N-vectors of 3-vectors to N×3 matrices
    verts′, faces′ = _mesh_to_cartesian(verts, faces, Rs)
   
    return verts′, faces′
end

function _mesh_to_cartesian(verts::AbstractVector, faces::AbstractVector, Rs::AbstractBasis{3})
    Nᵛᵉʳᵗˢ = length(verts); Nᶠᵃᶜᵉˢ = length(faces)
    verts′ = Matrix{Float64}(undef, Nᵛᵉʳᵗˢ, 3)
    @inbounds @simd for j in (1,2,3) # Cartesian xyz-coordinates
        R₁ⱼ, R₂ⱼ, R₃ⱼ = Rs[1][j], Rs[2][j], Rs[3][j]
        for i in OneTo(Nᵛᵉʳᵗˢ) # vertices
            verts′[i,j] = verts[i][1]*R₁ⱼ + verts[i][2]*R₂ⱼ + verts[i][3]*R₃ⱼ
        end
    end
    # convert `faces` from Nᶠᵃᶜᵉˢ-vector of 3-vectors to Nᶠᵃᶜᵉˢ×3 matrix
    faces′ = [faces[i][j] for i in OneTo(Nᶠᵃᶜᵉˢ), j in (1,2,3)]
    return verts′, faces′
end

# ::KVec ----------------------------------------------------------------------------------

# Plotting a single `KVec` or `RVec`
function plot(kv::AbstractVec{D}, 
              ax=plt.axes(projection = D==3 ? (using3D(); "3d") : "rectilinear")) where D

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

# Plotting a vector of `KVec`s or `RVec`s
function plot(vs::AbstractVector{<:AbstractVec{D}}) where D
    ax = plt.axes(projection= D==3 ? (using3D(); "3d") : "rectilinear")
    for v in vs
        plot(v, ax)
    end
    return ax
end
# Plotting a `LittleGroup`
plot(lgs::AbstractVector{<:LittleGroup}) = plot(position.(lgs))