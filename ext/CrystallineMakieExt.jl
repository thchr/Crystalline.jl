module CrystallineMakieExt
using Crystalline
using StaticArrays
using Makie
using Crystalline: AbstractFourierLattice

# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# ::AbstractFourierLattice, ::DirectBasis

# ---------------------------------------------------------------------------------------- #

# ↓ implementations of `_create_isosurf_plot_data`, `mesh_3d_levelsetlattice`, 
# `isocaps_3d_levelsetlattice`, & `facets(::AbstractBasis)`
include("levelsetplot_meshing.jl")

# ---------------------------------------------------------------------------------------- #

"""
    levelsetplot(flat::AbstractFourierLattice, Rs::DirectBasis; kws...)

Plot a level-set lattice `flat::AbstractFourierLattice` with lattice vectors given by 
`Rs::DirectBasis`. 

## Keyword arguments

The following keyword arguments affect the isosurface generation and the structure of the
level-set lattice:

- `samples`: number of samples used to build isosurface (default: `100` in 2D, `20` in 3D).
- `filling`: determine isovalue from relative filling fraction (default: `0.5`)
- `isoval`: isovalue used for isosurface (default: `nothing`); overrides `filling` if
  provided.

In addition, several keyword arguments control the appearance of the plot, including
`color`, `alpha`, `isocap_alpha`, `strokecolor`, `strokewidth`, `unitcell_linecolor`, &
`unitcell_linewidth`.

## Examples

Compute a modulated `FourierLattice` for a lattice in plane group 16 and plot it via
GLMakie.jl:

```julia-repl
julia> using Crystalline, GLMakie
julia> flat = modulate(levelsetlattice(16, Val(2)))
julia> Rs = directbasis(16, Val(2)) 
julia> plot(flat, Rs)
```
"""
@recipe LevelSetPlot (flat, Rs) begin
    "number of samples used to create the isosurface (default, `nothing`; i.e., automatically set based on dimensionality (30 for 3D, 100 for 2D))"
    samples = nothing
    "determine isovalue from relative filling fraction (default, `0.5`); overruled if `isoval` is set"
    filling = 0.5
    "level-set surface isovalue (default, `nothing`, i.e., inferred from `filling`; overrules `filling` if set)"
    isoval = nothing
    "mesh patch color (default, `:grey80`)"
    color = :grey80
    "mesh patch opacity (default, `1.0`)"
    alpha = 1.0
    "opacity of isocap faces (default, `1.0`); only relevant in 3D"
    isocap_alpha = 1.0
    "color of boundary lines (default, `:grey40`)"
    strokecolor = :grey40
    "line width of boundary lines (default, `2`)"
    strokewidth = 2
    "color of unit cell boundary lines (default, `:black`)"
    unitcell_linecolor = :black
    "line width of unit cell boundary (default, `2`)"
    unitcell_linewidth = 2
    visible = true # Added manually cf. https://github.com/MakieOrg/Makie.jl/pull/5248 (TODO: remove when merged & released)
    # TODO: add `repeat`?
end

function Makie.plot!(p::LevelSetPlot{<:Tuple{<:AbstractFourierLattice{3}, DirectBasis{3}}})
    # initialize parameters
    flat = p.flat[]::AbstractFourierLattice{3}
    Rs = p.Rs[]::DirectBasis{3}
    samples = (isnothing(p.samples[]) ? 30 : p.samples[])::Int
      
    # isosurface
    _, vals, isoval = _create_isosurf_plot_data(
        flat; N = samples, filling = p.filling[], isoval = p.isoval[])

    verts′, faces′ = mesh_3d_levelsetlattice(vals, isoval, Rs)
    Makie.mesh!(p, verts′, faces′; color = p.color, alpha = p.alpha)

    # isocaps: mesh
    patches_vs, patches_fs, patches_contours = isocaps_3d_levelsetlattice(vals, isoval, Rs)
    for (vs, fs) in zip(patches_vs, patches_fs)
        Makie.mesh!(p, vs, fs; color = p.color[], alpha=p.isocap_alpha[])
    end

    # isocaps: contours / boundary lines
    for vs in patches_contours
        Makie.lines!(p, vs; color = p.strokecolor, linewidth=p.strokewidth)
    end
    
    # unit cell
    uc_facets = facets(Rs)
    for facet in uc_facets
        Makie.lines!(p, facet; color=p.unitcell_linecolor, linewidth=p.unitcell_linewidth)
    end

    return p
end

function Makie.plot!(p::LevelSetPlot{<:Tuple{<:AbstractFourierLattice{2}, DirectBasis{2}}})
    # initialize parameters
    flat = p.flat[]::AbstractFourierLattice{2}
    Rs = p.Rs[]::DirectBasis{2}
    samples = (isnothing(p.samples[]) ? 100 : p.samples[])::Int
    color = p.color[]
    color_rgba = parse(Makie.RGBA, color)

    # isosurface data
    xy, vals, isoval = _create_isosurf_plot_data(
        flat; N = samples, filling = p.filling[], isoval = p.isoval[])

    # filled isocontour
    X = [x*Rs[1][1] + y*Rs[2][1] for x in xy, y in xy] # in Cartesian coords, in format
    Y = [x*Rs[1][2] + y*Rs[2][2] for x in xy, y in xy] # needed by `contourf`
    Makie.contourf!(
        p, X, Y, vals;
        levels = [isoval], colormap = [color_rgba],
        extendlow=:auto, # do a single-contour, and fill in all values below it
    )
    #= model = Makie.Mat4f( # TODO: implement via `model` transform rather than `X` & `Y`.
        Rs[1][1], Rs[1][2], 0.0, 0.0,
        Rs[2][1], Rs[2][2], 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    ) =#
    
    # isocontour lines (cannot currently do this in `contourf`, so at the moment, we just
    # rerun the contour calculation entirely - ugh: TODO: fix if/when Makie supports this)
    Makie.contour!(
        p, X, Y, vals; 
        levels = [isoval], color = p.strokecolor, linewidth = p.strokewidth
    )

    # unit cell
    uc_facets = facets(Rs)
    for facet in uc_facets
        Makie.lines!(p, facet; color=p.unitcell_linecolor, linewidth=p.unitcell_linewidth)
    end

    return p
end


# alias `levelsetplot(!)` to `plot(!)`
Makie.plottype(::AbstractFourierLattice{D}, ::DirectBasis{D}) where D = LevelSetPlot

function Makie.args_preferred_axis(
    ::Type{<:LevelSetPlot},
    ::AbstractFourierLattice{D},
    ::DirectBasis{D},
) where D
    return D == 3 ? Axis3 : Axis
end

# ---------------------------------------------------------------------------------------- #

function Makie.plot(
    flat::Union{Observable{<:L}, L} where L <: AbstractFourierLattice{D},
    Rs::Union{Observable{DirectBasis{D}}, DirectBasis{D}};
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis...)
    p = levelsetplot!(ax, flat, Rs; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

# ---------------------------------------------------------------------------------------- #
# default axis settings in 3D and 1D+2D

function _default_bare_axis!(f, ::Val{3}; axis_kws...)
    f[1,1] = ax = Makie.Axis3(f;
        aspect=:data,
        viewmode=:free,
        axis_kws...)
    Makie.hidedecorations!(ax)
    ax.protrusions[] = 0 # cf. https://github.com/MakieOrg/Makie.jl/issues/2259
    Makie.hidespines!(ax)

    return ax
end

function _default_bare_axis!(f, ::Union{Val{1}, Val{2}}; axis_kws...)
    f[1,1] = ax = Makie.Axis(f;
        aspect=Makie.DataAspect(),
        axis_kws...)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end

# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# ::DirectBasis

"""
    basisplot(Vs::AbstractBasis; kws...)

Plot a `Vs::AbstractBasis` and its parallepiped unit cell.

## Keyword arguments

The following keyword arguments affect the isosurface generation and the structure of the
level-set lattice:
 
- `origin :: AbstractVector`: a shift of the lattice origin (default: zero-vector)

## Examples

```julia-repl
julia> using Crystalline, GLMakie
julia> Rs = directbasis(16, Val(2))
julia> plot(Rs)
```
"""
@recipe BasisPlot (Vs,) begin
    "shift of the lattice origin (default: zero-vector)"
    origin = nothing
    visible = true # Added manually cf. https://github.com/MakieOrg/Makie.jl/pull/5248 (TODO: remove when merged & released)
end

function Makie.plot!(
    p::BasisPlot{<:Tuple{<:AbstractBasis{D}}}
) where D
    Vs = p.Vs[]::AbstractBasis{D}
    origin = p.origin[]
    if !isnothing(p.origin[])
        origin isa AbstractVector || error("`origin` must be an AbstractVector")
        length(origin) == D || error("`origin` must have length $D")
        Vs = Vs .+ Ref(SVector{D}(origin))
        origin = SVector{D, eltype(eltype(Vs))}(origin)
    else
        origin = zeros(SVector{D, eltype(eltype(Vs))})
    end :: SVector{D, eltype(eltype(Vs))}

    if D ≠ 1
        uc_lower_corner = sum(Vs) ./ 2
        for facet in facets(Vs; origin=uc_lower_corner)
            Makie.lines!(p, facet; color="black", linewidth = 2)
        end
    end

    for V in Vs
        start = D == 1 ? Point2f(origin[1], 0) : Point(origin)
        destination = D == 1 ? Point2f(V[1], 0) : Point(V)
        Makie.arrows2d!(p, start, destination; color=:blue, shaftwidth = 3)
    end
    
    return p
end

Makie.plottype(::AbstractBasis{D}) where D = BasisPlot # alias `basisplot(!)` to `plot(!)`

function Makie.args_preferred_axis(
    ::Type{<:BasisPlot},
    ::AbstractBasis{D},
) where D
    return D == 3 ? Axis3 : Axis
end

function Makie.plot(
    Vs::Union{Observable{AbstractBasis{D}}, AbstractBasis{D}};
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    f = Makie.Figure(; figure...)
    ax = _default_bare_axis!(f, Val(D); axis...)
    p = basisplot!(ax, Vs; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

# ---------------------------------------------------------------------------------------- #
end # module CrystallineMakieExt