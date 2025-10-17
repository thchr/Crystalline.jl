using Meshing: MarchingCubes, isosurface
using Statistics: quantile
using Crystalline: calcfouriergridded, DirectBasis, AbstractFourierLattice
using StaticArrays: SVector
import Contour as C
 # extend rather than define so it can be called externally via `Crystalline.(...)`
import Crystalline:
    _create_isosurf_plot_data,
    mesh_3d_levelsetlattice,
    isocaps_3d_levelsetlattice

# ---------------------------------------------------------------------------------------- #

function _create_isosurf_plot_data(
            flat::AbstractFourierLattice{D};
            N::Integer=(D==2 ? 100 : 30),
            filling::Union{Real, Nothing}=0.5,
            isoval::Union{Real, Nothing}=nothing) where D
 
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    if isnothing(isoval)
        isnothing(filling) && error(ArgumentError("`filling` and `isoval` cannot both be `nothing`"))
        # we don't want to "double count" the BZ edges - so to avoid that, exclude the last 
        # index of each dimension (same approach as in `MPBUtils.filling2isoval`)
        isoidxs = Base.OneTo(N-1)
        vals′ = D == 2 ? (@view vals[isoidxs, isoidxs]) :
                D == 3 ? (@view vals[isoidxs, isoidxs, isoidxs]) : error("unhandled dimension")

        isoval = quantile(Iterators.flatten(vals′), filling)
    end

    return xyz, vals, isoval
end

# Meshing utilities (for ::AbstractFourierLattice) ----------------------------------------

function mesh_3d_levelsetlattice(vals, isoval::Real, Rs::DirectBasis{3})
    # marching cubes algorithm to find isosurfaces (using Meshing.jl)
    method = MarchingCubes(iso=isoval)
    XYZ = (-0.5, 0.5)
    verts, faces = isosurface(vals, method, XYZ, XYZ, XYZ)    
    
    # transform to Cartesian basis & convert from N-vectors of 3-vectors to N×3 matrices
    verts′ = _mesh_to_cartesian(verts, Rs)
    faces′ = _face_vector_2_face_matrix(faces)

    return verts′, faces′
end
function mesh_3d_levelsetlattice(flat::AbstractFourierLattice{3}, isoval::Real, 
                                 Rs::DirectBasis{3};
                                 N = 20)
    xyz = range(-.5, .5, length=N)
    vals = calcfouriergridded(xyz, flat, N)
    return mesh_3d_levelsetlattice(vals, isoval, Rs)
end

function _mesh_to_cartesian(verts::AbstractVector, Rs::AbstractBasis{3})
    Nᵛᵉʳᵗˢ = length(verts)
    verts′ = Matrix{Float64}(undef, Nᵛᵉʳᵗˢ, 3)
    @inbounds @simd for j in (1,2,3) # Cartesian xyz-coordinates
        R₁ⱼ, R₂ⱼ, R₃ⱼ = Rs[1][j], Rs[2][j], Rs[3][j]
        for i in Base.OneTo(Nᵛᵉʳᵗˢ) # vertices
            verts′[i,j] = verts[i][1]*R₁ⱼ + verts[i][2]*R₂ⱼ + verts[i][3]*R₃ⱼ
        end
    end

    return verts′
end

function _face_vector_2_face_matrix(faces::AbstractVector)
    # convert `faces` from Nᶠᵃᶜᵉˢ-vector of 3-vectors to Nᶠᵃᶜᵉˢ×3 matrix
    Nᶠᵃᶜᵉˢ = length(faces)
    faces′ = [faces[i][j] for i in Base.OneTo(Nᶠᵃᶜᵉˢ), j in (1,2,3)]

    return faces′
end

# ---------------------------------------------------------------------------------------- #

function isocaps_3d_levelsetlattice(vals, isoval::Real, Rs::DirectBasis{3})
    N = size(vals, 1)
    xyz = range(-.5, .5, length=N)
    xyz_pad = push!(pushfirst!(collect(xyz), -0.5), 0.5)

    # pad our coordinates by 1 on each side, to allow us to artificially set the exterior to
    # something extremely large (to get contour from the boundary of the unit cell also)
    vals_slice_pad = Matrix{Float64}(undef, N+2, N+2)
    vals_slice_pad[1, :] .= 1e20
    vals_slice_pad[end, :] .= 1e20
    vals_slice_pad[:, 1] .= 1e20
    vals_slice_pad[:, end] .= 1e20

    patches_vs = Matrix{Float64}[]
    patches_fs = Matrix{Int}[]
    patches_contours = Vector{Point3f}[]
    for ss = 1:3
        vals_slice = ss == 1 ? (@view vals[1,:,:]) : 
                     ss == 2 ? (@view vals[:,1,:]) : 
                     ss == 3 ? (@view vals[:,:,1]) : error("unreachable")
        vals_slice_pad[2:end-1, 2:end-1] .= vals_slice

        cnts = C.contour(xyz_pad, xyz_pad, vals_slice_pad, isoval)
        cnts_as_points = Vector{Point2f}[]
        for line in C.lines(cnts)
            xs, ys = C.coordinates(line)
            allequal(xs) && allequal(ys) && continue # skip spurious empty lines
            push!(cnts_as_points, Point2f.(xs, ys))
        end

        # get translation for "counter-facing" face; note hack due to # hack due to ±0.5
        # difference in `convert_2d_to_3d_isocaps_reduced_coords`
        translation = (-1)^ss * Rs[ss]

        exteriors, interiors = find_exteriors_and_interiors(cnts_as_points)
        for (exterior, interior) in zip(exteriors, interiors)
            # use GeometryBasics.Polygon + faces to create a mesh from the exterior/interior
            # contours; this is why we bother to find which contours are exterior and which
            # are interior: this "splitting" is needed by `GeometryBasics.faces`
            poly = Makie.GeometryBasics.Polygon{2, Float64}(exterior, interior)
            fs = Makie.GeometryBasics.faces(poly)
            vs2D = vcat(exterior, reduce(vcat, interior; init=Vector{Point2f}[]))
            # convert 2D slice to 3D coordinates for current face (set by `ss`); still in
            # reduced coordinates
            vs3D = convert_2d_to_3d_isocaps_reduced_coords(ss, vs2D)
            vs = _mesh_to_cartesian(vs3D, Rs) # convert to Cartesian coordinates and matrix
            fs_as_matrix = _face_vector_2_face_matrix(fs)

            # now we have a mesh for this isocap face; add it
            push!(patches_vs, vs)
            push!(patches_fs, fs_as_matrix)

            # copy counter-facing sides
            push!(patches_vs, vs .+ reshape(translation, 1, 3))
            push!(patches_fs, reverse(fs_as_matrix, dims=2)) # `reverse` to flip normals
        end

        # get the boundaries of the isocaps (for plotting the "boundary line")
        for vs2D in cnts_as_points
            vs3D = convert_2d_to_3d_isocaps_reduced_coords(ss, vs2D)
            vs = _mesh_to_cartesian(vs3D, Rs) # convert to Cartesian coordinates and matrix
            vs_as_pts = Point3f.(vs[:,1], vs[:,2], vs[:,3]) # ugh; convert back again
            push!(patches_contours, vs_as_pts)
            push!(patches_contours, vs_as_pts .+ Ref(translation))
        end
    end

    return patches_vs, patches_fs, patches_contours
end

"""
    find_exteriors_and_interiors(lines)

Classify a collection of closed polygonal lines into patches. Each patch is identified with
an exterior boundary and a set of (possibly empty) interior holes.

# Strategy
The identification of exterior and interior boundaries is determined by containment
hierachy. More specifically:
1. Build a parent tree where each line's parent is the innermost line containing it
2. Compute the nesting depth of each line in this containment tree
3. Lines at even depths (0, 2, 4, ...) are patch boundaries (exteriors)
4. Lines at odd depths (1, 3, 5, ...) are holes in their parent patches
5. Each exterior's holes are its direct children in the tree

This strategy explicitly allows handling of cases where we have multiply nested lines,
corresponding to patches within holes of other patches.

# Returns
- `exteriors`: vector whose elements `exteriors[i]` give the boundary of the `i`th patch.
- `interiors`: vector whose elements `interiors[i]` give all hole boundaries of the
  `i`th patch.
"""
function find_exteriors_and_interiors(
    lines::AbstractVector{<:AbstractVector{P}}
) where P<:AbstractVector{<:Real}
    n = length(lines)
    
    # build parent tree: parent[i] stores the index j of the line that immediately 
    # contains line i (i.e., the innermost line containing i), or -1 if line i is 
    # at the root level (not contained in any other line)
    parent = fill(-1, n)
    for i in eachindex(lines)
        sample_i = lines[i][1]  # test point from line i
        for j in eachindex(lines)
            i == j && continue
            if inpolygon(sample_i, lines[j])
                # j contains i; verify that j is the *immediate* parent by ensuring
                # no intermediate line k exists where i ∈ k ∈ j. Without this check,
                # we might incorrectly assign a grandparent as the parent.
                is_immediate = true
                for k in eachindex(lines)
                    (k == i || k == j) && continue
                    sample_k = lines[k][1]
                    # check if k is between i and j in the containment hierarchy
                    if inpolygon(sample_i, lines[k]) && inpolygon(sample_k, lines[j])
                        is_immediate = false
                        break
                    end
                end
                if is_immediate
                    parent[i] = j
                    break
                end
            end
        end
    end
    
    # compute nesting depth of each line by walking up the parent chain
    depth = zeros(Int, n)
    for i in eachindex(lines)
        d = 0
        p = parent[i]
        while p ≠ -1
            d += 1
            p = parent[p]
        end
        depth[i] = d
    end
    
    # exteriors have even depth; collect each exterior with its direct children as holes
    exteriors = Vector{Vector{P}}()
    interiors = Vector{Vector{Vector{P}}}()
    for i in eachindex(lines)
        if iseven(depth[i])
            push!(exteriors, lines[i])
            # holes are the direct children of this exterior in the parent tree
            holes = [lines[j] for j in eachindex(lines) if parent[j] == i]
            push!(interiors, holes)
        end
    end

    return exteriors, interiors
end

function inpolygon(
    p::AbstractVector{T},
    polygon::AbstractVector{<:AbstractVector{T}}
) where T<:Real
    length(p) == 2 || error("`p` must be a 2D point")
    length(polygon) < 3 && error("ill-defined polygon") # must have at least 3 vertices

    # perform the Ray Casting algorithm
    px = p[1]; py = p[2]  
    is_inside = false
    N = length(polygon)
    for i in eachindex(polygon)
        p1 = polygon[i]
        p2 = polygon[mod1(i + 1, N)]
        
        y1, y2 = p1[2], p2[2]

        # check if the horizontal ray from the point intersects the edge
        # the condition `(y1 ≤ py < y2) || (y2 ≤ py < y1)` handles vertices on the ray
        if ((y1 ≤ py < y2) || (y2 ≤ py < y1))
            # calculate the x-coordinate of the intersection of the ray and the edge
            # [derived from the line equation: x = x1 + (y - y1) * (x2 - x1) / (y2 - y1)]
            # only done if an intersection in y is possible, avoiding division by zero for
            # horizontal lines
            x_intersection = (py - y1) * (p2[1] - p1[1]) / (y2 - y1) + p1[1]

            # if the intersection is to the right of the point, toggle the state
            if px < x_intersection
                is_inside = !is_inside
            end
        end
    end

    return is_inside
end

function convert_2d_to_3d_isocaps_reduced_coords(ss, vs2D)
    # create 3D vertices from those we extracted in the 2D setting (still reduced
    # coordinates); we deliberately pick -0.5 for s==2 and +0.5 otherwise, because this
    # appears to get us outward-pointing normals consistently: importantly though, we need
    # to account for this difference when we eventually add the "counter-facing sides"
    # (by incorporating a (-1)^ss factor into the translation by `Rs[ss]`)
    if ss == 1        # fixed x-side
        return [(0.5, vs2D[i][1], vs2D[i][2]) for i in eachindex(vs2D)]
    elseif ss == 2    # fixed y-side
        return [(vs2D[i][1], -0.5, vs2D[i][2]) for i in eachindex(vs2D)]
    elseif ss == 3    # fixed z-side
        return [(vs2D[i][1], vs2D[i][2], 0.5) for i in eachindex(vs2D)]
    else
        error("obtained unexpected value of `ss` not in 1:3")
    end
end
# ---------------------------------------------------------------------------------------- #
# copied from Woodpile.jl/ext/WoodpileMakieExt.jl

function facets(Rs :: DirectBasis{3}; origin::SVector{3,Float64}=SVector(0.0,0.0,0.0))
    # returns a vector of vertices (i.e., faces, ordered counter-clockwise, i.e., with 
    # normal point "outward") of the unit cell's facets
    # - `origin`: where to place the center of the unit cell
    z = origin - sum(Rs) / 2
    A1 = [z, z + Rs[3], z + Rs[2] + Rs[3], z + Rs[2]]
    A2 = [z, z + Rs[1], z + Rs[3] + Rs[1], z + Rs[3]]
    A3 = [z, z + Rs[2], z + Rs[1] + Rs[2], z + Rs[1]]
    B1 = reverse(A1) .+ (Rs[1],)
    B2 = reverse(A2) .+ (Rs[2],)
    B3 = reverse(A3) .+ (Rs[3],)
    return [A1, A2, A3, B1, B2, B3]
end

function facets(Rs :: DirectBasis{2}; origin::SVector{2,Float64}=SVector(0.0,0.0))
    z = origin - sum(Rs) / 2
    A = [z, z + Rs[1], z + Rs[1] + Rs[2], z + Rs[2], z]
    return [A]
end