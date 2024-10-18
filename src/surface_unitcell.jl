using Bravais
using StaticArrays
using LinearAlgebra: dot, cross, norm, normalize, I
using .SmithNormalForm: smith

# what are the constraints we impose on surface-lattice points {r}? they are:
#    (1)  n ⋅ r = 0 for surface lattice vectors r
#    (2)  r = ∑ⱼ aⱼRⱼ
# with integer coefficients aⱼ (j = 1,2,3); combining (1) & (2), we have
# 	 n ⋅ ∑ⱼ aⱼRⱼ = 0  ⇔  ∑ⱼ aⱼ(n ⋅ Rⱼ) = 0
function surface_basis_from_miller_indices_via_snf(
			n :: AbstractVector{<:Integer}; 
			# surface's Miller indices (normal vector ∝ n[1]*Gs[1] + n[2]*Gs[2] + n[3]*Gs[3])
			)
	length(Rs) == 3 || error("only 3D implemented; number of basis vectors is not 3")
	all(R->length(R) == 3, Rs) || error("only 3D implemented; basis vectors are not 3D")
	length(n) == 3 || error("only 3D implemented; normal vector is not 3D")

	# build a matrix to represent the relation n ⋅ r = 0; since we specify n in the a basis
	# of the reciprocal lattice vectors, and r in the direct basis, the dot product is
	# simply n ⋅ r = nᵢrⱼGᵢRⱼ = 2πnᵢrᵢ; the factors of 2π do not matter for the 0-relation,
	# so we can represent the action of n ⋅ r simply by a "matrix" nᵀ; the corresponding 
	# surface basis is simply the (integer) null space of this matrix, i.e., of nᵀ. In turn,
	# this integer null-space can be obtained as the last columns of T⁻¹ from the Smith
	# normal form of nᵀ = SΛT (the columns associated with the with zero-columns of the
	# diagonal matrix Λ)
	nᵀ = Matrix(transpose(n))
	F = smith(nᵀ)
	if F.SNF ≠ [1]
		error(lazy"failed to determine surface basis; n may have non-unity factors (n=$n)")
	end
	rs = F.Tinv[:, 2:3] # these two columns are a basis for the surface lattice
	@assert size(rs,1) == 3 && size(rs,2) == 2

	return rs
end

signed_volume(R₁, R₂, R₃) = dot(R₁, cross(R₂, R₃)) # R₁⋅(R₂×R₃) = R₃⋅(R₁×R₂) = R₂⋅(R₃×R₁)
signed_volume(r₁, r₂) = r₁[1] * r₂[2] - r₁[2] * r₂[1] # z⋅(r₁×r₂) w/ r₁,r₂ lifted to ℝ³

function rotation_matrix_n_to_z(
			n :: AbstractVector{<:Real},
			z :: AbstractVector{<:Real} = (@SVector[0.0,0.0,1.0])
			)
	# purpose: determine a rotation matrix P such that P⁻¹n = z
	# strategy: rotate around an axis perpendicular to both `n` and ̂z by an angle θ (the
	# angle that `n` makes to ̂z). The rotation matrix for a rotation around an axis `u` at
	# an angle is given by Rodrigues' rotation formula

	length(n) == length(z) == 3 || error(lazy"input must be 3D vectors; got n=$n & z=$z")

	n = normalize(SVector{3,Float64}(n)) # ̂n
	z = normalize(SVector{3,Float64}(z)) # ̂z
	I₃ = SMatrix{3,3,Float64}(I)
	n ≈ z && return I₃ # already aligned with ̂z; no rotation needed

	u = normalize(cross(n, z)) # = ̂n × ̂z (rotation axis)
	θ = acos(dot(n, z))

	# Rodrigues' rotation formula (rotate around `u` by θ): P = I + sin(θ)K + (1 - cos(θ))K²
	K = SMatrix{3,3,Float64}([0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0])
	P = I₃ + sin(θ) * K + (1 - cos(θ)) * K^2

	# computed `P` acts on `n` according to `P*n = z`; however, we want `P` to act on
	# as `P⁻¹n = z` (i.e., `P` rotates `n` to `z`); so return P⁻¹ = Pᵀ (cf. P ∈ O(3))
	return P'
end

## --------------------------------------------------------------------------------------- #

"""
	surface_basis(Rs, n; cartesian=true)

Compute a basis for the surface-cell obtained from terminating a 3D lattice `Rs` over a
surface surface specified by its normal vector `n` (or Miller indices if `cartesian=false`).

## Output

The function returns a tuple `(rs³ᴰ, rs′²ᴰ, P)`, whose elements are described below:

- `rs³ᴰ`: a `DirectBasis{3}`, whose first two basis vectors lie in the plane of the surface,
  and whose third vector is (positively) aligned with the surface normal.
  All basis vectors correspond to points in the original direct lattice.

- `rs′²ᴰ`: a `DirectBasis{2}`, whose basis vectors are given in the local coordinate system
  of the surface unit cell; effectively, this is ``(x,y)``-components of `rs³ᴰ` after a
  rotation that takes `n` to ``\\hat{\\mathbf{z}}``. The first basis vector is aligned with
  the ``\\hat{\\mathbf{x}}``-direction of the local coordinate system.

- `P`: a rotation matrix that takes `rs³ᴰ` to the local `n`-to-``\\hat{\\mathbf{z}}``
  rotated coordinate system of `rs′²ᴰ`.
  In particular, defining `rs′³ᴰ = transform.(DirectPoint.(rs³ᴰ), Ref(P))`, the following 
  holds:
  	- `rs′³ᴰ[i] ≈ [rs′²ᴰ[i]..., 0]` for `i ∈ 1:2`,
  	- `rs′³ᴰ[3]` is (positively) aligned with `[0,0,1]`.
  To transform the other way, i.e., from surface-local to lattice-global coordinates, simply
  use `P⁻¹ = transpose(P)` instead.

The returned basis is right-handed.

## Keyword arguments
- `cartesian :: Bool = true`: whether the surface normal `n` is specified in Cartesian
  coordinates (`true`) or in the basis of the reciprocal lattice vectors (`false`), i.e.,
  corresponding to the Cartesian vector `n[1]*Gs[1] + n[2]*Gs[2] + n[3]*Gs` with
  `Gs` denoting the Cartesian representation of the reciprocal lattice vectors, i.e.,
  `Gs = reciprocalbasis(Rs)`.
  The latter case (`false`) is a specification of the surface in terms of its Miller
  indices: the coordinates of `n` can then equivalently be interpreted as the inverse of the
  surface's intercept with each of the axes spanned by `Rs`.

## Example
```jldoctest
julia> sgnum = 195; # P23, a cubic lattice

julia> Rs = directbasis(sgnum, 3)
DirectBasis{3} (cubic):
 [1.0, 0.0, 0.0]
 [0.0, 1.0, 0.0]
 [0.0, 0.0, 1.0]

julia> n = [1,1,1]; # project along 111 direction

julia> rs³ᴰ, rs′²ᴰ, P = surface_basis(Rs, n; cartesian=true);

julia> rs³ᴰ # the surface basis in the original coordinate system
DirectBasis{3} (hexagonal):
 [1.0, -1.5700924586837747e-16, -0.9999999999999998]
 [-0.9999999999999998, 0.9999999999999997, 0.0]
 [1.0, 1.0, 1.0

julia> rs′²ᴰ # the in-plane surface vectors in a local "surface coordinate system"
DirectBasis{2} (hexagonal):
 [1.414213562373095, 0.0]
 [-0.7071067811865475, 1.2247448713915887]

julia> DirectBasis(transform.(DirectPoint.(rs³ᴰ), Ref(P))) # from rs³ᴰ to rs′²ᴰ coordinates 
DirectBasis{3} (hexagonal):
 [1.414213562373095, -1.1102230246251563e-16, 0.0]
 [-0.7071067811865476, 1.2247448713915885, 0.0]
 [0.0, -1.1102230246251563e-16, 1.7320508075688772]

!!! warning
   This function will likely be moved from Crystalline to Bravais.jl at some point in the
   future.
```
"""
function surface_basis(
			Rs :: DirectBasis{3}, 
			n :: AbstractVector{<:Real}; 
			cartesian :: Bool = true
		)
	# TODO: Add tests
	# TODO: Move this from Crystalline to Bravais (but this is annoying, because we need
	#       SmithNormalForm.jl, which is not a dependency of Bravais - and not even
	#	    registered, but vendored by Crystalline). Messy :(

	length(n) == 3 || throw(DimensionMismatch(lazy"normal vector must be 3D, got $n"))

	# if `cartesian=true` (default), we assume `n` has been given in Cartesian coordinates;
	# otherwise, if `cartesian=false`, `n` is assumed to be given in the basis of 
	# `reciprocalbasis(Rs)`, i.e., to be a set of Miller indices for the surface
	local n_cartesian
	if cartesian
		n_cartesian = n
		n = latticize(n, reciprocalbasis(Rs))
	else
		n_cartesian = cartesianize(float.(n), reciprocalbasis(Rs))
	end

	# cast `n` as the smallest integer vector proportional to `n`, so we can do integer math
	# from here on out
	local n_int
	if !(eltype(n) <: Integer)
		d = float_gcd(n...; atol=1e-9)
		_n = n ./ d
		n_int = round.(Int, _n)
		if norm(_n - n_int) > 1e-9
			error(LazyString(
					"failed to convert `n` to an integer vector with precision 1e-9; ",
					" got n=", _n, " vs. n_int=", n_int, 
					"; if possible, try to give `n` as integer Miller indices"))
		end
	else
		d = gcd(n)
		n_int = isone(d) ? convert.(Int, n) : convert.(Int, div.(n, d))
	end

    # find a basis for the surface lattice, using the constraint n ⋅ rⱼ = 0 for each surface
	# basis vector rⱼ; the below returns two basis vectors, in matrix form, specified in the
	# basis of the direct lattice vectors of the bulk
    rs_snf = surface_basis_from_miller_indices_via_snf(n_int)

	_rs = stack(Rs) * rs_snf # equiv. to `cartesianize(r₃_snf, Rs)`, but need different type
	r₁ = DirectPoint{3}(_rs[1,1], _rs[2,1], _rs[3,1])
	r₂ = DirectPoint{3}(_rs[1,2], _rs[2,2], _rs[3,2])

	# rotate the surface unit cell to have a normal along z, so we can embed it in 2D by
    # ignoring the third dimension
    P_to_z = rotation_matrix_n_to_z(n_cartesian)
	r₁′ = transform(r₁, P_to_z)
	r₂′ = transform(r₂, P_to_z)
	atol = 1e-9 * maximum(norm, Rs)
	if abs(r₁′[3]) > atol && abs(r₂′[3]) > atol
		error(lazy"rotation to 2D plane failed; projection to z is appreciably nonzero")
	end

	# if the basis is already right-handed (with out-of-plane direction =z), we can proceed
	# directly but otherwise we swap r₁ and r₂ before proceeding (specifically, we require
	# that a positive signed volume, i.e. ̃n ⋅ (r₁ × r₂) > 0, for the final output)
	if signed_volume(r₁′, r₂′) < 0
		r₁′, r₂′ = r₂′, r₁′
	end
	tmp_rs′_2D = DirectBasis{2}(r₁′[SOneTo(2)], r₂′[SOneTo(2)])

	# compute the associated 2D Niggli cell of the rotated basis vectors
	rs′²ᴰ, _ = nigglibasis(tmp_rs′_2D)
	r₁′²ᴰ = DirectPoint{2}(rs′²ᴰ[1])
	r₂′²ᴰ = DirectPoint{2}(rs′²ᴰ[2])

	# in case we to use the 2D projection later (i.e., `rs′²ᴰ`), it would be nice if 
	# its first vector was aligned with its local x-axis; one more rotation to ensure this
	θ = atan(r₁′²ᴰ[2], r₁′²ᴰ[1])
	P_to_x = @SMatrix [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
	P_to_x_2D = P_to_x[SOneTo(2), SOneTo(2)]
	r₁′²ᴰ, r₂′²ᴰ = transform.((r₁′²ᴰ, r₂′²ᴰ), Ref(P_to_x_2D))
	rs′²ᴰ = DirectBasis{2}(r₁′²ᴰ, r₂′²ᴰ)

	# overall rotation matrix going back from 2D to 3D
	P = P_to_z * P_to_x # correct order since P acts inversely
    P⁻¹ = transpose(P) # P⁻¹ = Pᵀ

	# rotate back from the 2D projection to 3D
	r₁′³ᴰ = DirectPoint{3}(r₁′²ᴰ..., 0.0)
	r₂′³ᴰ = DirectPoint{3}(r₂′²ᴰ..., 0.0)
	r₁³ᴰ = transform(r₁′³ᴰ, P⁻¹)
	r₂³ᴰ = transform(r₂′³ᴰ, P⁻¹)

	# now: find the "last" or "3rd" basis vector r₃, which we want to be aligned with `n`,
	# i.e., is normal to the surface; this requires that n×r₃ = 0; to ensure this, we must
	# construct a matrix that represents the action of n× on r₃ (accounting for the fact 
	# that n specified in the basis of the reciprocal lattice vectors, and r₃ is specified
	# in the direct basis); we denote this matrix C, s.t. C*r₃ = n×r₃ = 0; we seek a basis
	# for the null-space of C and we will use the smith normal form to get it again
	_Gs = reciprocalbasis(Rs) ./ (2π) # factors of 2π are irrelevant & troublesome; remove
	C = reduce(hcat, (sum(cross(n_int[k]*_Gs[k], Rs[j]) for k in 1:3) for j in 1:3)) # (n×r₃)ᵢ = Cᵢⱼ(r₃)ⱼ
	d = float_gcd(C...; atol=1e-9) # convert to nearest integer matrix
	norm(d) ≈ 1 || (C = C ./ d)
	C_int = round.(Int, C)
	isapprox(C_int, C, atol=1e-8) || error("failed to convert to integer matrix")
	F = smith(Matrix(C_int))
	if !(length(F.SNF) == 3 && F.SNF[1] ≠ 0 && F.SNF[2] ≠ 0 && F.SNF[3] == 0)
		error("unexpected Smith normal form for 3rd surface-cell basis vector calculation")
	end
	r₃_snf = DirectPoint{3}(F.Tinv[:,3]) # 3rd basis vector for surface-cell, in `Rs`-basis
	r₃³ᴰ = cartesianize(r₃_snf, Rs)
	if dot(r₃³ᴰ, n_cartesian) < 0
		r₃³ᴰ = -r₃³ᴰ # turn "negative" parallel to "positive" parallel
	end

	# check that everything is still good: r₃ is aligned with `n_cartesian` and the basis is
	# right-handed
	rs³ᴰ = DirectBasis{3}(r₁³ᴰ, r₂³ᴰ, r₃³ᴰ)
	if (signed_volume(r₁′²ᴰ, r₂′²ᴰ) ≤ 0 || 
	    signed_volume(r₁³ᴰ, r₂³ᴰ, r₃³ᴰ) ≤ 0 ||
	    norm(cross(r₃³ᴰ, n_cartesian)) > 1e-9 ||
		dot(r₃³ᴰ, n_cartesian) < 0)
	   	error("something went wrong during the surface basis computation; signed volume \
	   		   is negative or the out-of-plane (third) surface-unit-cell basis vector is \
			   not (positively) oriented towards the surface normal vector")
	end

    return rs³ᴰ, rs′²ᴰ, P
end

float_gcd(x,y; atol=1e-9) = abs(y) < atol ? x : float_gcd(y, x % y; atol)
float_gcd(x, y, z...; kws...) = float_gcd(x, float_gcd(y, z...);  kws...)

##
# Example for visualizing the surface Brillouin zone, and a surface BZ path, together
# with the bulk BZ and the bulk path

# TODO: Add this either to a tutorial/example/docstring thing, or port parts of it to
#       Brillouin.jl
#=
function plot_bz_projection(sgnum, Rs, rs³ᴰ, rs′²ᴰ; axis=NamedTuple())
	Gs = reciprocalbasis(Rs)
	kp = irrfbz_path(sgnum, conventionalize(Rs, centering(sgnum)))
	
	gs³ᴰ = reciprocalbasis(rs³ᴰ)
	# pick least-symmetric space group number in 2D
	csys = crystalsystem(rs′²ᴰ)
	sgnum²ᴰ = csys == "hexagonal"   ? 13 #= p3 =#   :
			  csys == "square"      ? 10 #= p4 =#   :
			  csys == "rectangular" ? 3  #= p1m1 =# :
			  csys == "oblique"     ? 1  #= p1 =#   : error()
	
	kp′²ᴰ = irrfbz_path(sgnum²ᴰ, rs′²ᴰ)
	gs′²ᴰ = basis(kp′²ᴰ)
	c′²ᴰ = wignerseitz(gs′²ᴰ)

	vol_Gs = Crystalline.volume(Gs)
	vol_gs²ᴰ = Crystalline.volume(gs′²ᴰ)
	α = vol_Gs/vol_gs²ᴰ / norm(gs³ᴰ[3])

	kp³ᴰ = KPath{3}(
		Dict(klab => SVector{3,Float64}(kv..., 0) for (klab, kv) in kp′²ᴰ.points), # points
		kp′²ᴰ.paths, # paths
		gs³ᴰ, # basis
		Ref(Brillouin.LATTICE) # setting
	)
	kp³ᴰ_shift = KPath{3}(
		Dict(klab => SVector{3,Float64}(kv..., α/2) for (klab, kv) in kp′²ᴰ.points), # points
		kp′²ᴰ.paths, # paths
		gs³ᴰ, # basis
		Ref(Brillouin.LATTICE) # setting
	)
	#kp³ᴰ = latticize(cartesianize(kp³ᴰ), Gs)
	#kp³ᴰ_shift = latticize(cartesianize(kp³ᴰ_shift), Gs)
	c = wignerseitz(Gs)

	faxp = plot(c; axis)
	plot!(kp; linecolor=:red, markercolor=:red)
	proj_col = RGBf(.75, .75, .75)
	for f = (-1, 1)
		_c = Cell{3}([SVector(v²ᴰ..., f*α/2) for v²ᴰ in c′²ᴰ.verts],
					 c′²ᴰ.faces, gs³ᴰ, Ref(Brillouin.LATTICE))
		plot!(_c; color= f == 1 ? :black : proj_col)
	end
	plot!(kp³ᴰ_shift; linecolor=:green, markercolor=:green)
	for v²ᴰ in c′²ᴰ.verts
		lines!(cartesianize.([SVector(v²ᴰ..., -α/2), SVector(v²ᴰ..., α/2)], Ref(gs³ᴰ));
			color=proj_col)
	end

	return faxp, kp, kp³ᴰ, α
end

##
using Crystalline
using Brillouin
#using CairoMakie
#CairoMakie.activate!(type = "svg")
using GLMakie
GLMakie.activate!()

sgnum = 216
Rs = primitivize(directbasis(sgnum), centering(sgnum))
n = [1,1,1]
rs³ᴰ, rs′²ᴰ, P = surface_basis(Rs, n; cartesian=true)
rs³ᴰ

faxp, kp, kp³ᴰ, α=plot_bz_projection(sgnum, Rs, rs³ᴰ, rs′²ᴰ; 
		axis=(;elevation=0.4, azimuth=-4.94, perspectiveness=0.15))
faxp
#save("bz_projection_$(join(n)).pdf", faxp)
=#