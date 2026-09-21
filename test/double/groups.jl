using Crystalline, Test, LinearAlgebra, StaticArrays
using Crystalline: SU2, SU2_BY_ROTATION, SU2_BY_ROTATION_HEX, SU2_BINARY_AXES, su2

@testset "Double groups" begin

@testset "SU2 algebra" begin
    u = SU2(cis(π/3), 0)
    @test u * inv(u) == one(SU2)

    # an `SU2` is a 2×2 matrix; powers go through `*`, and stay `SU2`s
    @test u isa AbstractMatrix{ComplexF64} && size(u) == (2, 2)
    @test Matrix(u) == matrix(u) && u[1, 2] == u.b && u[2, 1] == -conj(u.b)
    @test u^2 isa SU2 && u^2 == u*u && u^3 ≈ u*u*u
    # `u` is a three-fold rotation: `u³` is the rotation by 2π, i.e. Ē
    @test u^3 ≈ -one(SU2) && u^6 ≈ one(SU2)
    @test u^-1 ≈ inv(u) && isone(u^0)
    @test det(u) ≈ 1 && tr(u) ≈ 2real(u.a)
    @test isapprox(SU2(matrix(u)), u)
    @test_throws DomainError SU2(0.5, 0.0)                    # not normalized, det ≠ 1
    @test_throws DomainError SU2(ComplexF64[1 0; 0 2])       # not of the form [a b; -b* a*]

    # products stay normalized: a long chain must not drift out of tolerance
    v = foldl(*, fill(u, 10_000))
    @test abs(abs2(v.a) + abs2(v.b) - 1) < Crystalline.DEFAULT_ATOL
end

@testset "Altmann's conventions" begin
    # Altmann & Herzig (1994), eq. 16: the binary rotations about x, y, z. Reproducing these
    # is what pins our convention to Bilbao's; `R̃(πz) = diag(-i, i)` is their `2₀₀₁` cell.
    Rs = Crystalline.directbasis(221, Val(3))     # cubic: a₁∥x, a₂∥y, a₃∥z
    @test matrix(su2(S"x,-y,-z", 221)) ≈ ComplexF64[0 -im; -im 0]     # R̃(πx)
    @test matrix(su2(S"-x,y,-z", 221)) ≈ ComplexF64[0  -1;   1 0]     # R̃(πy)
    @test matrix(su2(S"-x,-y,z", 221)) ≈ ComplexF64[-im 0;  0 im]     # R̃(πz)

    # spin is axial, so inversion maps to the identity (Altmann eq. 11; Elcoro line 121)
    @test su2(S"-x,-y,-z", 221) == one(SU2)

    # No tabulated element is the negative of another, and none has `real(a) < 0`. Together
    # these make the barring decidable from the SU(2) element alone, and hence `compose`
    # independent of the setting.
    S = collect(values(SU2_BY_ROTATION)) ∪ collect(values(SU2_BY_ROTATION_HEX))
    @test !any(isapprox(u, -v) for u in S, v in S)
    @test all(u -> real(u.a) > -Crystalline.DEFAULT_ATOL, S)
end

@testset "Tabulated SU(2) elements match Altmann eq. 8 where geometry fixes them" begin
    # For every rotation with θ ≠ π the SU(2) element is fixed outright by
    # `U = cos(θ/2)𝟙 - i sin(θ/2)(𝐧⋅𝛔)`, so the table must reproduce it. Two-fold rotations
    # and mirrors have θ = π, where ±U describe the same spatial operation and the choice is
    # convention; those are covered by the axis tests instead.
    function altmann_su2(op, Rs)
        W  = rotation(op)
        W′ = det(W) < 0 ? -W : W            # spin is axial: use the proper part
        A  = hcat(Rs...)
        R  = A * W′ * inv(A)                      # the proper rotation, in Cartesian
        c  = round(2*clamp((tr(R) - 1)/2, -1, 1))/2       # cos θ ∈ {1, ½, 0, -½, -1}
        ax = SVector{3,Float64}(R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])/2  # sin(θ)·𝐧
        norm(ax) < 0.1 && return nothing    # θ = 0 or π: not fixed by geometry
        n = ax/norm(ax)
        c₂, s₂ = sqrt((1+c)/2), sqrt((1-c)/2)
        return SU2(c₂ - im*s₂*n[3], -(n[2] + im*n[1])*s₂)
    end
    # Bilbao's Cartesian frame for hexagonal/trigonal settings is the mirror image of
    # `directbasis`'s through the x-axis; everywhere else the two agree.
    function bilbao_basis(sgnum)
        Rs = Crystalline.directbasis(sgnum, Val(3))
        crystalsystem(sgnum, 3) ∈ ("hexagonal", "trigonal") || return Rs
        return [Rs[1], SVector(Rs[2][1], -Rs[2][2], Rs[2][3]), Rs[3]]
    end

    ntested = 0
    for sgnum in 1:230
        Rs = bilbao_basis(sgnum)
        for op in spacegroup(sgnum, Val(3))
            u = altmann_su2(op, Rs)
            u === nothing && continue
            @test isapprox(su2(op, sgnum), u; atol=1e-8)
            ntested += 1
        end
    end
    # of the 4425 operations across the 230 space groups, 1888 have θ ≠ 0, π and are thus
    # fixed by geometry; the guard is against the loop silently testing nothing
    @test ntested == 1888
end

@testset "Double group structure" begin
    for sgnum in (1, 2, 75, 143, 186, 221, 225)
        sg, dsg = spacegroup(sgnum, Val(3)), spacegroup(sgnum, Val(3), Val(true))
        @test dsg isa DSpaceGroup{3}
        @test length(dsg) == 2*length(sg)             # a double group holds both cosets
        @test count(isbarred, dsg) == length(sg)
        @test all(d -> findfirst(q -> isapprox(q, d), dsg) !== nothing,
                  (a*b for a in dsg, b in dsg))
        # Ē: the 2π rotation, present, barred, of order 2
        Ē = findfirst(d -> isbarred(d) && isone(SymOperation(d)), dsg)
        @test Ē !== nothing && isone(dsg[Ē]^2) && !isone(dsg[Ē])
        @test all(d -> isone(d * inv(d)) && isone(inv(d) * d), dsg)
    end
    @test @inferred(spacegroup(221, Val(3), Val(true))) isa DSpaceGroup{3}
    # the spinless path must be untouched by the `Val(false)` route
    @test spacegroup(221, Val(3), Val(false)) == spacegroup(221, Val(3))
    @test spacegroup(221, 3, true) == spacegroup(221, Val(3), Val(true))  # unstable form
end

@testset "Barring" begin
    dsg = spacegroup(186, Val(3), Val(true))
    i = findfirst(d -> seitz(d) == "3₀₀₁⁺", dsg)
    @test i !== nothing
    @test seitz(dsg[i]^3) == "ᵈ1"          # a 2π rotation is Ē, not the identity
    @test !isone(dsg[i]^3) && isone(dsg[i]^6)

    # The case that no rule but Altmann's pole convention can separate: both are (unbarred)²
    # and both are `2₀₀₁` spatially, but they land on opposite SU(2) elements.
    sg75 = spacegroup(75, Val(3), Val(true))
    f = findfirst(d -> seitz(d) == "4₀₀₁⁺", sg75)
    m = findfirst(d -> seitz(d) == "4₀₀₁⁻", sg75)
    @test seitz(sg75[f]^2) == "2₀₀₁"
    @test seitz(sg75[m]^2) == "ᵈ2₀₀₁"

    @test isbarred(DSymOperation(S"x,y,z", one(SU2))) == false

    # A two-fold about an untabulated direction is outside Altmann's tables; we fall back
    # to our own convention there, which must still be a consistent one — exactly one of
    # `u`, `-u` barred, whatever the axis.
    for v in (SVector(1.0, 2.0, 3.0), SVector(-1.0, 0.3, 0.0), SVector(0.0, 0.0, -2.0))
        n = normalize(v)
        uₙ = SU2(-im*n[3], -(n[2] + im*n[1]))
        @test isbarred(uₙ) != isbarred(-uₙ)
    end
end

@testset "Composition is independent of the setting" begin
    # The barring must survive a change of basis, and composition must remain correct there:
    # a lookup keyed on the rotation part could not do this, since primitivizing changes it.
    sgnum = 225
    dsg  = spacegroup(sgnum, Val(3), Val(true))
    cntr = centering(sgnum, 3)
    prim = [DSymOperation(primitivize(SymOperation(d), cntr), su2(d)) for d in dsg]

    @test all(isbarred(p) == isbarred(d) for (p, d) in zip(prim, dsg))  # barring survives
    @test seitz(prim[2]) isa String              # `seitz` must not need the rotation part
    @test all(findfirst(q -> isapprox(q, a*b), prim) !== nothing for a in prim, b in prim)
    @test all(isbarred(dsg[i]*dsg[j]) == isbarred(prim[i]*prim[j])
              for i in eachindex(dsg), j in eachindex(dsg))
end

end # @testset "Double groups"
