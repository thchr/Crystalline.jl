using Test
using Crystalline: safe_idxmax_antichain

@testset "Fourier lattices" begin
    for (sgnum, D) in ((5,2), (110, 3))    # centering 'c' in 2D and 'I' in 3D (both body-centered)
        cntr = centering(sgnum, D)

        flat = levelsetlattice(sgnum, D) # Fourier lattice in basis of Rs (rectangular)
        flat = modulate(flat)            # modulate the lattice coefficients

        flat′ = primitivize(flat, cntr)  # Fourier lattice in basis of Rs′ (oblique)

        # test that a primitivize -> conventionalize cycle leaves the lattice unchanged
        @test flat ≈ conventionalize(flat′, cntr)
    end

    # convert-to-integer bug found by Ali (conversion to integer-vectors bugged for 
    # ntuple(_->i, Val(3)) with i > 1)
    sgnum = 146
    flat = levelsetlattice(sgnum, Val(3), ntuple(_->2, Val(3)))
    @test primitivize(flat, centering(sgnum)) isa typeof(flat) # = doesn't throw

    # `isoval2filling` and `filling2isoval`
    Norbits = length(flat.orbits)
    mflat = modulate(flat, collect(ComplexF64, Norbits:-1:1)./Norbits)
    for design_fill in 0.2:0.2:0.8
        iso = filling2isoval(flat, design_fill)
        fill = isoval2filling(flat, iso)
        @test fill ≈ design_fill atol=1e-3
        @test isoval2filling(flat, iso) ≈ fill atol=1e-3
        @test filling2isoval(flat, fill) ≈ iso atol=1e-3
    end
end

# the space groups whose default truncation exceeds the historical `ntuple(i->2, D)`
const HIGH_IDXMAX_SGNUMS = Dict(197=>3, 201=>3, 204=>3, 207=>3, 211=>3,
                                196=>4, 202=>4, 210=>4, 222=>4,
                                203=>5, 209=>5, 219=>5, 226=>5, 228=>5)

@testset "Symmetry-safe Fourier truncations" begin
    @testset "`default_idxmax` is symmetry-safe" begin
        # central invariant: the default truncation must never be symmetry-unsafe
        for D in 1:3
            Dᵛ = Val(D)
            for sgnum in 1:MAX_SGNUM[D]
                idxmax = default_idxmax(sgnum, Dᵛ)
                @test is_symmetry_safe(sgnum, idxmax, Dᵛ)
                @test allequal(idxmax)               # isotropic
                @test all(≥(2), idxmax)              # floored at the historical default
                # ... and it is the *smallest* isotropic safe truncation above that floor
                N = first(idxmax)
                if N > 2
                    @test !is_symmetry_safe(sgnum, ntuple(_->N-1, Dᵛ), Dᵛ)
                end
            end
        end
    end

    @testset "defaults differ from `(2,…,2)` only for tabulated 3D groups" begin
        for D in 1:2, sgnum in 1:MAX_SGNUM[D] # 1D & 2D are unaffected
            @test default_idxmax(sgnum, Val(D)) == ntuple(_->2, Val(D))
        end
        for sgnum in 1:MAX_SGNUM[3]
            N = get(HIGH_IDXMAX_SGNUMS, sgnum, 2)
            @test default_idxmax(sgnum, Val(3)) == (N, N, N)
            # groups needing a higher truncation are exactly those unsafe at `(2,2,2)`
            @test is_symmetry_safe(sgnum, (2,2,2), Val(3)) == (N == 2)
        end
    end

    @testset "safety is monotonic in `idxmax`" begin
        # `is_symmetry_safe` must be upward-closed: raising any component keeps it safe
        for D in 1:3, sgnum in 1:MAX_SGNUM[D]
            Dᵛ = Val(D)
            for idxmax in Iterators.product(ntuple(_->1:3, Dᵛ)...)
                is_symmetry_safe(sgnum, idxmax, Dᵛ) || continue
                for k in 1:D
                    idxmax′ = ntuple(j -> j == k ? idxmax[j]+1 : idxmax[j], Dᵛ)
                    @test is_symmetry_safe(sgnum, idxmax′, Dᵛ)
                end
            end
        end
    end

    @testset "anisotropic truncations" begin
        # the check must not collapse to the isotropic minimum: sg 27 is unsafe at (1,1,1)
        # (it generically gains an inversion center) but safe at the anisotropic (2,1,1)
        @test !is_symmetry_safe(27, (1,1,1))
        @test is_symmetry_safe(27, (2,1,1))
        @test default_idxmax(27, Val(3)) == (2,2,2)
        # tabulated minimal truncations are minimal: lowering any component is unsafe
        for D in 1:3, sgnum in 1:MAX_SGNUM[D]
            Dᵛ = Val(D)
            for p in safe_idxmax_antichain(sgnum, Dᵛ)
                @test is_symmetry_safe(sgnum, p, Dᵛ)
                for k in 1:D
                    p[k] == 1 && continue
                    @test !is_symmetry_safe(sgnum, ntuple(j -> j==k ? p[j]-1 : p[j], Dᵛ), Dᵛ)
                end
            end
        end
    end

    @testset "`minimal_idxmax`" begin
        for D in 1:3, sgnum in 1:MAX_SGNUM[D]
            Dᵛ = Val(D)
            m = minimal_idxmax(sgnum, Dᵛ)

            # every returned minimum is safe, and is minimal (lowering a component is not)
            for p in m.general
                @test is_symmetry_safe(sgnum, p, Dᵛ)
                for k in 1:D
                    p[k] == 1 && continue
                    p′ = ntuple(j -> j==k ? p[j]-1 : p[j], Dᵛ)
                    @test !is_symmetry_safe(sgnum, p′, Dᵛ)
                end
            end
            # the minima are mutually incomparable (an antichain): none dominates another
            for p in m.general, q in m.general
                p === q && continue
                @test !all(k -> p[k] ≤ q[k], 1:D)
            end

            # `isotropic` agrees with its derivation from `general` ...
            @test allequal(m.isotropic)
            @test is_symmetry_safe(sgnum, m.isotropic, Dᵛ)
            N = first(m.isotropic)
            @test N == minimum(maximum, m.general)
            N > 1 && @test !is_symmetry_safe(sgnum, ntuple(_->N-1, Dᵛ), Dᵛ)
            # ... and `default_idxmax` is just `isotropic`, floored at 2
            @test default_idxmax(sgnum, Dᵛ) == max.(2, m.isotropic)
        end

        # the returned vector must be a copy: mutating it cannot corrupt the tabulated data
        m = minimal_idxmax(228, Val(3))
        push!(empty!(m.general), (0,0,0))
        @test minimal_idxmax(228, Val(3)).general == safe_idxmax_antichain(228, Val(3))
        @test length(minimal_idxmax(228, Val(3)).general) == 6

        # an anisotropic minimum can be much cheaper than the isotropic one
        @test length(getorbits(levelsetlattice(228, Val(3), (1,3,5)))) <
              length(getorbits(levelsetlattice(228, Val(3), (5,5,5))))

        # `D::Integer` convenience method, and rejection of invalid dimensions
        @test minimal_idxmax(16, 2) == minimal_idxmax(16, Val(2))
        @test_throws DomainError minimal_idxmax(1, Val(4))
    end

    @testset "enantiomorphic partners agree" begin
        # enantiomorphic pairs are related by an improper transformation, so they must
        # require identical truncations
        for (sgnum¹, sgnum²) in ENANTIOMORPHIC_PAIRS
            @test default_idxmax(sgnum¹, Val(3)) == default_idxmax(sgnum², Val(3))
        end
    end

    @testset "`levelsetlattice` integration" begin
        # `:auto` is the default and warns for no group; explicit unsafe values do warn
        @test levelsetlattice(228, Val(3)) == levelsetlattice(228, Val(3), :auto)
        @test (@test_logs levelsetlattice(228, Val(3))) isa Crystalline.UnityFourierLattice{3}
        @test_logs (:warn,) levelsetlattice(228, Val(3), (2,2,2))
        @test_logs (:warn,) levelsetlattice(27, Val(3), (1,1,1))
        @test_logs levelsetlattice(27, Val(3), (2,1,1)) # safe, despite anisotropic
        @test_throws DomainError levelsetlattice(1, Val(3), :automatic)

        # invalid `sgnum` is rejected in every dimension (incl. 1D & the upper bounds)
        for D in 1:3
            @test_throws DomainError levelsetlattice(0, Val(D))
            @test_throws DomainError levelsetlattice(MAX_SGNUM[D]+1, Val(D))
        end

        # the `D::Integer` convenience methods agree with the `Val{D}` methods
        @test levelsetlattice(16, 2) == levelsetlattice(16, Val(2))
        @test levelsetlattice(16, 2, :auto) == levelsetlattice(16, Val(2))
        @test levelsetlattice(16, 2, (2,2)) == levelsetlattice(16, Val(2), (2,2))
        @test_throws DomainError levelsetlattice(16, 2, (2,2,2))

        # backwards compatibility: unaffected groups are unchanged by the new default
        @test length(getorbits(levelsetlattice(230, Val(3)))) == 3 # cf. `docs/src/lattices.md`
        @test levelsetlattice(230, Val(3)) == levelsetlattice(230, Val(3), (2,2,2))
        @test levelsetlattice(16, Val(2)) == levelsetlattice(16, Val(2), (2,2))
    end
end