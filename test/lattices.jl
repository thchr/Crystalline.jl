using Test

@testset "Fourier lattices" begin
    for (sgnum, D) in ((5,2), (110, 3))    # centering 'c' in 2D and 'I' in 3D (both body-centered)
        cntr = centering(sgnum, D)

        flat = levelsetlattice(sgnum, D) # Fourier lattice in basis of Rs (rectangular)
        flat = modulate(flat)            # modulate the lattice coefficients

        flat′ = primitivize(flat, cntr)  # Fourier lattice in basis of Rs′ (oblique)

        # test that a primitivize -> conventionalize cycle leaves the lattice unchanged
        @test flat ≈ conventionalize(flat′, cntr)
    end

    # convert-to-integer bug found by Ali
    flat = levelsetlattice(sgnum, Val(3), ntuple(_->2, Val(3)))
    @test primitivize(flat, centering(sgnum)) isa UnityFourierLattice{3} # = doesn't throw
end