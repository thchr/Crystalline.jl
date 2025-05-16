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