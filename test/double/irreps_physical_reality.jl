using Crystalline, Test, LinearAlgebra
using Crystalline: check_multtable_vs_ir

datafile = joinpath(pkgdir(Crystalline), "data", "irreps", "pgs", "3d",
                    "irreps_data_spinful.jld2")
if !isfile(datafile)
    @warn "spinful point group irrep data not found; skipping spinful `physical_realify` \
           tests" datafile
else

@testset "`physical_realify` for double-valued irreps" begin

@testset "`timereversal_unitary`" begin
    # spinless: the identity matrix
    pgir = first(physical_realify(pgirreps("4", Val(3))))
    @test timereversal_unitary(pgir) == I(irdim(pgir))

    # spinful: `J = iσʸ ⊗ 𝟙ₙ`, i.e., with the Kramers index as the outer index (and not
    # `𝟙ₙ ⊗ iσʸ`, which would interleave the two)
    dpgirs = physical_realify(pgirreps("m-3m", Val(3); spinful=Val(true)))
    dpgir = argmax(irdim, dpgirs)
    n = irdim(dpgir) ÷ 2
    @test n ≥ 2 # else the two orderings would coincide, and the test would be vacuous
    @test timereversal_unitary(dpgir) == kron([0 1; -1 0], I(n))
    @test timereversal_unitary(dpgir) == [zeros(n,n) I(n); -I(n) zeros(n,n)]

    # an odd-dimensional spinful irrep has no `J`; it must be glued to a partner first
    odd_ir = first(pgirreps("3", Val(3); spinful=Val(true)))
    @test isodd(irdim(odd_ir))
    @test_throws ErrorException timereversal_unitary(odd_ir)
end

@testset "Point group irreps" begin
    for pglab in Crystalline.PG_IUCs[3]
        irs    = pgirreps(pglab, Val(3); spinful=Val(true))
        re_irs = realify(irs)
        irs′   = physical_realify(irs)

        # equivalent to the ordinary coreps, and labelled identically
        @test characters(irs′) ≈ characters(re_irs)
        @test label.(irs′) == label.(re_irs)

        for ir in irs′
            Γ = timereversal_unitary(ir)
            @test all(D -> Γ*conj(D)*Γ' ≈ D, ir.matrices)
            @test all(D -> D'D ≈ I, ir.matrices) # still unitary

            isexplicitlyreal = all(D -> D ≈ real(D), ir.matrices)
            if reality(ir) == PSEUDOREAL
                # `realify` leaves these unglued; they have no explicitly real form
                @test !isexplicitlyreal
            else
                @test isexplicitlyreal
            end
        end

        # applying it again changes nothing
        @test all(((a, b),) -> a.matrices == b.matrices, zip(irs′, physical_realify(irs′)))
    end
end

@testset "Site symmetry irreps" begin
    for sgnum in 1:MAX_SGNUM[3]
        for siteg in sitegroups(spacegroup(sgnum, Val(3); spinful=Val(true)))
            irs    = siteirreps(siteg)
            re_irs = realify(irs)
            irs′   = physical_realify(irs)
            @test characters(irs′) ≈ characters(re_irs)

            mt = MultTable(siteg)
            n = length(siteg) ÷ 2
            for ir in irs′
                Γ = timereversal_unitary(ir)
                @test all(D -> Γ*conj(D)*Γ' ≈ D, ir.matrices)
                @test all(check_multtable_vs_ir(mt, ir))      # still a representation
                @test ir.matrices[n+1:2n] ≈ -ir.matrices[1:n] # still D(Ēg) = -D(g)

                isexplicitlyreal = all(D -> D ≈ real(D), ir.matrices)
                if reality(ir) == PSEUDOREAL
                    @test !isexplicitlyreal
                else
                    @test isexplicitlyreal
                end
            end
        end
    end
end

end # @testset "`physical_realify` for double-valued irreps"
end # if isfile(datafile)
