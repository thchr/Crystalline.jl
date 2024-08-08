using Crystalline, Test

@testset "Conjugacy classes" begin
@testset "Identical number of conjugacy classes and irreps" begin
    # The number of conjugacy classes should equal the number of inequivalent irreps
    # (see e.g. Inui, Section 4.6.1 or Bradley & Cracknel, Theorem 1.3.10) for "linear"
    # or ordinary irreps (i.e. not ray/projective irreps)
    @testset "Point groups" begin
        for Dᵛ in (Val(1), Val(2), Val(3))
            D = typeof(Dᵛ).parameters[1]
            for pgiuc in Crystalline.PG_IUCs[D]
                pg = pointgroup(pgiuc, Dᵛ)
                conj_classes = classes(pg)
                pgirs = pgirreps(pgiuc, Dᵛ)
                @test length(pgirs) == length(conj_classes)
            end
        end
    end

    @testset "Little groups" begin
        for Dᵛ in (Val(1), Val(2), Val(3))
            D = typeof(Dᵛ).parameters[1]
            for sgnum in 1:MAX_SGNUM[D]
                lgirsd = lgirreps(sgnum, Dᵛ)
                for (klab, lgirs) in lgirsd
                    lg = group(first(lgirs))
                    includes_ray_reps = any(lgir->Crystalline.israyrep(lgir)[1], lgirs)
                    if !includes_ray_reps
                        # the number of conjugacy classes is not generally equal to the
                        # number of inequivalent irreps if the irreps include ray (also
                        # called projective) representations as discussed by Inui, Sec. 5.4
                        # (instead, the number of "ray classes" is equal to the number of
                        # inequivalent irreps); so we only test if there are no ray reps
                        conj_classes = classes(lg)

                        @test length(lgirs) == length(conj_classes)
                    end
                end
            end
        end
    end
end

@testset "Primitive/conventional setting vs. conjugacy classes" begin
    # The number of conjugacy classes should be invariant of whether we choose to compute it
    # in a primitive or conventional setting (provided we specify the `cntr` argument
    # correctly to `classes`):
    for Dᵛ in (Val(1), Val(2), Val(3))
        D = typeof(Dᵛ).parameters[1]
        for sgnum in 1:MAX_SGNUM[D]
            # space groups
            sgᶜ = spacegroup(sgnum, Dᵛ) # conventional setting
            sgᵖ = primitivize(sgᶜ)      # primitive setting
            @test length(classes(sgᶜ)) == length(classes(sgᵖ, nothing))
            # little groups
            lgsᶜ = littlegroups(sgnum, Dᵛ)
            for lgᶜ in values(lgsᶜ)
                lgᵖ = primitivize(lgᶜ)
                @test length(classes(lgᶜ)) == length(classes(lgᵖ, nothing))
            end
        end
    end
end


# TODO
# @testset "Conjugacy classes and characters" begin
    # the characters of elements in the same conjugacy class must be identical for ordinary
    # (non-ray) irreps. I.e., for two conjugate elements ``a`` and ``b`` in ``G``, there
    # exists a ``g ∈ G`` s.t. ``gag⁻¹ = b``. As a consequence: ``χ(b) = χ(gag⁻¹) = 
    # Tr(D(gag⁻¹)) = Tr(D(g)D(a)D(g⁻¹)) = Tr(D(a)D(g⁻¹)D(g)) = Tr(D(a)) = χ(a)``.
# end
end # @testset "Conjugacy classes"

@testset "Abelian groups" begin 
    # Abelian groups only have one-dimensional irreps (unless there are ray-irreps, then the
    # rule doesn't hold in general)
    @testset "Point groups" begin
        for Dᵛ in (Val(1), Val(2), Val(3))
            D = typeof(Dᵛ).parameters[1]
            for pgiuc in Crystalline.PG_IUCs[D]
                pg = pointgroup(pgiuc, Dᵛ)
                if is_abelian(pg)
                    pgirs = pgirreps(pgiuc, Dᵛ)
                    @test all(pgir -> Crystalline.irdim(pgir) == 1, pgirs)
                end
            end
        end
    end

    @testset "Little groups" begin
        for Dᵛ in (Val(1), Val(2), Val(3))
            D = typeof(Dᵛ).parameters[1]
            for sgnum in 1:MAX_SGNUM[D]
                lgirsd = lgirreps(sgnum, Dᵛ)
                for (klab, lgirs) in lgirsd
                    lg = group(first(lgirs))
                    if is_abelian(lg)
                        includes_ray_reps = any(lgir->Crystalline.israyrep(lgir)[1], lgirs)
                        if !includes_ray_reps
                            @test all(lgir -> Crystalline.irdim(lgir) == 1, lgirs)
                        end
                    end
                end
            end
        end
    end
end # @testset "Abelian groups"
