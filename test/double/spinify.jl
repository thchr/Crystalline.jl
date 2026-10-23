using Crystalline, Test
using Crystalline: MAX_SGNUM

@testset "Spinifying irreps and band representations" begin
    @testset "Textbook: p-orbitals in m-3m split into j = ½ and j = 3/2" begin
        pgirs  = pgirreps("m-3m", Val(3))
        pgirsˢ = pgirreps("m-3m", Val(3); spinful=Val(true))
        T₁ᵤ = pgirs[findfirst(ir -> mulliken(ir) == "T₁ᵤ", pgirs)]
        @test mulliken.(spinify(T₁ᵤ, pgirsˢ)) == ["E₁ᵤˢ", "Fᵤˢ"] # Γ₆⁻ ⊕ Γ₈⁻ (Koster)
        # the spinification of every irrep doubles its dimension
        for ir in pgirs
            @test sum(irdim, spinify(ir, pgirsˢ)) == 2irdim(ir)
        end
    end

    @testset "Little group irreps (incl. nonsymmorphic & k-lines)" begin
        for sgnum in (4, 19, 62, 194, 227, 230), timereversal in (false, true)
            lgirsd  = lgirreps(sgnum, Val(3))
            lgirsdˢ = lgirreps(sgnum, Val(3); spinful=Val(true))
            for (klab, lgirs) in lgirsd
                lgirsˢ = lgirsdˢ[klab]
                if timereversal
                    lgirs, lgirsˢ = realify(lgirs), realify(lgirsˢ)
                end
                for lgir in lgirs
                    lgirsˢ′ = spinify(lgir, lgirsˢ)
                    @test lgirsˢ′ isa Vector{DLGIrrep{3}}
                    @test sum(irdim, lgirsˢ′) == 2irdim(lgir)
                end
            end
        end
    end

    @testset "Irreps at different k-points or Wyckoff positions are rejected" begin
        # symmorphic: Γ and R have isomorphic little groups (m-3m) with identical operations,
        # so a mismatch would otherwise go unnoticed
        lgirsd, lgirsdˢ = lgirreps(221, Val(3)), lgirreps(221, Val(3); spinful=Val(true))
        @test_throws ErrorException spinify(lgirsd["Γ"][1], lgirsdˢ["R"])
        @test spinify(lgirsd["Γ"][1], lgirsdˢ["Γ"]) isa Vector{DLGIrrep{3}}
        # nonsymmorphic, & via a symmetry vector whose k-points are misaligned with `brsˢ`
        brs, brsˢ = bandreps(230, Val(3)), bandreps(230, Val(3); spinful=Val(true))
        @test_throws ErrorException spinify(SymmetryVector(brs[1]), reverse(irreps(brsˢ)))
        # site irreps: 1a and 1b in Pm-3m (221) have isomorphic site groups (m-3m)
        sitegs, sitegsˢ = sitegroups(221, Val(3)), sitegroups(221, Val(3); spinful=Val(true))
        i1a = findfirst(g -> Crystalline.positionlabel(g) == "1a", sitegs)
        i1b = findfirst(g -> Crystalline.positionlabel(g) == "1b", sitegsˢ)
        siteir = first(siteirreps(sitegs[i1a]))
        @test_throws ErrorException spinify(siteir, siteirreps(sitegsˢ[i1b]))
        @test spinify(siteir, siteirreps(sitegsˢ[i1a])) isa Vector{DSiteIrrep{3}}
    end

    @testset "Inconsistent inputs are rejected" begin
        # time reversal applied to `irsˢ` but not to `ir`: a complex irrep's spinification
        # cannot be decomposed into co-reps
        pgirs, pgirsˢ = pgirreps("3", Val(3)), pgirreps("3", Val(3); spinful=Val(true))
        @test_throws ErrorException spinify(pgirs[2], realify(pgirsˢ))
        # `brsˢ` missing a spinful band rep that the spinified band rep requires
        brs, brsˢ = bandreps(230, Val(3)), bandreps(230, Val(3); spinful=Val(true))
        i = findfirst(!iszero, spinify(brs[1], brsˢ).coefs)
        @test_throws ErrorException spinify(brs[1], Collection(deleteat!(collect(brsˢ), i)))
    end

    @testset "Site irreps agree with band representations" begin
        for (D, sgnum) in ((2, 10), (2, 17), (3, 62), (3, 183), (3, 221)),
                timereversal in (false, true)
            brs  = bandreps(sgnum, D; timereversal)
            brsˢ = bandreps(sgnum, D; spinful=true, timereversal)
            for br in brs
                siteirsˢ = Collection([brˢ.siteir for brˢ in brsˢ
                                       if position(brˢ) ≈ position(br)])
                siteirsˢ′ = spinify(br.siteir, siteirsˢ)
                cbrˢ = spinify(br, brsˢ)
                @test sort(label.(siteirsˢ′)) == sort([label(brˢ.siteir)
                    for (c, brˢ) in zip(cbrˢ.coefs, brsˢ) for _ in 1:Int(c)])
            end
        end
    end

    @testset "Graphene: SOC can gap the Dirac point at K" begin
        brs  = bandreps(17, Val(2))
        brsˢ = bandreps(17, Val(2); spinful=Val(true))
        # p_z orbitals on the honeycomb sites: (2b|A₁)
        br = brs[findfirst(br -> label(position(br)) == "2b" && label(br.siteir) == "A₁", brs)]
        cbrˢ = spinify(br, brsˢ)
        @test count(!iszero, cbrˢ.coefs) == 1 # a single, elementary spinful band rep
        nˢ = spinify(SymmetryVector(br), brsˢ)
        @test occupation(nˢ) == 4
        # the 2D irrep K₃ of the Dirac point spinifies to 1D+1D+2D irreps (Kane–Mele gap)
        iK = findfirst(==("K"), klabels(nˢ))
        lgirsKˢ = irreps(nˢ)[iK]
        @test sort(irdim.(lgirsKˢ)[findall(!iszero, multiplicities(nˢ)[iK])]) == [1, 1, 2]
    end
    
    @testset "Composite band representations" begin
        brs, brsˢ = bandreps(2, Val(3)), bandreps(2, Val(3); spinful=Val(true))
        cbr = CompositeBandRep([2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], brs)
        cbrˢ = spinify(cbr, brsˢ)
        @test cbrˢ isa CompositeBandRep{3, DLGIrrep{3}, DSiteIrrep{3}}
        @test cbrˢ.coefs == 2spinify(brs[1], brsˢ).coefs - spinify(brs[3], brsˢ).coefs +
                            spinify(brs[16], brsˢ).coefs
        @test cbrˢ == spinify(SymmetryVector(cbr), brsˢ)
    end

    @testset "Every spinful rep reachable by spinifying a spinless rep (D=$D)" for D in 1:3
        for sgnum in 1:MAX_SGNUM[D], timereversal in (true, false)
            brs  = bandreps(sgnum, D; spinful=false, timereversal)
            brsˢ = bandreps(sgnum, D; spinful=true,  timereversal)
            coefsˢ = zeros(Rational{Int}, length(brsˢ))
            for br in brs
                n = SymmetryVector(br)
                nˢ   = spinify(n, brsˢ)  # spinification of `n` as a `SymmetryVector`
                cbrˢ = spinify(br, brsˢ) # -------||------ `br` as a `CompositeBandRep`
                @test nˢ == cbrˢ
                @test occupation(nˢ) == 2occupation(br)
                @test all(c -> c ≥ 0 && isinteger(c), cbrˢ.coefs)
                @test all(≥(0), nˢ)
                coefsˢ += cbrˢ.coefs # summing, intentionally
            end
            # every spinful band rep is reached by spinifying some spinless band rep
            @test all(!iszero, coefsˢ)
        end
    end
end
