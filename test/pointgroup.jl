using Crystalline, Test
using Crystalline: corep_orthogonality_factor
using LinearAlgebra: dot

@testset "Point groups" begin

@testset "Find a point group for every space group" begin
    for D in 1:3
        for sgnum in 1:MAX_SGNUM[D]
            G = spacegroup(sgnum)
            pg = Crystalline.find_parent_pointgroup(G)
            @test pg !== nothing
        end
    end
end

@testset "Matching operator-sorting between symmorphic space groups & point groups" begin
    for D in 1:3
        for sgnum in 1:MAX_SGNUM[D]
            G = spacegroup(sgnum)
            # get isogonal point group of G: strip away translational parts and retain 
            # unique rotation parts; does not change ordering (except for stripping out
            # non-unique rotational parts) and retains conventional lattice vectors
            isogonal_G = pointgroup(G)
            # find a matching point group from those stored obtained directly from
            # pointgroup(iuclab, D), via find_parent_pointgroup
            pg = Crystalline.find_parent_pointgroup(G)
            # compare operator sorting and setting; equivalent in the database
            @test pg == isogonal_G
        end
    end
end

@testset "Schoenflies notation: space vs. point group" begin
    for sgnum in 1:MAX_SGNUM[3]
        sg = spacegroup(sgnum, Val(3))
        pg = Crystalline.find_parent_pointgroup(sg)
        sglab = schoenflies(sgnum)
        pglab = schoenflies(pg)
        @test pglab == sglab[1:lastindex(pglab)]
    end
end

@testset "Generators" begin
    for (Dᵛ, gtype) in ((Val(1), PointGroup{1}), (Val(2), PointGroup{2}), (Val(3), PointGroup{3}))
        D = typeof(Dᵛ).parameters[1]
        for iuclab in Crystalline.PG_IUCs[D]
            ops1 = sort!(pointgroup(iuclab, Dᵛ))
            ops2 = sort!(generate(generators(iuclab, gtype)))
            @test ops1 ≈ ops2
        end
        for pgnum in 1:length(Crystalline.PG_NUM2IUC[D])
            ops1 = sort!(pointgroup(pgnum, Dᵛ))
            ops2 = sort!(generate(generators(pgnum, gtype)))
            @test ops1 ≈ ops2
        end
    end
end

@testset "Irrep orthogonality" begin
## 2ⁿᵈ orthogonality theorem (characters) [automatically includes 1ˢᵗ orthog. theorem also]: 
#       ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦfNₒₚ⁽ᵃ⁾  
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same point group (with i running over the 
# Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements). `f` incorporates a multiplicative factor due to the
# conversionn to "physically real" irreps; see `corep_orthogonality_factor(..)`.
@testset "2ⁿᵈ orthogonality theorem" begin
for D in 1:3
    for pgiuc in PG_IUCs[D]
        pgirs′ = pgirreps(pgiuc, Val(D))
        Nₒₚ = order(first(pgirs′))
        # check both "ordinary" irreps and "physically real" irreps (coreps)
        for irtype in (identity, realify)
            pgirs = irtype(pgirs′)
            for (a, pgir⁽ᵃ⁾) in enumerate(pgirs)
                χ⁽ᵃ⁾ = characters(pgir⁽ᵃ⁾)

                f = corep_orthogonality_factor(pgir⁽ᵃ⁾)
                for (β, pgir⁽ᵝ⁾) in enumerate(pgirs)
                    χ⁽ᵝ⁾ = characters(pgir⁽ᵝ⁾)
                    # ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ == δₐᵦf|g|
                    @test (dot(χ⁽ᵃ⁾, χ⁽ᵝ⁾) ≈ (a==β)*f*Nₒₚ)  atol=1e-12
                end
            end
        end
    end
end
end # @testset "2ⁿᵈ orthogonality theorem"

## Great orthogonality theorem of irreps:
#       ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ = δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same point group (with i running over the
# Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements)
# NB: cannot test this for coreps (see notes in `test/irreps_reality.jl`)
@testset "Great orthogonality theorem" begin
    αβγ = nothing
    for D in 1:3
        for pgiuc in PG_IUCs[D]
            pgirs = pgirreps(pgiuc, Val(D))
            Nₒₚ = order(first(pgirs))
            for (a, pgir⁽ᵃ⁾) in enumerate(pgirs) 
                D⁽ᵃ⁾   = pgir⁽ᵃ⁾(αβγ)      # vector of irreps in (a)
                dim⁽ᵃ⁾ = size(first(D⁽ᵃ⁾),1)

                for (β, pgir⁽ᵝ⁾) in enumerate(pgirs)
                    D⁽ᵝ⁾ = pgir⁽ᵝ⁾(αβγ)    # vector of irreps in (β)
                    dim⁽ᵝ⁾ = size(first(D⁽ᵝ⁾),1)
                    δₐᵦ = (a==β)
                    for n in Base.OneTo(dim⁽ᵃ⁾), j in Base.OneTo(dim⁽ᵝ⁾)     # rows of each irrep
                        δₐᵦδₙⱼ = δₐᵦ*(n==j)
                        for m in Base.OneTo(dim⁽ᵃ⁾), k in Base.OneTo(dim⁽ᵝ⁾) # cols of each irrep
                            δₐᵦδₙⱼδₘₖ = δₐᵦδₙⱼ*(m==k)

                            # compute ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ
                            g_orthog = sum(conj(D⁽ᵃ⁾[i][n,m])*D⁽ᵝ⁾[i][j,k] for i in Base.OneTo(Nₒₚ)) 

                            # test equality to δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
                            @test isapprox(g_orthog, δₐᵦδₙⱼδₘₖ*Nₒₚ/dim⁽ᵃ⁾, atol=1e-12)
                        end
                    end
                end
            end
        end
    end
end # @testset "Great orthogonality theorem"
end # @testset "Irrep orthogonality"
end # @testset "Point groups