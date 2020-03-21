using SGOps, Test

@testset "Find a point group for every space group" begin
    for D = 1:3
        for sgnum in 1:MAX_SGNUM[D]
            G = spacegroup(sgnum)
            pg = SGOps.find_parent_pointgroup(G)
            @test pg !== nothing
        end
    end
end

## Great orthogonality theorem of irreps:
#       ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ = δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same little group (with 
# i running over the Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements)
@testset "Great orthogonoality theorem (point group irreps)" begin
    αβγ = nothing
    for D in 1:3
        for pgiuc in SGOps.PGS_IUCs[D]
            pgirs = get_pgirreps(pgiuc, Val(D))
            Nₒₚ = length(operations(group(first(pgirs))))
            for (a, pgir⁽ᵃ⁾) in enumerate(pgirs) 
                D⁽ᵃ⁾ = irreps(pgir⁽ᵃ⁾,αβγ)        # vector of irreps in (a)
                dim⁽ᵃ⁾ = size(first(D⁽ᵃ⁾),1)

                for (β, pgir⁽ᵝ⁾) in enumerate(pgirs)
                    D⁽ᵝ⁾ = irreps(pgir⁽ᵝ⁾,αβγ)    # vector of irreps in (β)
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
end