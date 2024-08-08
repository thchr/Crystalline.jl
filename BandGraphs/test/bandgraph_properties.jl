using Pkg
(dirname(Pkg.project().path) == @__DIR__) || Pkg.activate(@__DIR__)

using Test
using BandGraphs
using Crystalline
using SymmetryBases

using Crystalline: irdim

@testset "Consistency of band graph adjacency matrix blocks" begin
@testset for D in (2,3)
    Dᵛ = Val(D)
    @testset for sgnum in 1:MAX_SGNUM[D]
        # skip 147 and 123 since it takes a long time to do the compatibility basis there
        sgnum ∈ (47, 123) && continue
        checked_multiplicity = false
        for timereversal in (false, true)
            subts = subduction_tables(sgnum, Dᵛ; timereversal)
            sb, _ = compatibility_basis(sgnum, D; timereversal)
            lgirsd = lgirreps(sgnum, Dᵛ)
            timereversal && realify!(lgirsd)
            for _n in sb
                n = SymVector(_n, sb.irlabs, lgirsd)
                bandg = build_subgraphs(n, subts, lgirsd)
                partitions, subgraphs = bandg.partitions, bandg.subgraphs
                
                # band occupation & grand-sums of all subgraphs must be constant & equal
                occupations = [sum(s.A) for s in subgraphs]
                @test all(==(n.μ), occupations)
                A = assemble_adjacency(bandg)
                @test all(b -> iszero(b) || sum(b) == n.μ, A.blocks)

                # column-wise and row-wise sums of adjacency blocks must give corresponding
                # irrep dimensions of columns and rows, respectively
                nonmax_irdims = [irdim.(s.p_nonmax.lgirs) for s in subgraphs]
                max_irdims    = [irdim.(s.p_max.lgirs) for s in subgraphs]
                colwise_sums = [vec(sum(s.A; dims=1)) for s in subgraphs] # sums per column (across rows)
                rowwise_sums = [vec(sum(s.A; dims=2)) for s in subgraphs] # sums per row (across columns)

                @test colwise_sums == nonmax_irdims
                @test rowwise_sums == max_irdims

                # every non-maximal irrep must appear at least twice
                if !checked_multiplicity # no need to check every symmetry vector; one is enough
                    nonmax_multiplicity = [count(s->s.p_nonmax == p, subgraphs) 
                                                        for p in partitions if !p.maximal]
                    @test all(≥(2), nonmax_multiplicity)
                    if !isempty(subgraphs)
                        @test !isempty(nonmax_multiplicity)
                    end
                    check_multiplicity = true
                end
            end
        end
    end
end
end # @testset
nothing