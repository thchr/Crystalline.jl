using Crystalline, Test

if !isdefined(Main, :LGIRS′)
    include(joinpath(@__DIR__, "../build/ParseIsotropy.jl"))
    using .ParseIsotropy
    LGIRS′ = parselittlegroupirreps() # parsed explicitly from ISOTROPY's data (see ParseIsotropy)
end
if !isdefined(Main, :LGIRS)
    LGIRS  = lgirreps.(1:MAX_SGNUM[3], Val(3)) # loaded from our saved .jld2 files
end

@testset "Test equivalence of parsed and loaded LGIrreps" begin
    addition_klabs = ("WA", "HA", "KA", "PA")
    for sgnum in 1:MAX_SGNUM[3]
        lgirsd′ = LGIRS′[sgnum] # parsed variant
        lgirsd  = LGIRS[sgnum]  # loaded variant; may contain more irreps than parsed ISOTROPY (cf. `manual_lgirrep_additions.jl`)

        # test we have same k-points in `lgirsd` and `lgirsd′` (or, possibly more in
        # k-points `lgirsd`, if it contains manual additions)
        klabs  = sort!(collect(keys(lgirsd)))
        klabs′ = sort!(collect(keys(lgirsd′)))
        if klabs == klabs′
            @test length(lgirsd) == length(lgirsd′)
        else
            # in case loaded irreps contain manually added 'xA'-type k-points (either "HA",
            # "HA", "KA", or "PA", see `addition_klabs`) from Φ-Ω
            extra_klabidxs = findall(klab -> klab ∉ klabs′, klabs)::Vector{Int}
            extra_klabs = klabs[extra_klabidxs]
            @test length(lgirsd) > length(lgirsd′) && all(∈(addition_klabs), extra_klabs)
        end

        for (kidx, (klab, lgirs)) in enumerate(lgirsd)
            if !haskey(lgirsd′, klab) && klab ∈ addition_klabs
                # `lgirs` is a manual addition that is not in ISOTROPY; nothing to compare
                continue
            end
            
            lgirs′ = lgirsd′[klab]
            @test length(lgirs) == length(lgirs′)
            for (iridx, lgir) in enumerate(lgirs)
                lgir′ = lgirs′[iridx]
                # test that labels agree
                @test label(lgir) == label(lgir′)
                # test that little groups agree
                @test position(lgir) ≈ position(lgir′)
                @test all(operations(lgir) .== operations(lgir′))
                # test that irreps agree
                for αβγ in (nothing, Crystalline.TEST_αβγ)
                    @test lgir(αβγ) == lgir′(αβγ)
                end
            end
        end
    end
end