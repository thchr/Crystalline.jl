using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS = lgirreps.(1:MAX_SGNUM[3], Val(3)) # loaded from our saved .jld2 files
end

@testset "k-vectors required by BandRepSet analysis" begin
allpaths = false
spinful  = false
debug = false

@testset "Complex (no TR) irreps" begin
# --- test complex-form irreps (not assuming time-reversal symmetry) ---
for (sgnum, lgirsd) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=false)
    irlabs_BRS = irreplabels(BRS)
    klabs_BRS = klabels(BRS)

    irlabs_ISO = [Crystalline.formatirreplabel(label(lgir)) for lgirs in values(lgirsd) for lgir in lgirs]
    klabs_ISO = keys(lgirsd)

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab = klabel(irlab_BRS)
        @test irlab_BRS ∈ irlabs_ISO
        if debug && irlab_BRS ∉ irlabs_ISO
            @info "Cannot find complex irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
        end

        # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
        kidx_BRS = findfirst(==(klab), klabs_BRS)
        @test BRS.kvs[kidx_BRS] == kvec(first(lgirsd[klab]))
        if debug && BRS.kvs[kidx_BRS] ≠ kvec(first(lgirsd[klab]))
            println("Different definitions of k-point labels in space group ", sgnum)
            println("   BRS, ", klab, ": ",string(BRS.kvs[kidx_BRS]))
            println("   ISO, ", klab, ": ",string(kvec(first(lgirsd[klab]))), "\n")
        end

    end
end
end

@testset "Physically irreducible irreps/co-reps (with TR)" begin
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
for (sgnum, lgirsd) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
    irlabs_BRS = irreplabels(BRS)
    klabs_BRS = klabels(BRS)

    irlabs_ISO = Vector{String}()
    realirlabs_ISO = Vector{String}()
    klabs_ISO = keys(lgirsd)
    klabs_ISO = Vector{String}(undef, length(lgirsd))
    for lgirs in values(lgirsd)
        append!(irlabs_ISO,     [label(lgir) for lgir in lgirs])
        append!(realirlabs_ISO, label.(realify(lgirs)))
    end
    irlabs_ISO = Crystalline.formatirreplabel.(irlabs_ISO)
    realirlabs_ISO = Crystalline.formatirreplabel.(realirlabs_ISO)

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab = klabel(irlab_BRS)
        @test irlab_BRS ∈ realirlabs_ISO
        if debug && irlab_BRS ∉ realirlabs_ISO
            @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
        end

        # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
        kidx_BRS = findfirst(==(klab), klabs_BRS)
        @test BRS.kvs[kidx_BRS] == kvec(first(lgirsd[klab]))
    end
end
end
end

@testset "BandRepSet and BandRep" begin
    BRS = bandreps(230)
    # iterated concatenation of vectors of `BRS` should give `matrix`
    @test matrix(BRS) == matrix(BRS; includedim=true) == hcat(BRS...)
    # length of BandRep as vectors should be = number of irreps + 1 (i.e. includes filling)
    @test length(BRS[1]) == length(BRS[1].irvec)+1
    @test BRS[1] == vcat(BRS[1].irvec, dim(BRS[1]))
end