using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS = get_lgirreps.(1:MAX_SGNUM[3], Val(3)) # loaded from our saved .jld2 files
end

@testset "k-vectors required by BandRepSet analysis" begin
allpaths = false
spinful  = false
verbose = false

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
        if irlab_BRS ∉ irlabs_ISO
            if verbose
                @info "Cannot find complex irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false
            kidx_BRS = findfirst(==(klab), klabs_BRS)

            # test that for each of the (P,K,W,H)A k-label variants, that the associated k-vector is 
            # just equal to minus the kvector in the (P,K,W,H) variant (i.e. without the 'A' postscript)
            @test -BRS.kvs[kidx_BRS] == kvec(first(lgirsd[klab[1:1]]))

            # TODO: to get the irreps of these variants, we need to follow the prescription 
            # detailed in CDML p. 69-73 (though that won't work, presumably, for sgnum=205)
            # This is the subject of Cracknell & Davies 1976b (On the completeness of tables
            # of irreducible representations of the classical space groups)
        else
            kidx_BRS = findfirst(==(klab), klabs_BRS)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsd[klab]))
            if verbose
                if BRS.kvs[kidx_BRS] ≠ kvec(first(lgirsd[klab]))
                    println("Different definitions of k-point labels in space group ", sgnum)
                    println("   BRS, ", klab, ": ",string(BRS.kvs[kidx_BRS]))
                    println("   ISO, ", klab, ": ",string(kvec(first(lgirsd[klab]))))
                    println()
                end
            end
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
        if irlab_BRS ∉ realirlabs_ISO
            if verbose
                @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false

        else
            kidx_BRS = findfirst(==(klab), klabs_BRS)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsd[klab]))
        end
    end
end
end
end