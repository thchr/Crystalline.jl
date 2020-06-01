using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS = get_all_lgirreps(Val(3))  # loaded from our saved .jld2 files
end

@testset "k-vectors required by BandRepSet analysis" begin
allpaths = false
spinful  = false
verbose = false

@testset "Complex (no TR) irreps" begin
# --- test complex-form irreps (not assuming time-reversal symmetry) ---
for (sgnum, lgirsvec) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=false)
    irlabs_BRS = irreplabels(BRS)
    klabs_BRS = klabels(BRS)

    irlabs_ISO = [Crystalline.formatirreplabel(label(lgir)) for lgirs in lgirsvec for lgir in lgirs]
    klabs_ISO = [klabel(lgirs[1]) for lgirs in lgirsvec]

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab_BRS = klabel(irlab_BRS)
        if irlab_BRS ∉ irlabs_ISO
            if verbose
                @info "Cannot find complex irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false
            kidx_BRS = findfirst(==(klab_BRS), klabs_BRS)
            kidx_ISO_related = findfirst(==(klab_BRS[1:1]), klabs_ISO)

            # test that for each of the (P,K,W,H)A k-label variants, that the associated k-vector is 
            # just equal to minus the kvector in the (P,K,W,H) variant (i.e. without the 'A' postscript)
            @test -BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO_related]))

            # TODO: to get the irreps of these variants, we need to follow the prescription 
            # detailed in CDML p. 69-73 (though that won't work, presumably, for sgnum=205)
            # This is the subject of Cracknell & Davies 1976b (On the completeness of tables
            # of irreducible representations of the classical space groups)
        else
            kidx_BRS = findfirst(==(klab_BRS), klabs_BRS)
            kidx_ISO = findfirst(==(klab_BRS), klabs_ISO)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO]))
            if verbose
                if BRS.kvs[kidx_BRS] ≠ kvec(first(lgirsvec[kidx_ISO]))
                    println("Different definitions of k-point labels in space group ", sgnum)
                    println("   BRS, ", klab_BRS, ": ",string(BRS.kvs[kidx_BRS]))
                    println("   ISO, ", klabs_ISO[kidx_ISO], ": ",string(kvec(first(lgirsvec[kidx_ISO]))))
                    println()
                end
            end
        end
    end
end
end

@testset "Physically irreducible irreps/co-reps (with TR)" begin
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
for (sgnum, lgirsvec) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
    irlabs_BRS = irreplabels(BRS)
    klabs_BRS = klabels(BRS)

    irlabs_ISO = Vector{String}()
    realirlabs_ISO = Vector{String}()
    klabs_ISO = Vector{String}(undef, length(lgirsvec))
    for (kidx, lgirs) in enumerate(lgirsvec)
        append!(irlabs_ISO,     [label(lgir) for lgir in lgirs])
        append!(realirlabs_ISO, label.(realify(lgirs)))
        klabs_ISO[kidx] = klabel(first(lgirs))
    end
    irlabs_ISO = Crystalline.formatirreplabel.(irlabs_ISO)
    realirlabs_ISO = Crystalline.formatirreplabel.(realirlabs_ISO)

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab_BRS = klabel(irlab_BRS)
        if irlab_BRS ∉ realirlabs_ISO
            if verbose
                @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false

        else
            kidx_BRS = findfirst(==(klab_BRS), klabs_BRS)
            kidx_ISO = findfirst(==(klab_BRS), klabs_ISO)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO]))
        end
    end
end
end
end