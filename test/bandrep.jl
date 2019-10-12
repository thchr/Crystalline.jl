using SGOps, Test

if !isdefined(Main, :LGIRS)
    LGIRS = get_all_lgirreps(3)  # loaded from our saved .jld2 files
end

@testset "k-vectors required by BandRepSet analysis" begin
allpaths = false
spinful  = false
showmissing = true

@testset "Complex (no TR) irreps" begin
# --- test complex-form irreps (not assuming time-reversal symmetry) ---
for (sgnum, lgirsvec) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths, spinful, "Elementary")
    irlabs_BRS = BRS.irreplabs
    klabs_BRS = BRS.klabs

    irlabs_ISO = [SGOps.formatirreplabel(label(lgir)) for lgirs in lgirsvec for lgir in lgirs]
    klabs_ISO = [klabel(lgirs[1]) for lgirs in lgirsvec]

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab_BRS = klabel(irlab_BRS)
        # unfortunately, in our ISO dataset the general point is Ω; in 
        # the band rep set it is GP: put a bandaid on that for now
        klab_BRS′ = klab_BRS == "GP" ? "Ω" : klab_BRS
        irlab_BRS′ = occursin("GP", irlab_BRS) ? "Ω₁" : irlab_BRS
        if irlab_BRS′ ∉ irlabs_ISO
            if showmissing
                @info "Cannot find complex irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false
            kidx_BRS = findfirst(x->x==klab_BRS, klabs_BRS)
            kidx_ISO_related = findfirst(x->x==klab_BRS′[1:1], klabs_ISO)

            # test that for each of the (P,K,W,H)A k-label variants, that the associated k-vector is 
            # just equal to minus the kvector in the (P,K,W,H) variant (i.e. without the 'A' postscript)
            @test -BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO_related]))

            # Note that (Z, ZA) points are not simply the plus/minus pairs of 
            # eachother; rather, Z≡[α,0.5,0] and ZA≡[0.5,α,0]; this affects 
            # space groups 195, 198, 200, and 201, but ZA only features in the
            # allpaths bandreps. 

            # TODO:
            # ... to get the irreps of these variants, we need to follow the prescription 
            # detailed in CDML p. 69-73 (though that won't work, presumably, for sgnum = 205)
            # Pretty sure this is the subject of Cracknell & Davies 1976b (On the completeness
            # of tables of irreducible representations of the classical space groups)
        else
            kidx_BRS = findfirst(x->x==klab_BRS, klabs_BRS)
            kidx_ISO = findfirst(x->x==klab_BRS′, klabs_ISO)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO]))
            if showmissing
                if BRS.kvs[kidx_BRS] ≠ kvec(first(lgirsvec[kidx_ISO]))
                    println("Different definitions of k-point labels in space group ", sgnum)
                    println("   BRS, ", klab_BRS, ": ",string(BRS.kvs[kidx_BRS]))
                    println("   ISO, ", klabs_ISO[kidx_ISO], ": ",string(kvec(first(lgirsvec[kidx_ISO]))))
                    println()
                end
            end
            # Note that the definition of Σ appears to be different in the BRSs
            # (i.e. in Bilbao) and in ISOTROPY. Specifically, 
            #   Σ_BRS≡[α,α,0] while Σ_ISO≡[α,-2α,0]
            # This impacts space groups 146, 148, 155, 160, 161, 166, & 167
            # Ah, this is actually a parsing error on our part; it seems that 
            # we did not parse (u, -2*u, 0) correctly due to the multiplication
            # operator '*'. Need to fix that and recrawl/parse. 
            # TODO: The root cause of this problem is that KVec(::String) doesn't
            #       handle * multiply signs: as a result, we have 
            #            KVec("u,-2*u,0") = [α, α, 0.0]
            #       We could probably just strip out '*' chars early on in the method.
            #       Done! ✓
            # TODO: Need to recrawl bandreps after fixing this; and we should probably
            #       do a conversion of "GP" to "Ω" to be consistent with the notation 
            #       we used for ISOTROPY's bandreps.
        end
    end
end
end

@testset "Physically irreducible irreps/co-reps (with TR)" begin
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
for (sgnum, lgirsvec) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    irlabs_BRS = BRS.irreplabs
    klabs_BRS = BRS.klabs

    irlabs_ISO = Vector{String}()
    realirlabs_ISO = Vector{String}()
    klabs_ISO = Vector{String}(undef, length(lgirsvec))
    for (kidx, lgirs) in enumerate(lgirsvec)
        append!(irlabs_ISO,     [label(lgir) for lgir in lgirs])
        append!(realirlabs_ISO, realify(lgirs)[2])
        klabs_ISO[kidx] = klabel(first(lgirs))
    end
    irlabs_ISO = SGOps.formatirreplabel.(irlabs_ISO)
    realirlabs_ISO = SGOps.formatirreplabel.(realirlabs_ISO)

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab_BRS = klabel(irlab_BRS)
        # unfortunately, in our ISO dataset the general point is Ω; in 
        # the band rep set it is GP: put a bandaid on that for now
        klab_BRS′ = klab_BRS == "GP" ? "Ω" : klab_BRS
        irlab_BRS′ = occursin("GP", irlab_BRS) ? "Ω₁" : irlab_BRS
        if irlab_BRS′ ∉ realirlabs_ISO
            if showmissing
                @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            @test_broken false

        else
            kidx_BRS = findfirst(x->x==klab_BRS, klabs_BRS)
            kidx_ISO = findfirst(x->x==klab_BRS′, klabs_ISO)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO]))
        end
    end
end
end
end