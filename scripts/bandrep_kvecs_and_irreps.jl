using SGOps, Test

if !isdefined(Main, :LGIRS)
    LGIRS = get_all_lgirreps(3)  # loaded from our saved .jld2 files
end

allpaths = false
spinful = false
showmissing = true
@testset begin "Copy of test in test/bandrep.jl"
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
incomplete_sgs = Vector{Int64}()
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
        if irlab_BRS ∉ realirlabs_ISO
            if showmissing
                @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            sgnum ∉ incomplete_sgs && push!(incomplete_sgs, sgnum)
            @test_broken false
            
        else
            kidx_BRS = findfirst(==(klab_BRS), klabs_BRS)
            kidx_ISO = findfirst(==(klab_BRS), klabs_ISO)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test BRS.kvs[kidx_BRS] == kvec(first(lgirsvec[kidx_ISO]))
        end
    end
end
if !isempty(incomplete_sgs)
    print(stdout, "The following space groups have incomplete little group (TR-invariant) irreps in ISOTROPY:\n\t")
    join(stdout, incomplete_sgs, ", "); println()
end
end