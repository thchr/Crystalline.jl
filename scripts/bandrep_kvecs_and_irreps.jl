using SGOps, Test

if !isdefined(Main, :LGIRS)
    LGIRS = get_all_lgirreps(3)  # loaded from our saved .jld2 files
end

allpaths = true
spinful  = false
brtype   = "Elementary TR"
showmissing = true
# --- test physically irreducible irreps/co-reps ---
@testset "Are all the k-vectors from bandreps in ISOTROPY?" begin
incomplete_sgs = Vector{Int64}()
incomplete_lgirs = Dict{Int64,Array{String}}()
for (sgnum, lgirsvec) in enumerate(LGIRS)
    BRS = bandreps(sgnum, allpaths, spinful, brtype)
    irlabs_BRS = BRS.irreplabs
    klabs_BRS = BRS.klabs

    irlabs_ISO = Vector{String}()
    klabs_ISO = Vector{String}(undef, length(lgirsvec))
    for (kidx, lgirs) in enumerate(lgirsvec)
        if brtype == "Elementary"
            append!(irlabs_ISO, [label(lgir) for lgir in lgirs])
        elseif brtype == "Elementary TR"
            append!(irlabs_ISO, realify(lgirs)[2])
        end
        klabs_ISO[kidx] = klabel(first(lgirs))
    end
    irlabs_ISO = SGOps.formatirreplabel.(irlabs_ISO)


    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        klab_BRS = klabel(irlab_BRS)
        if irlab_BRS ∉ irlabs_ISO
            if showmissing
                @info "Cannot find real irrep $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"
            end
            sgnum ∉ incomplete_sgs && push!(incomplete_sgs, sgnum)
            if !haskey(incomplete_lgirs, sgnum)
                incomplete_lgirs[sgnum] = [irlab_BRS]
            else
                push!(incomplete_lgirs[sgnum], irlab_BRS)
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

if !isempty(incomplete_sgs)
    print(stdout, "The following space groups have incomplete little group (TR-invariant) irreps in ISOTROPY:\n\t")
    join(stdout, incomplete_sgs, ", "); println()
    sorted_incomplete_lgirs = sort(keys(incomplete_lgirs) .=> values(incomplete_lgirs), by=x->getfield(x,1))
    display(sorted_incomplete_lgirs)
end
end