using SGOps, SmithNormalForm, Test


allpaths = false
spinful  = false
for sgnum = 1#1:230
    global BRS = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    display(BRS); 
    display(smith(matrix(BRS)).SNF)
    display(classification(BRS))
    #display(SGOps.humanreadable.(SGOps.reps(BRS)))
end

#LGIRS=parselittlegroupirreps();
println("\n\n\n\n")
# --- test complex-form irreps (not assuming time-reversal symmetry) ---
for sgnum = 1:230
    BRS = bandreps(sgnum, allpaths, spinful, "Elementary")
    irlabs_BRS = BRS.irreplabs
    klabs_BRS = BRS.klabs

    irlabs_ISOTROPY = [SGOps.formatirreplabel(label(LGIRS[sgnum][kidx][irridx])) for kidx = 1:length(LGIRS[sgnum]) for irridx in 1:length(LGIRS[sgnum][kidx])]
    klabs_ISOTROPY = [klabel(LGIRS[sgnum][kidx][1]) for kidx = 1:length(LGIRS[sgnum])]

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        if irlab_BRS ∉ irlabs_ISOTROPY
            @warn "Cannot find $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"

            kidx_BRS = findfirst(x->x==klabel(irlab_BRS), klabs_BRS)
            kidx_ISOTROPY_nearest = findfirst(x->x==string(first(klabel(irlab_BRS))), klabs_ISOTROPY)

            # test that for each of the (P,K,W,H)A k-label variants, that the associated k-vector is 
            # just equal to minus the kvector in the (P,K,W,H) variant (i.e. without the 'A' postscript)
            @test string(-BRS.kvs[kidx_BRS]) == string(kvec(LGIRS[sgnum][kidx_ISOTROPY_nearest][1]))

            # ... to get the irreps of these variants, we need to follow the prescription 
            # detailed in CDML p. 69-73 (though that won't work, presumably, for sgnum = 205)
            # Pretty sure this is the subject of Cracknell & Davies 1976b (On the completeness
            # of tables of irreducible representations of the classical space groups)
        else
            kidx_BRS = findfirst(x->x==klabel(irlab_BRS), klabs_BRS)
            kidx_ISOTROPY = findfirst(x->x==klabel(irlab_BRS), klabs_ISOTROPY)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test string(BRS.kvs[kidx_BRS]) == string(kvec(LGIRS[sgnum][kidx_ISOTROPY][1]))
        end
    end
end

println("\n\n\n\n")
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
for sgnum = 1:230
    BRS = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    irlabs_BRS = BRS.irreplabs
    klabs_BRS = BRS.klabs

    irlabs_ISOTROPY = [SGOps.formatirreplabel(label(LGIRS[sgnum][kidx][irridx])) for kidx = 1:length(LGIRS[sgnum]) for irridx in 1:length(LGIRS[sgnum][kidx])]
    klabs_ISOTROPY = [klabel(LGIRS[sgnum][kidx][1]) for kidx = 1:length(LGIRS[sgnum])]

    for (iridx_BRS, irlab_BRS) in enumerate(irlabs_BRS)
        if irlab_BRS ∉ irlabs_ISOTROPY
            @warn "Cannot find $(irlab_BRS) in ISOTROPY dataset (sgnum = $sgnum)"

            kidx_BRS = findfirst(x->x==klabel(irlab_BRS), klabs_BRS)
            kidx_ISOTROPY_nearest = findfirst(x->x==string(first(klabel(irlab_BRS))), klabs_ISOTROPY)
        else
            kidx_BRS = findfirst(x->x==klabel(irlab_BRS), klabs_BRS)
            kidx_ISOTROPY = findfirst(x->x==klabel(irlab_BRS), klabs_ISOTROPY)

            # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
            @test string(BRS.kvs[kidx_BRS]) == string(kvec(LGIRS[sgnum][kidx_ISOTROPY][1]))
        end
    end
end