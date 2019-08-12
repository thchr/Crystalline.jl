using SGOps

BRS = Vector{BandRepSet}(undef, 230)
allpaths = false
spinful  = true
for sgnum = 230#1:230
    BRS[sgnum] = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    display(BRS[sgnum]); println()
end

