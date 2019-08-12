using SGOps


allpaths = false
spinful  = false
for sgnum = 212#1:230
    global BRS = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    display(BRS); 
    display(matrix(BRS))
    #display(SGOps.humanreadable.(SGOps.reps(BRS)))
end

