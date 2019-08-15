using SGOps


allpaths = false
spinful  = false
for sgnum = 1#1:230
    global BRS = bandreps(sgnum, allpaths, spinful, "Elementary TR")
    display(BRS); 
    display(smith(matrix(BRS)).SNF)
    display(classification(BRS))
    #display(SGOps.humanreadable.(SGOps.reps(BRS)))
end

