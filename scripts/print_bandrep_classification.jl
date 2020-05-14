allpaths = false
spinful  = false
for sgnum = 1:230
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
    display(BRS); 
    #display(smith(matrix(BRS)).SNF)
    display(classification(BRS))
    #display(Crystalline.humanreadable.(Crystalline.reps(BRS)))
end