allpaths = false
spinful  = false
for sgnum = 1:230
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
    display(BRS)
    display(classification(BRS))
end