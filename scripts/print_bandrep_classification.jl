using Pkg: activate, project
dirname(project().path) == (@__DIR__) || activate(@__DIR__)

using Crystalline

allpaths = false
spinful  = false
timereversal = true
for sgnum = 1:230
    BRS = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=timereversal)
    display(BRS)
    display(classification(BRS))
end