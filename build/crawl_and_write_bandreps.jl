using SGOps

function writebandreps(sgnum, allpaths, timereversal=true)
    paths_str = allpaths ? "allpaths" : "maxpaths"
    brtype_str = timereversal ? "ElementaryTR" : "Elementary"

    BR_dlm = SGOps.html2dlm(crawlbandreps(sgnum, allpaths, timereversal), '⊕', SGOps.BandRepTrait())

    filename = (@__DIR__)*"/../data/bandreps/3d/$(brtype_str)/$(paths_str)/$(string(sgnum)).csv"
    open(filename; write=true, create=true, truncate=true) do io
        write(io, BR_dlm)
    end 
end

# run to crawl everything... (takes about 10-15 min)
for allpaths in [true, false]
    for timereversal in [true, false]
        for sgnum in 1:230
            display(sgnum)
            writebandreps(sgnum, allpaths, timereversal)
        end
    end
end