using SGOps

function writebandreps(sgnum, allpaths, brtype="Elementary TR")
    paths_str = allpaths ? "allpaths" : "maxpaths"

    BR_dlm = SGOps.html2dlm(crawlbandreps(sgnum, allpaths, brtype), "⊕", SGOps.BandRepTrait())

    filename = (@__DIR__)*"/../data/bandreps/3d/$(filter(!isspace, brtype))/$(paths_str)/$(string(sgnum)).csv"
    open(filename; write=true, create=true, truncate=true) do io
        write(io, BR_dlm)
    end 
end

# run to crawl everything...
for allpaths = false# [true, false]
    for brtype = ["Elementary"]#["Elementary TR", "Elementary"]
        for sgnum in 1:230
            display(sgnum)
            writebandreps(sgnum, allpaths, brtype)
        end
    end
end