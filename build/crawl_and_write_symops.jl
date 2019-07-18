using JSON2, SGOps

#= Small convenience script to crawl and subsequently write all the short-
   hand forms of the symmetry operations of the 230 three-dimensional
   space-groups. Enables us to just read the symmetry data from the hard-
   disk rather than constantly querying the Bilbao server =#
for i = 1:230
    symops_str = SGOps.crawl_symops_shorthand(i)
    filename = (@__DIR__)*"/../json/symops/3d/sg"*string(i)*".json"
    open(filename; write=true, create=true, truncate=true) do io
        JSON2.write(io, symops_str)
    end 
end

