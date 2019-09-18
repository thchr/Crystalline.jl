using SGOps, Test, Crayons

#LGIRS=parselittlegroupirreps();

# Compare our calculation of the Herring criterion with the tabulated 
# reality types in ISOTROPY
@testset "Herring criterion" begin
    #= error_count = 0 =#
    for sgnum = 1:230
        sgops = operations(LGIRS[sgnum][1][1]) # this is important: we may _not_ use trivial repeated sets, i.e. get_symops would not work generally
        for kidx = 1:length(LGIRS[sgnum])
            for iridx = 1:length(LGIRS[sgnum][kidx])
                iso_rawtype = type(LGIRS[sgnum][kidx][iridx])
                iso_type = iso_rawtype == 1 ? 1 : (iso_rawtype == 2 ? -1 : 0) # map ISOTROPY's types from {1,2,3} to {1,-1,0} = {real, pseudoreal, complex}
                #= try =#
                herring_sum, herring_norm = herring(LGIRS[sgnum][kidx][iridx], sgops)
                herring_type = herring_sum ≠ 0 ? herring_sum/herring_norm : 0
                @test iso_type ≈ herring_type
                #= catch err # debugging printouts
                    if true
                        println(sgnum, ", ", 
                                SGOps.subscriptify.(label(LGIRS[sgnum][kidx][iridx])), ", ", 
                                centering(sgnum,3), ", ",
                                issymmorph(sgnum))
                    end
                    error_count += 1
                    println(err)
                end =#
            end
        end
    end

    #= 
    if error_count>0
        println(Crayon(foreground=:red, bold=true), "Outright errors: ", error_count)
    else
        println(Crayon(foreground=:green, bold=true), "No errors")
    end
    =#
end
    

linelen = 60
for sgnum = 1:230
    println("\n\nSpace group $(sgnum)\n", "──┬", '─'^linelen)
    for kidx = 1:length(LGIRS[sgnum])
        realify(LGIRS[sgnum][kidx], true)
        if kidx ≠ length(LGIRS[sgnum]) 
            println("──┼", '─'^linelen)
        end
    end
    println("──┴", '─'^linelen)
end