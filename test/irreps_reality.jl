using Crystalline, Test #, Crayons

if !isdefined(Main, :LGIRS)
    LGIRS = parselittlegroupirreps()
end

# Compare our calculation of the Herring criterion with the tabulated 
# reality types in ISOTROPY
@testset "Herring criterion" begin
    #= error_count = 0 =#       # for debugging
    for (sgnum, lgirsd) in enumerate(LGIRS)
        sgops = operations(first(lgirsd["Γ"])) # this is important: we may _not_ use trivial repeated sets, i.e. spacegroup(..) would not work generally
        for lgirs in values(lgirsd)
            for lgir in lgirs
                iso_rawtype = type(lgir)
                iso_type = iso_rawtype == 1 ? 1 : (iso_rawtype == 2 ? -1 : 0) # map ISOTROPY's types from {1,2,3} to {1,-1,0} = {real, pseudoreal, complex}
                #= try =#       # for debugging
                herring_type = herring(lgir, sgops)
                @test iso_type ≈ herring_type
                #= catch err    # for debugging
                    if true
                        println(sgnum, ", ", 
                                Crystalline.subscriptify.(label(lgir)), ", ", 
                                centering(sgnum,3), ", ",
                                issymmorph(sgnum))
                    end
                    error_count += 1
                    println(err)
                end =#
            end
        end
    end

    #=                          # for debugging
    if error_count>0
        println(Crayon(foreground=:red, bold=true), "Outright errors: ", error_count)
    else
        println(Crayon(foreground=:green, bold=true), "No errors")
    end
    =#
end