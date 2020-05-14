using Crystalline, Test, PrettyTables

debug = false
#LGIRS = parselittlegroupirreps()
@testset "Order of space group, Bilbao vs. ISOTROPY, in primitivized basis" begin
    for sgnum = 1:230
        cntr = centering(sgnum, 3)
        # sgops from Bilbao (BCD)
        ops = operations(spacegroup(sgnum))
        reduce_ops
        # go to primitive basis ⇒ reduce to unique set ⇒ go back to conventional basis
        ops′ = reduce_ops(ops, cntr, true)
        if debug
            if length(ops) ≠ length(ops′) # a trivial translation set was then removed
                println(sgnum, ": ", length(ops), " ≠ ", length(ops′))
            end
        end
        
        Nops_ISO = length(operations(LGIRS[sgnum][1][1]))
        @test Nops_ISO == length(ops′) # test that ISOTROPY spacegroup(..) indeed excludes trivial translation sets
    end
end


let count = 0, failures = Int[]
    for sgnum = 1:230
        cntr = centering(sgnum, 3)
        # sgops from Bilbao (BCD)
        ops_BCD  = operations(spacegroup(sgnum))
        ops_BCD′ = reduce_ops(ops_BCD, cntr, true) # go to primitive basis ⇒ reduce to unique 
                                                   # set ⇒ go back to conventional basis

        # sgops from ISOTROPY (via Γ point)
        ops_ISO = operations(LGIRS[sgnum][1][1])
        
        # sorting according to seitz notation
        sorted_idx_BCD = sortperm(seitz.(ops_BCD′))
        sorted_idx_ISO = sortperm(seitz.(ops_ISO))
        ops_BCD′ = ops_BCD′[sorted_idx_BCD]
        ops_ISO = ops_ISO[sorted_idx_ISO]

        # extracting various (sorted) metrics of the sgops
        seitz_BCD  = seitz.(ops_BCD′)
        seitz_ISO  = seitz.(ops_ISO)
        matrix_BCD = matrix.(ops_BCD′)
        matrix_ISO = matrix.(ops_ISO)
        τ_BCD      = translation.(ops_BCD′)
        τ_ISO      = translation.(ops_ISO)

        # comparisons of (sorted) sgops across BCD and ISO
        BCD_vs_ISO = seitz_BCD .== seitz_ISO                          # 6 disagreements (trivial differences of **primitive** lattice translations)

        # print some stuff if BCD and ISO sgops sets not equivalent
        if any(!, BCD_vs_ISO) 
            count += 1
            push!(failures, sgnum)
            dτ = τ_BCD .- τ_ISO
            P = Crystalline.primitivebasismatrix(cntr,3)
            dτ_primitive = [P\dτ_i for dτ_i in dτ]
            if debug
                println("\nsgnum=$(sgnum) ($(bravaistype(sgnum, 3))):")
                pretty_table([seitz_BCD seitz_ISO BCD_vs_ISO dτ dτ_primitive], 
                            ["BCD", "ISOTROPY", "==?", "dτ (conventional basis)", "dτ (primitive basis)"],
                            highlighters=Highlighter((d,i,j)->!d[i,3],
                                                    Crayon(background =:red)))
            end
        end
    end
    print("\n$(count)/230 disagreements, for space groups: ")
    join(stdout, failures, ", "); println()
end