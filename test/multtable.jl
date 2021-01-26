using Crystalline, Test

if !isdefined(Main, :LGIRSDIM)
    LGIRSDIM = Tuple(get_lgirreps.(1:MAX_SGNUM[D], Val(D)) for D in 1:3)
end

@testset "Multiplication tables" begin

# check that operations are sorted identically across distinct irreps for fixed sgnum, and k-label
@testset "Space groups (consistent operator sorting order in ISOTROPY)" begin
for LGIRS in LGIRSDIM # ... D in 1:3
    for lgirsd in LGIRS
        for lgirs in values(lgirsd)
            ops = operations(first(lgirs))
            for lgir in lgirs[2:end] # don't need to compare against itself
                ops′ = operations(lgir)
                @test all(ops.==ops′)
                # if !all(ops.==ops′)
                #     println("sgnum = $(num(lgir)) at $(klabel(lgir))")
                #     display(ops.==ops′)
                #     display(ops)
                #     display(ops′)
                #     println('\n'^3)
                # end
            end
        end
    end
end # for LGIRS in LGIRSDIM
end # @testset "Space groups (consistent ..."

@testset "Little groups: group property" begin
#numset=Int64[]
for (D, LGIRS) in enumerate(LGIRSDIM)
    for lgirsd in LGIRS
        for lgirs in values(lgirsd)
            sgnum = num(first(lgirs))
            cntr = centering(sgnum, D)
            ops = operations(first(lgirs))          # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr) # ops in primitive basis
            mt = MultTable(primitive_ops)
            @test mt.isgroup
            
            # for debugging
            #mt.isgroup && union!(numset, num(first(lgirs))) # collect info about errors, if any exist
        end
    end
end # for LGIRS in LGIRSDIM
# for debugging
#!isempty(numset) && println("The multiplication tables of $(length(numset)) little groups are faulty:\n   # = ", numset)
end # @testset "Little groups: ..."

for D in 1:3
    @testset "Space groups (Bilbao: $(D)D)" begin
        for sgnum in 1:MAX_SGNUM[D]
            cntr = centering(sgnum, D)
            ops = operations(spacegroup(sgnum, Val(D)))  # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr)      # ops in primitive basis
            mt = MultTable(primitive_ops)
            @test mt.isgroup
        end
    end
end


@testset "Complex LGIrreps" begin
#failcount = 0
for (D, LGIRS) in enumerate(LGIRSDIM)
    for lgirsd in LGIRS
        for lgirs in values(lgirsd)
            sgnum = num(first(lgirs))
            cntr = centering(sgnum, D)
            ops = operations(first(lgirs))          # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr) # ops in primitive basis
            mt = MultTable(primitive_ops)

            for lgir in lgirs
                for αβγ in (nothing, Crystalline.TEST_αβγs[D]) # test finite and zero values of αβγ
                    checkmt = Crystalline.check_multtable_vs_ir(mt, lgir, αβγ; verbose=false)
                    @test all(checkmt)
                    #if !all(checkmt); failcount += 1; end
                end
            end
        end
    end
end # for LGIRS in LGIRSDIM
#println("\nFails: $(failcount)\n\n")
end # @testset "Complex LGIrreps"

@testset "Complex PGIrreps" begin
    for D in 1:3
        for pgiuc in PGS_IUCs[D]
            pgirs = get_pgirreps(pgiuc, D)
            pg = group(first(pgirs))
            for pgir in pgirs
                mt = MultTable(operations(pg))
                @test mt.isgroup

                checkmt = Crystalline.check_multtable_vs_ir(mt, pgir, nothing; verbose=false)
                @test all(checkmt)
            end
        end
    end
end
end

# TODO: Test physically irreducible irreps (co-reps)