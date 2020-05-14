using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS = parselittlegroupirreps()
end

@testset "Multiplication tables" begin

# check that operations are sorted identically across distinct irreps for fixed sgnum, and k-label
@testset "Space groups (consistent operator sorting order in ISOTROPY)" begin 
    for lgirs in LGIRS
        for lgirvec in lgirs
            ops = operations(first(lgirvec))
            for lgir in lgirvec[2:end] # don't need to compare against itself
                ops′ = operations(lgir)
                @test all(ops.==ops′)
                if !all(ops.==ops′)
                    println("sgnum = $(num(lgir)) at $(klabel(lgir))")
                    display(ops.==ops′)
                    display(ops)
                    display(ops′)
                    println('\n'^3)
                end
            end
        end
    end
end

@testset "Space groups (ISOTROPY: 3D)" begin
    numset=Int64[]
    for lgirs in LGIRS
        for lgirvec in lgirs
            sgnum = num(first(lgirvec)); cntr = centering(sgnum, 3);
            ops = operations(first(lgirvec))              # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr) # ops in primitive basis
            mt = multtable(primitive_ops)
            checkmt = @test isgroup(mt) 
            
            # for debugging
            #isgroup(mt) && union!(numset, num(first(lgirvec))) # collect info about errors, if any exist
        end
    end
    # for debugging
    #!isempty(numset) && println("The multiplication tables of $(length(numset)) little groups are faulty:\n   # = ", numset)
end

for D in 1:3
    @testset "Space groups (Bilbao: $(D)D)" begin
        for sgnum in 1:MAX_SGNUM[D]
            cntr = centering(sgnum, D);
            ops = operations(spacegroup(sgnum, Val(D)))  # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr)      # ops in primitive basis
            mt = multtable(primitive_ops)
            checkmt = @test isgroup(mt) 
        end
    end
end


@testset "Complex LGIrreps" begin
    #failcount = 0
    for lgirs in LGIRS
        for lgirvec in lgirs
            sgnum = num(first(lgirvec)); cntr = centering(sgnum, 3);
            ops = operations(first(lgirvec))        # ops in conventional basis
            primitive_ops = primitivize.(ops, cntr) # ops in primitive basis
            mt = multtable(primitive_ops)

            for lgir in lgirvec
                for αβγ = [nothing, [1,1,1]*1e-1] # test finite and zero values of αβγ
                    checkmt = checkmulttable(mt, lgir, αβγ; verbose=false)
                    @test all(checkmt)
                    #if !all(checkmt); failcount += 1; end
                end
            end
        end
    end
    #println("\nFails: $(failcount)\n\n")
end

@testset "Complex PGIrreps" begin
    for D in 1:3
        for pgiuc in PGS_IUCs[D]
            pgirs = get_pgirreps(pgiuc, D)
            pg = group(first(pgirs))
            for pgir in pgirs
                mt = multtable(operations(pg))
                @test isgroup(mt)

                checkmt = checkmulttable(mt, pgir, nothing; verbose=false)
                @test all(checkmt)
            end
        end
    end
end
end

# TODO: Test physically irreducible irreps (co-reps)