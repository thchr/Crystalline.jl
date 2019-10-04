using SGOps, Test

if !isdefined(Main, :LGIRS)
    LGIRS = parselittlegroupirreps()
end

# check that operations are sorted identically across distinct irreps for fixed sgnum, and k-label
@testset "Operations, matching sorting (ISOTROPY)" begin 
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

@testset "Multiplication table, groups (ISOTROPY: 3D)" begin
    numset=Int64[]
    for lgirs in LGIRS
        for lgirvec in lgirs
            sgnum = num(first(lgirvec)); cntr = centering(sgnum, 3);
            ops = operations(first(lgirvec))              # ops in conventional basis
            primitive_ops = SGOps.primitivize.(ops, cntr) # ops in primitive basis
            mt = multtable(primitive_ops)
            checkmt = @test isgroup(mt) 
            
            # for debugging
            #isgroup(mt) && union!(numset, num(first(lgirvec))) # collect info about errors, if any exist
        end
    end
    # for debugging
    #!isempty(numset) && println("The multiplication tables of $(length(numset)) little groups are faulty:\n   # = ", numset)
end


sgnums = (1:17, 1:230)
for dim in 2:3
    @testset "Multiplication table, groups (Bilbao: $(dim)D)" begin
        for sgnum in sgnums[dim-1]
            cntr = centering(sgnum, dim);
            ops = operations(get_sgops(sgnum, dim))      # ops in conventional basis
            primitive_ops = SGOps.primitivize.(ops, cntr) # ops in primitive basis
            mt = multtable(primitive_ops)
            checkmt = @test isgroup(mt) 
        end
    end
end

@testset "Multiplication table, irreps" begin
    for lgirs in LGIRS
        for lgirvec in lgirs
            sgnum = num(first(lgirvec)); cntr = centering(sgnum, 3);
            ops = operations(first(lgirvec))              # ops in conventional basis
            primitive_ops = SGOps.primitivize.(ops, cntr) # ops in primitive basis
            mt = multtable(primitive_ops)

            for lgir in lgirvec
                for αβγ = [nothing, [1,1,1]*1e-1] # test finite and zero values of αβγ
                    checkmt = checkmulttable(mt, lgir, αβγ)
                    @test all(checkmt)
                end
            end
        end
    end
end

# TODO: Test physically irreducible irreps (co-reps)