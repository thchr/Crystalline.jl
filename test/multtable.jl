using SGOps, Test

if !isdefined(Main, :LGIR)
    LGIR = parselittlegroupirreps.()
end
@testset "Operations, matched sorting" begin # check that operations are sorted identically across distinct irreps for fixed sgnum, and k-label
    for lgirs in LGIR
        for lgirvec in lgirs
            ops = operations(first(lgirvec))
            for lgir in lgirvec
                ops′ = operations(lgir)
                @test all(ops.==ops′)
            end
        end
    end
end

@testset "Multiplication table, groups" begin
    numset=Int64[]
    for lgirs in LGIR
        for lgirvec in lgirs
            ops = operations(first(lgirvec));
            mt = multtable(ops)
            #@test isgroup(mt)
            
            if !isgroup(mt)
                #println(num(first(lgirvec)), " ",
                #        centering(num(first(lgirvec))), " ", 
                #        klabel(first(lgirvec)), " " )
                union!(numset, num(first(lgirvec)))
            end
        end
    end
    @test_broken false
    println("The multiplication tables of __$(length(numset))__ little groups are faulty:\n   ", numset)
end

@testset "Multiplication table, irreps" begin
    for lgirs in LGIR
        for lgirvec in lgirs
            ops = operations(first(lgirvec));
            mt = multtable(ops)
            if isgroup(mt)
                for lgir in lgirvec
                    for αβγ = [nothing, [1,1,1]*1e-1] # test finite and zero values of αβγ
                        checkmt = checkmulttable(mt, lgir, αβγ)
                        @test all(checkmt)
                    end
                end
            else
                @test_broken false # doesn't currently work for groups with trivial fractional "translation sets"
            end 
        end
    end
end

