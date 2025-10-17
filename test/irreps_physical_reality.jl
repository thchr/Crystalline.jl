using LinearAlgebra: eigvals
using Crystalline
using Test

function _eigsortby(λ::Number)
    r, i = reim(round(λ, digits=13)) # round slightly to avoid sorting issues
    # normalize -0.0 to 0.0 for similar reasons (don't want (-0.0, 1) to sort differently
    # than (0.0, 1))
    return ((iszero(r) ? zero(r) : r), (iszero(i) ? zero(i) : i))
end
@testset "physical_realify" begin
    @testset "pgirreps" begin
        for D in 1:3
            for pglab in Crystalline.PG_IUCs[D]
                irs = pgirreps(pglab, D)
                re_irs  = realify(irs)
                re_irs′ = physical_realify(irs)

                @test characters(re_irs′) ≈ characters(re_irs)
                @test all(ir -> all(x -> x ≈ real(x), ir.matrices), re_irs′) # D ≈ Re(D)
                @test all(zip(re_irs, re_irs′)) do (re_ir, re_ir′) # remains similar
                    all(zip(re_ir.matrices, re_ir′.matrices)) do (X, Y)
                        # check similarity by checking identical eigenvalues: this is a bit
                        # annoying to do robustly with complex floating point numbers,
                        # because two close complex numbers might be sorted into different
                        # positions due to rounding errors; so we hack it by rounding down
                        # to fewer digits, and then compare after re-sorting
                        sort!(eigvals(X), by=_eigsortby) ≈ sort!(eigvals(Y), by=_eigsortby)
                    end
                end

                # invariant to before/after realification
                re_irs′′ = physical_realify(re_irs)
                @test all(zip(re_irs′, re_irs′′)) do (re_ir′, re_ir′′)
                    re_ir′.matrices == re_ir′′.matrices
                end
                @test group(first(re_irs′)) == group(first(re_irs′′))
                @test label.(re_irs′) == label.(re_irs′′)
            end
        end
    end # @testset "pgirreps"

    @testset "siteirreps" begin
        for D in 1:3
            for sgnum in MAX_SGNUM[D]
                wps = wyckoffs(sgnum, D)
                sg = spacegroup(sgnum, D)
                for wp in wps
                    siteg = sitegroup(sg, wp)
                    irs = siteirreps(siteg)
                    re_irs  = realify(irs)
                    re_irs′ = physical_realify(irs)

                    @test characters(re_irs′) ≈ characters(re_irs)
                    @test all(ir -> all(x -> x ≈ real(x), ir.matrices), re_irs′) # D ≈ Re(D)
                    @test all(zip(re_irs, re_irs′)) do (re_ir, re_ir′) # remains similar
                        all(zip(re_ir.matrices, re_ir′.matrices)) do (X, Y)
                            sort!(eigvals(X), by=_eigsortby) ≈ sort!(eigvals(Y), by=_eigsortby)
                        end
                    end
        
                    # invariant to before/after realification
                    re_irs′′ = physical_realify(re_irs)
                    @test all(zip(re_irs′, re_irs′′)) do (re_ir′, re_ir′′)
                        re_ir′.matrices == re_ir′′.matrices
                    end
                    @test group(first(re_irs′)) == group(first(re_irs′′))
                    @test label.(re_irs′) == label.(re_irs′′)
                end
            end
        end
    end # @testset "siteirreps"
end # @testset "physical_realify"
