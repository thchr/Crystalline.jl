using Crystalline, Test

@testset "Every time-reversal invariant momentum is listed" begin
    # the time-reversal invariant momenta (TRIMs) in 1D, 2D, and 3D, in a primitive basis
    function get_trims(::Val{D}) where D
        if D == 1
            ((0.0,), (0.5,))
        elseif D == 2
            ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5))
        elseif D == 3
            ((0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 0.5),
             (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.5, 0.5, 0.5))
        else
            error("invalid dimension $D")
        end
    end

    # with time-reversal, a TRIM is special even if its little group is trivial (e.g.,
    # hosting Kramers pairs for spinful particles): so every TRIM must be in the star of a
    # listed k-point, also in the hand-made 2D and 1D listings
    function covers_all_trims(lgs::Dict{String, LittleGroup{D}}, cntr::Char) where D
        covered = Set{NTuple{D, Float64}}()
        for lg in values(lgs)
            isspecial(lg) || continue
            for kv in orbit(lg)
                k = primitivize(kv, cntr)()
                if all(x -> isapprox(2x, round(2x); atol=1e-8), k)
                    push!(covered, NTuple{D, Float64}(mod.(round.(2k) ./ 2, 1)))
                end
            end
        end
        return all(∈(covered), get_trims(Val(D)))
    end
    # the irreps, single- and double-valued, are listed at the k-points of `littlegroups`
    for D in 1:3, sgnum in 1:MAX_SGNUM[D]
        lgs = littlegroups(sgnum, Val(D))
        @test covers_all_trims(lgs, centering(sgnum, D))
    end

    # the TRIMs with a trivial little group in 2D, which were initially missing
    @test issubset(["Y", "B", "A"], keys(littlegroups(1, Val(2))))  # p1
    @test "S" ∈ keys(littlegroups(5, Val(2)))                       # c1m1
    @test "M" ∈ keys(littlegroups(13, Val(2)))                      # p3
end
