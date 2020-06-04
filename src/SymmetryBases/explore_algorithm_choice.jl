using Crystalline, PyCall
import Base: OneTo
if !isdefined(Main, :PyNormaliz)
    const PyNormaliz = pyimport("PyNormaliz")
end

algos = ("DualMode", "PrimalMode")
cumtimes = zeros(2)
times = zeros(2)
for sgnum in 1:230
    global cumtimes, times, algos
    print("SG", sgnum, ": ")
    BRS = bandreps(sgnum, spinful=false, timereversal=true)
    
    B = matrix(BRS, true) # Matrix with columns of EBRs.

    F = Crystalline.smith(B)    
    dᵇˢ = count(!iszero, F.SNF)

    if sgnum ∈ (47, 123, 191, 131, 221)
        println("   ... skipping - too time consuming (many conditions: ", dᵇˢ,")!")
        continue
    end

    S = F.S[:,OneTo(dᵇˢ)]                 # These are all the nontrivial conditions on zⱼ 
    C = PyNormaliz.Cone(inequalities = S)
    default_time = @elapsed C.HilbertBasis()

    for (i, mode) in enumerate(algos)
        C′ = PyNormaliz.Cone(inequalities = S)
        times[i] = @elapsed C′.Compute("HilbertBasis", mode)
    end
    cumtimes = cumtimes .+ times

    _, match = findmin(abs.(times .- default_time))
    notmatch = match == 1 ? 2 : 1

    # Choice expected based on         
    r = C.Rank()# monoid rank
    e = C.EmbeddingDim() # embedding dimension
    s = size(B, 1) # number of constraints = number of irreps
    chooses_dualmode = r + 50/r ≤ s ≤ 2e # see Normaliz manual, sec. 6.5

    # Print results/info
    println("Default:     ", algos[match], "\t",
            "[predicted ", chooses_dualmode ? "DualMode" : "PrimalMode", "]")
    print(" "^(4+ndigits(sgnum)), "Non-default: ", algos[notmatch], "\t")
    if times[notmatch] < times[match] 
        println("Speedup: ×", round(times[match]/times[notmatch], digits=1))
    else
        println("Slowdown: ×", round(times[notmatch]/times[match], digits=1))
    end
    println()
end

# Cumulative timings
println("Cumulative timings:")
println.(Ref("   "), algos, Ref(":\t"), round.(cumtimes, digits=1), Ref(" s"));