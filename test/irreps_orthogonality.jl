using Crystalline, Test, LinearAlgebra

if !isdefined(Main, :LGIRS) # load complex little groups, if not already loaded
    LGIRS = get_all_lgirreps(Val(3)) # ≡ parselittlegroupirreps()
end
#Crystalline.add_ΦnotΩ_lgirs!.(LGIRS)

@testset "Irrep orthogonality (complex little groups)" begin

## 1ˢᵗ orthogonality theorem (characters): 
#       ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾
# for each irrep
# Dᵢ⁽ᵃ⁾ with i running over the Nₒₚ elements of the little group 
@testset "1ˢᵗ orthogonality theorem" begin
    for lgirs in LGIRS       # loop over space groups: lgirs contains _all_ little groups and their irreps
        for lgirvec in lgirs # loop over little group: lgirvec contains all the associated irreps
            Nₒₚ = order(first(lgirvec)) # number of elements in little group
            for lgir in lgirvec # specific irrep {Dᵢ⁽ᵃ⁾} of the little group
                χ = characters(lgir) # characters χᵢ⁽ᵃ⁾ of every operation
                @test sum(abs2, χ) ≈ Nₒₚ # check ∑ᵢ|χᵢ⁽ᵃ⁾|² = Nₒₚ⁽ᵃ⁾ 
            end
        end
    end
end

## 2ⁿᵈ orthogonality theorem (characters): 
#       ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾ = δₐᵦNₒₚ⁽ᵃ⁾  
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same little group (with 
# i running over the Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements)
@testset "2ⁿᵈ orthogonality theorem" begin
    for lgirs in LGIRS          # lgirs: vectors of little group irrep collections
        for lgirvec in lgirs    # lgirvec:  tuples of distinct little group irreps
            Nₒₚ = order(first(lgirvec))    
            for (a, lgir⁽ᵃ⁾) in enumerate(lgirvec) 
                χ⁽ᵃ⁾ = characters(lgir⁽ᵃ⁾)
                for (β, lgir⁽ᵝ⁾) in enumerate(lgirvec)
                    χ⁽ᵝ⁾ = characters(lgir⁽ᵝ⁾)
                    orthog2nd = dot(χ⁽ᵃ⁾, χ⁽ᵝ⁾) # dot conjugates the first vector automatically, 
                                                # i.e. this is just ∑ᵢχᵢ⁽ᵃ⁾*χᵢ⁽ᵝ⁾
                    @test (orthog2nd ≈ (a==β)*Nₒₚ)  atol=1e-12
                end
            end
        end
    end
end


## Great orthogonality theorem of irreps: 
#       ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ = δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
# for irreps Dᵢ⁽ᵃ⁾ and Dᵢ⁽ᵝ⁾ in the same little group (with 
# i running over the Nₒₚ = Nₒₚ⁽ᵃ⁾ = Nₒₚ⁽ᵝ⁾ elements)
@testset "Great orthogonality theorem (little group irreps)" begin
    αβγ = nothing#[1,1,1]*1e-1
    debug = false# true
    count = total = 0 # counters
    for lgirs in LGIRS          # lgirs: vectors of little group irrep collections
        for lgirvec in lgirs    # lgirvec: vector of distinct little group irreps
            Nₒₚ = order(first(lgirvec))
            for (a, lgir⁽ᵃ⁾) in enumerate(lgirvec) 
                D⁽ᵃ⁾ = irreps(lgir⁽ᵃ⁾,αβγ)      # vector of irreps in (a)
                dim⁽ᵃ⁾ = size(first(D⁽ᵃ⁾),1)

                for (β, lgir⁽ᵝ⁾) in enumerate(lgirvec)
                    D⁽ᵝ⁾ = irreps(lgir⁽ᵝ⁾,αβγ)  # vector of irreps in (β)
                    dim⁽ᵝ⁾ = size(first(D⁽ᵝ⁾),1)
                    δₐᵦ = (a==β)
                    if a == β && debug
                        display(label(lgir⁽ᵃ⁾))
                        display(label(lgir⁽ᵝ⁾))
                        g_orthog_kron = sum(kron(conj.(D⁽ᵃ⁾[i]), D⁽ᵝ⁾[i]) for i in Base.OneTo(Nₒₚ))
                        display(g_orthog_kron)
                    end
                    for n in Base.OneTo(dim⁽ᵃ⁾), j in Base.OneTo(dim⁽ᵝ⁾)     # rows of each irrep
                        δₐᵦδₙⱼ = δₐᵦ*(n==j)
                        for m in Base.OneTo(dim⁽ᵃ⁾), k in Base.OneTo(dim⁽ᵝ⁾) # cols of each irrep
                            δₐᵦδₙⱼδₘₖ = δₐᵦδₙⱼ*(m==k)

                            # compute ∑ᵢ[Dᵢ⁽ᵃ⁾]ₙₘ*[Dᵢ⁽ᵝ⁾]ⱼₖ
                            g_orthog = sum(conj(D⁽ᵃ⁾[i][n,m])*D⁽ᵝ⁾[i][j,k] for i in Base.OneTo(Nₒₚ)) 
                            
                            # test equality to δₐᵦδₙⱼδₘₖNₒₚ⁽ᵃ⁾/dim(D⁽ᵃ⁾)
                            check_g_orthog = isapprox(g_orthog, δₐᵦδₙⱼδₘₖ*Nₒₚ/dim⁽ᵃ⁾, atol=1e-12)
                            if !check_g_orthog
                                if debug
                                    println("sgnum = $(num(lgir⁽ᵃ⁾)), ",
                                            "{α = $(label(lgir⁽ᵃ⁾)), β = $(label(lgir⁽ᵝ⁾))}, ",
                                            "(n,m) = ($n,$m), (j,k) = ($j,$k)"
                                            )
                                    println("   $(g_orthog) vs. $(δₐᵦδₙⱼδₘₖ*Nₒₚ/dim⁽ᵃ⁾)")
                                end
                                count += 1
                            end
                            total += 1
                            @test check_g_orthog
                        end
                    end
                end
            end
        end
	end
	if debug
        println("\n\n", count, "/", total, " errored",
                " (", round(100*count/total, digits=1), "%)\n")
	end
end
end

if false
	println() 

	## Comparison with Bilbao's REPRES for P1 and P3 of sg 214
	D°s = Dict( # polar representation (r,ϕ) with ϕ in degrees
			:P1 => [[(1.0,  0.0)  (0.0,   0.0); (0.0,   0.0) (1.0,   0.0)],
					[(1.0, 90.0)  (0.0,   0.0); (0.0,   0.0) (1.0, 270.0)],
					[(0.0,  0.0)  (1.0, 180.0); (1.0,   0.0) (0.0,   0.0)],
					[(0.0,  0.0)  (1.0,  90.0); (1.0,  90.0) (0.0,   0.0)],
					[(1/√2,255.0) (1/√2,345.0); (1/√2, 75.0) (1/√2,345.0)],
					[(1/√2,345.0) (1/√2,255.0); (1/√2,165.0) (1/√2,255.0)],
					[(1/√2,345.0) (1/√2, 75.0); (1/√2,345.0) (1/√2,255.0)],
					[(1/√2, 75.0) (1/√2,345.0); (1/√2, 75.0) (1/√2,165.0)],
					[(1/√2,105.0) (1/√2,285.0); (1/√2, 15.0) (1/√2, 15.0)],
					[(1/√2,195.0) (1/√2,195.0); (1/√2,105.0) (1/√2,285.0)],
					[(1/√2,285.0) (1/√2,285.0); (1/√2, 15.0) (1/√2,195.0)],
					[(1/√2, 15.0) (1/√2,195.0); (1/√2,105.0) (1/√2,105.0)]],
			:P3 => [[(1.0,  0.0)  (0.0,   0.0); (0.0,   0.0) (1.0,   0.0)],
					[(1.0, 90.0)  (0.0,   0.0); (0.0,   0.0) (1.0, 270.0)],
					[(0.0,  0.0)  (1.0, 180.0); (1.0,   0.0) (0.0,   0.0)],
					[(0.0,  0.0)  (1.0,  90.0); (1.0,  90.0) (0.0,   0.0)],
					[(1/√2, 15.0) (1/√2,105.0); (1/√2,195.0) (1/√2,105.0)],
					[(1/√2,105.0) (1/√2, 15.0); (1/√2,285.0) (1/√2, 15.0)],
					[(1/√2,105.0) (1/√2,195.0); (1/√2,105.0) (1/√2, 15.0)],
					[(1/√2,195.0) (1/√2,105.0); (1/√2,195.0) (1/√2,285.0)],
					[(1/√2,345.0) (1/√2,165.0); (1/√2,255.0) (1/√2,255.0)],
					[(1/√2, 75.0) (1/√2, 75.0); (1/√2,345.0) (1/√2,165.0)],
					[(1/√2,165.0) (1/√2,165.0); (1/√2,255.0) (1/√2, 75.0)],
					[(1/√2,255.0) (1/√2, 75.0); (1/√2,345.0) (1/√2,345.0)]]
			)
	polar2complex(vs) = vs[1]*(cosd(vs[2])+1im*sind(vs[2]))
	Ds = Dict(label=>[polar2complex.(m) for m in D°] for (label, D°) in pairs(D°s))
	for (label, D) in pairs(Ds)
		Nₒₚ = length(D)
		dimD = size(first(D),1)

		kroncheck = sum(kron(conj.(D[i]), D[i]) for i in Base.OneTo(Nₒₚ))
		println("\nIrrep $(string(label)):")
		println("   |G|/l = $(Nₒₚ/dimD)")
		display(kroncheck)    
	end
end