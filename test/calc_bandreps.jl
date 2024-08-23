using Test
using Crystalline
using Crystalline: dlm2struct

@testset "calc_bandreps" begin

# defines `is_exceptional_br` to check if a BR induced by a maximal Wyckoff position is an
# "exceptional" BR, i.e., is in fact not an elementary BR (EBR)
include("ebr_exceptions.jl")

# ---------------------------------------------------------------------------------------- #

@testset "2D: `calc_bandreps` vs. (tabulated) `bandreps`" begin
    # test that `bandreps` agrees with `calc_bandreps` in 2D (since the former is generated
    # by the latter; basically check that the data behind the former is up-to-date)
    for sgnum in 1:17
        for timereversal in (false, true)
            @test calc_bandreps(sgnum, Val(2); timereversal, allpaths = false) ==
                  bandreps(sgnum, 2; timereversal, allpaths = false)
        end
    end
end

@testset "3D: every reference EBR must have a match in a calculated BR" begin
    debug = false
    error_counts = Dict{String, Int}()
    for sgnum in 1:230
        had_sg_error = false
        for timereversal in (false, true)
            had_tr_error = false

            brsᶜ = calc_bandreps(sgnum, Val(3); timereversal, allpaths = false)
            brsʳ = bandreps(sgnum, 3; timereversal, allpaths = false)
                      
            # find a permutation of the irreps in `brsʳ` that matches `brsᶜ`'s sorting, s.t.
            # `brsʳ.irlabs[irʳ²ᶜ_perm] == brsᶜ.irlabs`
            irʳ²ᶜ_perm = [something(findfirst(==(irᶜ), brsʳ.irlabs)) for irᶜ in brsᶜ.irlabs]
            append!(irʳ²ᶜ_perm, length(brsᶜ.irlabs)+1) # append occupation number

            seen_wp = Dict{String, Bool}()
            for brʳ in brsʳ
                wpʳ = brʳ.wyckpos
                idx = findfirst(brᶜ -> brᶜ.wyckpos == wpʳ && brᶜ.label == brʳ.label, brsᶜ)
                haskey(seen_wp, wpʳ) || (seen_wp[wpʳ] = false)
                if isnothing(idx)
                    # this can sometimes happen spuriously because Bilbao (i.e. `brsʳ`)
                    # lists the EBRs with unabbreviated Mulliken labels (e.g., A₁ rather 
                    # than A or ¹E²E rather than E) even when it is possible to abbreviate
                    # unambiguously (what `mulliken` does and what `brsᶜ` consequently
                    # references); we account for that by doing a subsequent, slighly looser
                    # check (below)
                    idx = findfirst(brsᶜ) do brᶜ
                        brʳ.wyckpos == brᶜ.wyckpos || return false
                        labʳ = replace(brʳ.label, "↑G"=>"") 
                        labᶜ = replace(brᶜ.label, "↑G"=>"")
                        length(labᶜ) == 1 && only(labᶜ) == first(labʳ) && return true # e.g., A ~ A₁
                        labʳ′ = replace(labʳ, "¹"=>"", "²"=>"")    # e.g., ¹E²E ~ E
                        labʳ′ = replace(labʳ′, "EE" => "E", 
                                               "EgEg" => "Eg", "EᵤEᵤ" => "Eᵤ",
                                               "E₁E₁" => "E₁", "E₂E₂" => "E₂",
                                               "E′′E′′" => "E′′", "E′E′" => "E′",
                                               "E₁gE₁g" => "E₁g", "E₂gE₂g" => "E₂g",
                                               "E₁ᵤE₁ᵤ" => "E₁ᵤ", "E₂ᵤE₂ᵤ" => "E₂ᵤ")
                        labʳ′ == labᶜ && return true
                    end
                end
                
                @test idx !== nothing # test that an associated BR exists in brsʳ (it must)
                idx = something(idx)

                # there are a number of cases where it is impossible to test uniquely
                # against Bilbao, because the assignment of irrep label to the site symmetry
                # group is ambiguous (e.g., for site groups isomorphic to 222; but also for 
                # mmm and mm2): in this case, the assignment of irreps depends on a choice
                # of coordinate system, which we cannot be uniquely determined but
                # so, for these cases, we cannot do a one-to-one comparison (but we can at
                # least test whether the reference Bilbao vector occurs in the set of
                # computed band representation vectors)
                if brʳ.sitesym ∉ ("222", "mmm", "mm2", "-4m2")
                    brᶜ = brsᶜ[idx]
                    @test Vector(brᶜ) == brʳ[irʳ²ᶜ_perm]
                    @test dim(brᶜ) == dim(brʳ)
                elseif brʳ.sitesym == "-4m2"
                    # this is a special case: in principle, it is possible to uniquely
                    # assign the irrep labels, but we currently assign ones that are
                    # different from those in Bilbao (see issue #59)
                    # FIXME: remove this here and above once #59 is resolved
                    @test brʳ[irʳ²ᶜ_perm] ∈ Vector.(brsᶜ)
                    continue
                else
                    # simply test that the reference Bilbao _vector_ is in the set of
                    # computed EBRs
                    @test brʳ[irʳ²ᶜ_perm] ∈ Vector.(brsᶜ)
                    continue
                end
                
                if debug
                    brᶜ = brsᶜ[idx]
                    if !(Vector(brᶜ) == brʳ[irʳ²ᶜ_perm])
                        had_sg_error || (print("sgnum = $sgnum\n"); had_sg_error = true)
                        had_tr_error || (println("  timereversal = $timereversal"); had_tr_error = true)
                        if !seen_wp[wpʳ]
                            print("    site = $wpʳ (")
                            printstyled("$(brʳ.sitesym)"; color=:yellow)
                            println(")")
                            seen_wp[wpʳ] = true
                        end
                        println("      $(replace(brʳ.label, "↑G"=>""))")
                        brʳ′ = deepcopy(brʳ)
                        brʳ′.irvec .= brʳ.irvec[irʳ²ᶜ_perm[1:end-1]]
                        brʳ′.irlabs .= brʳ.irlabs[irʳ²ᶜ_perm[1:end-1]]
                        println("        ref  = $brʳ′")
                        println("        calc = $brᶜ")
                        if haskey(error_counts, brʳ.sitesym)
                            error_counts[brʳ.sitesym] += 1
                        else
                            error_counts[brʳ.sitesym] = 1
                        end
                    end
                end
            end
        end
        debug && had_sg_error && println()
    end
end

@testset "Plane groups vs. parent space groups" begin
    for (sgnum²ᴰ, sgnum³ᴰ) in enumerate(Crystalline.PLANE2SPACE_NUMS)
        for timereversal in (false, true)
            brs²ᴰ = bandreps(sgnum²ᴰ, 2; timereversal, allpaths = false)
            brs³ᴰ = bandreps(sgnum³ᴰ, 3; timereversal, allpaths = false)

            # dimensions
            @test sort!(dim.(brs²ᴰ)) == sort!(dim.(brs³ᴰ))

            # topological classification
            @test classification(brs²ᴰ) == classification(brs³ᴰ)
        end
    end
end
end # @testset "calc_bandreps" begin