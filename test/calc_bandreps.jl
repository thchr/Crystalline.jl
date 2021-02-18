using Test
using Crystalline
using Crystalline: dlm2struct

# ---------------------------------------------------------------------------------------- #
# load .csv-style data from .jl file
include("calc_bandreps_csvdata.jl") # defines `reference_csv::Dict`

# create `BandRepSet`s from `reference_csv` and `calc_bandreps`
key_T = NamedTuple{(:num, :tr, :allpaths),Tuple{Int64,Bool,Bool}}
reference_brs = Dict{key_T,BandRepSet}()
calc_brs = Dict{key_T,BandRepSet}()
for num in (1, 16, 17)
    for tr in (false, true)
        for allpaths in (false, true)
            num == 1 && !(!tr && allpaths)  && continue # csv not stored

            key = (num = num, tr = tr, allpaths = allpaths)
            io = IOBuffer(reference_csv[key])
            reference_brs[key] = dlm2struct(io, num, allpaths, #=spinful=# false, tr)
            calc_brs[key] = calc_bandreps(num, Val(2), allpaths=allpaths, timereversal=tr)
        end
    end
end

# ---------------------------------------------------------------------------------------- #

@testset "2D: calc_bandreps vs. bandreps" begin
    # test that `bandreps` agrees with `calc_bandreps` in 2D (since the former is generated
    # by the latter; basically check that the data behind the former is up-to-date)
    for sgnum in 1:17
        for tr in (false, true)
            for allpaths in (false, true)
                @test calc_bandreps(sgnum, Val(2), timereversal=tr, allpaths=allpaths) ==
                      bandreps(sgnum, 2, timereversal=tr, allpaths=allpaths)
            end
        end
    end
end

# TODO: Increase maximum 3D sgnum to 230 and add list of exceptions for "exceptional" band
#       representations induced by maximal wyckoff positions
@testset "3D: Checking Wyckoff Position Sets" begin
    for sgnum in 1:87
        for tr in (false, true)
            for allpaths in (true)
                brsᶜ = calc_bandreps(sgnum, Val(3), timereversal=tr, allpaths=allpaths)
                brsʳ = bandreps(sgnum, timereversal=tr, allpaths=allpaths)

                # wyckoff labels
                @test unique!(sort!(wyck.(brsᶜ))) == unique!(sort!(wyck.(brsʳ)))

                # dimensions
                #println(sgnum, " | tr = ", tr, " | allpaths = ", allpaths)
                @test sort!(dim.(brsᶜ)) == sort!(dim.(brsʳ))  # fails

                if allpaths # too many differences between ISOTROPY's and Bilbao's inclusion
                            # "non-special" k-points; don't bother comparing
                    @test sort(irreplabels(brsᶜ)) ⊆ sort(irreplabels(brsʳ))
                end
            end
        end
    end
end

@testset "2D: Checking dims against known bandreps for allpaths=true, timereversal=false" begin
    for sgnum in (1, 16, 17)
        for tr in (false, true)
            for allpaths in (false, true)
                sgnum == 1 && !(!tr && allpaths) && continue # csv not stored
                
                key = (num = sgnum, tr = tr, allpaths = allpaths)
                brsʳ = reference_brs[key] # reference version
                brsᶜ = calc_brs[key]          # calculated version
                # check dimensions             
                @test sort!(dim.(brsʳ)) == sort!(dim.(brsᶜ))

                # check wyckoff labels
                @test sort!(wyck.(brsʳ)) == sort!(wyck.(brsᶜ))

                # check irrep labels
                if !(sgnum == 17 && allpaths) # there's a T point in 2D that Bilbao didn't include in 3D; so we skip
                    # only test subset-equality, since there's unfortunately some mismatch
                    # between the included irreps/k-points
                    @test sort(irreplabels(brsᶜ)) ⊆ sort(irreplabels(brsʳ))
                end
            end
        end
    end
end

irvec(br::BandRep) = br.irvec
@testset "2D: Checking irvecs against verified csv files" begin
    for sgnum in (16, 17)
        for tr in (false, true)
            key  = (num=sgnum, tr=tr, allpaths=false)
            brsʳ = reference_brs[key]
            brsᶜ = calc_brs[key]
            
            # TODO: This is too loose a check most likely. Could do better by searching &
            #       matching klabels + wyckoff-labels & site-symmetry labels
            for irvecᶜ in irvec.(brsᶜ)
                found_irvec = false
                for irvecʳ in irvec.(brsʳ)
                    if Set(irvecʳ) == Set(irvecᶜ)
                        found_irvec = true
                        break
                    end
                end
                @test found_irvec
            end
        end
    end
end


@testset "2D: Checking irvecs through irlab permutation" begin
    for sgnum in (16, 17)
        for tr in (false, true)
            key  = (num = sgnum, tr = tr, allpaths = false)
            
            brsʳ    = reference_brs[key]
            irlabsʳ = irreplabels(brsʳ)
            irvecsʳ = irvec.(brsʳ)

            brsᶜ    = calc_brs[key]
            irlabsᶜ = irreplabels(brsᶜ)
            irvecsᶜ = irvec.(brsᶜ)

            # permuted indices of the irlabels (indices of the "reference" version for each
            # calculated version's index in order)
            irlabs_permᶜ²ʳ = Int[]
            for irlabᶜ in irlabsᶜ
                append!(irlabs_permᶜ²ʳ, findfirst(==(irlabᶜ), irlabsʳ))
            end
            println(irlabsʳ[irlabs_permᶜ²ʳ])
            println(irlabsᶜ)
            
            println(irvecsᶜ[1][irlabs_permᶜ²ʳ])
            for irvecʳ in irvecsʳ
                @test findfirst(==(irvecʳ[irlabs_permᶜ²ʳ]), irvecsᶜ) != nothing
            end
        end
    end
end


@testset "3D: Checking irvecs in 3D" begin
    for sgnum in 1:87
        for tr in (false, true)
            brsᶜ = calc_bandreps(sgnum, Val(3), allpaths = false, timereversal = tr)
            brsʳ = bandreps(sgnum, allpaths = false, timereversal = tr)

            # TODO: Too loose; see comments in analogous 2D check
            for irvecᶜ in irvec.(brsᶜ)
                found_irvec = false
                for irvecʳ in irvec.(brsʳ)
                    if Set(irvecʳ) == Set(irvecᶜ)
                        found_irvec = true
                        break
                    end
                end
                @test found_irvec
            end
        end
    end
end

@testset "Plane groups vs. parent space groups" begin
    parent³ᴰ_nums = [1, 3, 6, 7, 8, 25, 28, 32, 35, 75, 99, 100, 143, 156, 157, 168, 183]
    for (sgnum²ᴰ, sgnum³ᴰ) in zip(1:17, parent³ᴰ_nums)
        for tr in (false, true)
            for allpaths in (false, true)
                brs²ᴰ = bandreps(sgnum²ᴰ, 2, timereversal=tr, allpaths=allpaths)
                brs³ᴰ = bandreps(sgnum³ᴰ, 3, timereversal=tr, allpaths=allpaths)

                # dimensions
                @test sort!(dim.(brs²ᴰ)) == sort!(dim.(brs³ᴰ))

                # topological classification
                @test classification(brs²ᴰ) == classification(brs³ᴰ)
            end
        end
    end
end