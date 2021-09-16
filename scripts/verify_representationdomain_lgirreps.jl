using Crystalline, LinearAlgebra, Test

# load `LGIRS_add` with some special-point lgirs
include(joinpath((@__DIR__), "../data/stopgap_missing_lgirreps.jl"))

# Choose one of 195, 198, 200, 201
#for sgnum in [195, 198, 200, 201] # works for 195 and 200; not for 198 and 201
#sgnum = 201
#target_klab = "ZA"
#sgnum = 24
#target_klab = "WA"
sgnum = 82
target_klab = "PA"

lgirsd = get_lgirreps(sgnum, Val(3))
haskey(lgirsd, target_klab) && throw("`get_lgirreps` already returns this k-label")
target_lgirs = LGIRS_add[sgnum][target_klab]

added_lgirs = add_ΦnotΩ_lgirs!(deepcopy(lgirsd), true)[target_klab];

#=
println("\nOriginal: ", string(target_lgirs[1].lg.kv))
display.(target_lgirs)
println("\nComputed: ", string(added_lgirs[1].lg.kv))
display(added_lgirs[1].lg.kv)
display.(added_lgirs)
println()
=#
##

# Check operator sorting and k-vector
#@test target_lgirs[1].lg.kv == added_lgirs[1].lg.kv
@test all(operations(group(first(target_lgirs))) .== operations(group(first(added_lgirs))))

# Print some info
println('\n','-'^25, ' ', sgnum, ' ', '-'^25, "\nk = ", string(kvec(group(first(target_lgirs)))), '\n')
print("ops = ")
join(stdout, seitz.(operations(group(first(target_lgirs)))), ", ")
println('\n')

# Difference
αβγ = [0.3,0.2,0.4]

δχ = characters.(target_lgirs, Ref(αβγ)) .- characters.(added_lgirs, Ref(αβγ))
δP = [lgir(αβγ) for lgir in target_lgirs] .- [lgir(αβγ) for lgir in added_lgirs]

println.(label.(added_lgirs), Ref(": |δχ| = "), norm.(δχ), ", |δP| = ", norm.(δP)); println()

τs = getfield.(added_lgirs, :translations)
println.(label.(added_lgirs), Ref(": τs = "), τs)
println('-'^55)
#end