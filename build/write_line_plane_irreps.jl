# Write the little group irreps of the plane (2D) and line (1D) groups, derived from their
# parent space groups (see `LinePlaneIrreps.jl`).
#
# Currently, only the double-valued irreps are written. The single-valued ones are still
# those of `setup_2d_littlegroup_irreps.jl` (2D) and `write_littlegroup_irreps.jl` (1D),
# whose labels in places differ from the parents' (see `compare_line_plane_irreps.jl`); to
# replace them, uncomment the second loop.
#
#   julia --project=. build/write_line_plane_irreps.jl

include(joinpath(@__DIR__, "line_plane_irreps.jl"))
using .LinePlaneIrreps

for D in (2, 1)
    write_lgirreps(Val(D); spinful=Val(true))
end
# for D in (2, 1)
#     write_lgirreps(Val(D); spinful=Val(false))
# end
