using Pkg: activate, project
dirname(project().path) == (@__DIR__) || activate(@__DIR__)

using Crystalline

# "lift" a 2D operation into three dimensions
f(x) = SymOperation(hcat(vcat(hcat(rotation(x), [0; 0]), [0 0 1]),
                         vcat(translation(x), [0])))

println("Identical sorting and setting of point group operators in 3D vs. 2D?")
for pgiuc in PGS_IUCs[2]
    ops3d = operations(pointgroup(pgiuc, Val(3)))
    ops2d = operations(pointgroup(pgiuc, Val(2)))
    ops2d′ = f.(ops2d)

    same_setting = all(xyzt.(ops3d) .== xyzt.(ops2d′))
    println("   ", pgiuc, ":\t", same_setting)
    if !same_setting
        println.(Ref("      "), seitz.(ops3d), Ref(" vs. "), seitz.(ops2d′))
    end
end

# the above results imply that we can employ a whole-sale transfer of irreps from 3D point
# groups to 2D point groups, because the 2D. vs. 3D settings never differ in a material way
# (only "2" and "m" point groups have such differences; but there, it doesn't impact the 
# transfer, because there's only two operations in each of these groups)