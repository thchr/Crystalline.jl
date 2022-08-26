let subgroup_data = JLD2.load_object("data/spacegroup_subgroups_data.jld2")
    const global SUBGROUPSD_1D = Dict{Int, GroupRelations{1, SpaceGroup{1}}}()
    const global SUBGROUPSD_2D = Dict{Int, GroupRelations{2, SpaceGroup{2}}}()
    const global SUBGROUPSD_3D = Dict{Int, GroupRelations{3, SpaceGroup{3}}}()
    for (D, subgroupsd) in zip(1:3, (SUBGROUPSD_1D, SUBGROUPSD_2D, SUBGROUPSD_3D))
        for (numᴳ, data) in enumerate(subgroup_data[D])
            subgroups = Vector{GroupRelation{D,SpaceGroup{D}}}(undef, length(data))
            for (i, (numᴴ, index, kindstr, Pps)) in enumerate(data)
                transforms = Vector{ConjugacyTransform{D}}(undef, length(Pps))
                for (j, (P, p)) in enumerate(Pps)
                    P′ = isnothing(P) ? P : SqSMatrix{D,Float64}(P)
                    p′ = isnothing(p) ? p : SVector{D,Float64}(p)
                    transforms[j] = ConjugacyTransform(P′, p′)
                end

                kind = kindstr == :t ? TRANSLATIONENGLEICHE : KLASSENGLEICHE
                subgroups[i] = GroupRelation{SpaceGroup{D}}(numᴳ, numᴴ, index, kind, transforms)
            end
            subgroupsd[numᴳ] = GroupRelations{SpaceGroup{D}}(numᴳ, reverse!(subgroups))
        end
    end
end