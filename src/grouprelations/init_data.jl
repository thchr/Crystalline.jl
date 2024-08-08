# ---------------------------------------------------------------------------------------- #
# Dictionaries of maximal subgroup relationships in 1D, 2D, & 3D
const MAXSUB_RELATIONSD_1D = Dict{Int, GroupRelations{1, SpaceGroup{1}}}()
const MAXSUB_RELATIONSD_2D = Dict{Int, GroupRelations{2, SpaceGroup{2}}}()
const MAXSUB_RELATIONSD_3D = Dict{Int, GroupRelations{3, SpaceGroup{3}}}()
let subgroup_data = JLD2.load_object(joinpath((@__DIR__), "..", "..", "data", 
                                              "spacegroup_subgroups_data.jld2"))
    for (D, maxsubsd) in 
            zip(1:3, (MAXSUB_RELATIONSD_1D, MAXSUB_RELATIONSD_2D, MAXSUB_RELATIONSD_3D))
        for (numᴳ, data) in enumerate(subgroup_data[D])
            grouprels = Vector{GroupRelation{D,SpaceGroup{D}}}(undef, length(data))
            for (i, (numᴴ, index, kindstr, Pps)) in enumerate(data)
                classes = Vector{ConjugacyTransform{D}}(undef, length(Pps))
                for (j, (P, p)) in enumerate(Pps)
                    P′ = isnothing(P) ? P : SqSMatrix{D,Float64}(P)
                    p′ = isnothing(p) ? p : SVector{D,Float64}(p)
                    classes[j] = ConjugacyTransform(P′, p′)
                end

                kind = kindstr == :t ? TRANSLATIONENGLEICHE : KLASSENGLEICHE
                grouprels[i] = GroupRelation{SpaceGroup{D}}(numᴳ, numᴴ, index, kind, classes)
            end
            maxsubsd[numᴳ] = GroupRelations{SpaceGroup{D}}(numᴳ, reverse!(grouprels))
        end
    end
end

# ---------------------------------------------------------------------------------------- #
# Dictionaries of minimal supergroup relationships in 1D, 2D, & 3D
const MINSUP_RELATIONSD_1D = Dict{Int, GroupRelations{1, SpaceGroup{1}}}()
const MINSUP_RELATIONSD_2D = Dict{Int, GroupRelations{2, SpaceGroup{2}}}()
const MINSUP_RELATIONSD_3D = Dict{Int, GroupRelations{3, SpaceGroup{3}}}()
for (D, minsupsd, maxsubsd) in 
        zip(1:3,
            (MINSUP_RELATIONSD_1D, MINSUP_RELATIONSD_2D, MINSUP_RELATIONSD_3D),
            (MAXSUB_RELATIONSD_1D, MAXSUB_RELATIONSD_2D, MAXSUB_RELATIONSD_3D))
    for numᴴ in 1:MAX_SGNUM[D]
        grouprels = Vector{GroupRelation{D,SpaceGroup{D}}}()
        for numᴳ in 1:MAX_SGNUM[D]
            # NB: the loop-range above cannot be reduced to `numᴴ:MAX_SGNUM[D]` since
            # that would exclude some klassengleiche relationships where `numᴴ > numᴳ`
            rsᴳᴴ = maxsubsd[numᴳ]
            idxs = findall(rᴳᴴ -> rᴳᴴ.num==numᴴ, rsᴳᴴ)
            for i in idxs
                rᴳᴴ = rsᴳᴴ[i]
                rᴴᴳ = GroupRelation{SpaceGroup{D}}(numᴴ, numᴳ,
                        rᴳᴴ.index, rᴳᴴ.kind, 
                        rᴳᴴ.classes #= NB: note, per Bilbao convention, we do not 
                                        invert transformations in `.classes` =#)
                push!(grouprels, rᴴᴳ)
            end
        end
        minsupsd[numᴴ] = GroupRelations{SpaceGroup{D}}(numᴴ, grouprels)
    end
end
