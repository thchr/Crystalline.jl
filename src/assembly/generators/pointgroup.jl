function generators(iuclab::String, ::Type{PointGroup{D}}=PointGroup{3}) where D
    @boundscheck _check_valid_pointgroup_label(iuclab, D)
    codes = PG_GENS_CODES_Ds[D][iuclab]

    # convert `codes` to `SymOperation`s and add to `operations`
    operations = Vector{SymOperation{D}}(undef, length(codes))
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(get_indexed_rotation(code, Val{D}()), zero(SVector{D,Float64}))
        operations[n] = op
    end

    return operations
end
function generators(pgnum::Integer, ::Type{PointGroup{D}}, setting::Integer=1) where D
    iuclab = pointgroup_num2iuc(pgnum, Val(D), setting)
    return generators(iuclab, PointGroup{D})
end