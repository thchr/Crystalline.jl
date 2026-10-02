function generators(
    iuclab::String,
    ::Type{T}
) where {D, T <: AbstractPointGroup{D}}
    S = isspinful(T)
    S && (D == 3 || _only_3d(D))
    @boundscheck _check_valid_pointgroup_label(iuclab, D)
    codes = PG_GENS_CODES_Ds[D][iuclab]
    # `hexagonal` is only relevant for spinful
    hexagonal = S ? _ishexagonal_pg(pointgroup_iuc2num(iuclab, D), Val(D)) : false

    # convert `codes` to `SymOperation`s and add to `operations`
    operations = Vector{eltype(T)}(undef, length(codes) + S) # extra slot for Ē if spinful
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(get_indexed_rotation(code, Val{D}()), zero(SVector{D,Float64}))
        operations[n] = _maybe_attach_su2(eltype(T), op, hexagonal)
    end
    S && (operations[end] = DSymOperation{D}(one(SymOperation{D}), -one(SU2))) # add on Ē

    return operations
end
# NB: The method below exists separately, rather than putting `::Type{T} = PointGroup{3}` in
#     the signature above, because Julia otherwise warns about an unused `D` parameter in
#     the method signature, due to the automatically generated `generators(iuclab)` method.
generators(iuclab::String) = generators(iuclab, PointGroup{3})

# ---------------------------------------------------------------------------------------- #

function generators(
    pgnum::Integer,
    ::Type{T},
    setting::Integer=1
) where {D, T <: AbstractPointGroup{D}}
    iuclab = pointgroup_num2iuc(pgnum, Val(D), setting)
    return generators(iuclab, T)
end