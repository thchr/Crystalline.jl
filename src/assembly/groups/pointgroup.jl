"""
    pointgroup(iuclab::String, ::Union{Val{D}, Integer}=Val(3))  -->  PointGroup{D}

Return the symmetry operations associated with the point group identified with label
`iuclab` in dimension `D` as a `PointGroup{D}`.
"""
function pointgroup(iuclab::AbstractString, Dᵛ::Val{D}=Val(3)) where D
    @boundscheck _check_valid_pointgroup_label(iuclab, D)
    pgnum = pointgroup_iuc2num(iuclab, D) # this is not generally a particularly well-established numbering
    return _pointgroup(iuclab, pgnum, Dᵛ)
end
@inline pointgroup(iuclab::String, D::Integer) = pointgroup(iuclab, Val(D))

"""
    pointgroup(pgnum::Integer, ::Union{Val{D}, Integer}=Val(3), setting::Integer=1)
                                                                      -->  PointGroup{D}

Return the symmetry operations associated with the point group identfied with canonical
number `pgnum` in dimension `D` as a `PointGroup{D}`. The connection between a point group's
numbering and its IUC label is enumerated in `Crystalline.PG_NUM2IUC[D]` and
`Crystalline.IUC2NUM[D]`.

Certain point groups feature in multiple setting variants: e.g., IUC labels 321 and 312 both
correspond to `pgnum = 18` and correspond to the same group structure expressed in two
different settings. The `setting` argument allows choosing between these setting variations.
"""
function pointgroup(pgnum::Integer, Dᵛ::Val{D}=Val(3), setting::Integer=1) where D
    iuclab = pointgroup_num2iuc(pgnum, Dᵛ, setting) # also checks validity of `(pgnum, D)`
    return _pointgroup(iuclab, pgnum, Dᵛ)
end
@inline function pointgroup(pgnum::Integer, D::Integer, setting::Integer=1)
    return pointgroup(pgnum, Val(D), setting)
end

function _pointgroup(iuclab::String, pgnum::Integer, Dᵛ::Val{D}) where D
    codes = PG_CODES_Ds[D][iuclab]

    Nop = (length(codes)+1) # number of operations
    operations = Vector{SymOperation{D}}(undef, Nop)
    operations[1] = one(SymOperation{D})
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(get_indexed_rotation(code, Dᵛ), 
                             zero(SVector{D,Float64}))
        operations[n+1] = op
    end

    return PointGroup{D}(pgnum, iuclab, operations)
end

function _check_valid_pointgroup_label(iuclab::AbstractString, D)
    D ∉ (1,2,3) && _throw_invalid_dim(D)
    if iuclab ∉ PG_IUCs[D]
        throw(DomainError(iuclab, 
                "iuc label not found in database (see possible labels in `PG_IUCs[D]`)"))
    end
end