"""
    pointgroup(iuclab::String, ::Union{Val{D}, Integer}=Val(3),
               ::Union{Val{S}, Bool}=Val(false))  -->  PointGroup{D} or DPointGroup{D}

Return the symmetry operations associated with the point group identified with label
`iuclab` in dimension `D` as a `PointGroup{D}`.

If `S` is `true`, the double group is returned instead, as a `DPointGroup{D}` (currently
supported in 3D only).
"""
function pointgroup(
    iuclab::AbstractString,
    Dᵛ::Val{D}=Val(3),
    spinfulᵛ::Val{S}=Val(false)
) where {D, S}
    @boundscheck _check_valid_pointgroup_label(iuclab, D)
    pgnum = pointgroup_iuc2num(iuclab, D) # this is not generally a particularly well-established numbering
    return _pointgroup(iuclab, pgnum, Dᵛ, spinfulᵛ)
end
@inline function pointgroup(iuclab::String, D::Integer, spinful::Bool=false)
    return pointgroup(iuclab, Val(D), Val(spinful))
end

"""
    pointgroup(pgnum::Integer, ::Union{Val{D}, Integer}=Val(3), setting::Integer=1,
               ::Union{Val{S}, Bool}=Val(false))  -->  PointGroup{D} or DPointGroup{D}

Return the symmetry operations associated with the point group identfied with canonical
number `pgnum` in dimension `D` as a `PointGroup{D}`. The connection between a point group's
numbering and its IUC label is enumerated in `Crystalline.PG_NUM2IUC[D]` and
`Crystalline.IUC2NUM[D]`.

Certain point groups feature in multiple setting variants: e.g., IUC labels 321 and 312 both
correspond to `pgnum = 18` and correspond to the same group structure expressed in two
different settings. The `setting` argument allows choosing between these setting variations.

If `S` is `true`, the double group is returned instead, as a `DPointGroup{D}` (currently
supported in 3D only).
"""
function pointgroup(
    pgnum::Integer,
    Dᵛ::Val{D}=Val(3),
    setting::Integer=1,
    spinfulᵛ::Val{S}=Val(false)
) where {D, S}
    iuclab = pointgroup_num2iuc(pgnum, Dᵛ, setting) # also checks validity of `(pgnum, D)`
    return _pointgroup(iuclab, pgnum, Dᵛ, spinfulᵛ)
end
@inline function pointgroup(pgnum::Integer, D::Integer, setting::Integer=1,
                            spinful::Bool=false)
    return pointgroup(pgnum, Val(D), setting, Val(spinful))
end

function _pointgroup(iuclab::String, pgnum::Integer, Dᵛ::Val{D}, ::Val{S}) where {D, S}
    S && D ≠ 3 && _only_3d(D)
    codes = PG_CODES_Ds[D][iuclab]

    Nop = (length(codes)+1) # number of operations
    operations = Vector{SymOperation{D}}(undef, Nop)
    operations[1] = one(SymOperation{D})
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(get_indexed_rotation(code, Dᵛ), 
                             zero(SVector{D,Float64}))
        operations[n+1] = op
    end

    pg = PointGroup{D}(pgnum, iuclab, operations)
    return S ? doublegroup(pg) : pg
end

function _check_valid_pointgroup_label(iuclab::AbstractString, D)
    D ∉ (1,2,3) && _throw_invalid_dim(D)
    if iuclab ∉ PG_IUCs[D]
        throw(DomainError(iuclab, 
                "iuc label not found in database (see possible labels in `PG_IUCs[D]`)"))
    end
end