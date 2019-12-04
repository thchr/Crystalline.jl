function read_pgops_xyzt(label::String, D::Integer=2)
    if D ∉ (2,3); throw(DomainError(D, "dimension D must be 2 or 3")); end

    filepath = (@__DIR__)*"/../data/pgops/"*string(D)*"d/"*label*".json"
    symops_str = open(filepath) do io
        JSON2.read(io)
    end
    return symops_str
end

function get_pgops(label::String, D::Integer=2)
    D ≠ 2 && throw(DomainError("Only 2D point groups are implemented at this point"))
    pgnum = pointgroup_label2num(label)
    sgops_str = read_pgops_xyzt(label, D)
    
    return PointGroup(pgnum, label, SymOperation.(sgops_str), D)
end

function get_pgops(pgnum::Integer, D::Integer=2)
    D ≠ 2 && throw(DomainError("Only 2D point groups are implemented at this point"))
    label = pointgroup_num2label(pgnum)
    sgops_str = read_pgops_xyzt(label, D)

    return PointGroup(pgnum, label, SymOperation.(sgops_str), D)
end

# TODO: These labels can be mapped to Schoenflies notation, but the mapping is not one-to-one
#       but rather many-to-one; e.g. 3m1 and 31m maps to C₃ᵥ but correspond to different axis
#       choices. 
#       When we get to the 3D point groups this is even more important, as the axis systems 
#       then factor in as well: e.g. unique axis a vs. unique axis c and rhombohedral vs. hexagonal
const POINTGROUP_NUM2LABEL = ("1", "2", "m", "mm2", "4", "4mm", "3", "3m1", "31m", "6", "6mm")
pointgroup_num2label(num::Int64) = POINTGROUP_NUM2LABEL[num]
function pointgroup_label2num(label::String)
    # could use Unrolled.jl to make this as fast as a direct lookup 
    # (not sure if it can unroll through `findfirst`)
    @inbounds for (idx, label′) in pairs(POINTGROUP_NUM2LABEL)
        label′ == label && return idx
    end
    return nothing
end

function find_parent_pointgroup(ops::AbstractVector{SymOperation})
    D = dim(first(ops))
    D ≠ 2 && throw(DomainError("Only 2D point groups are implemented at this point"))
    xyzt_pgops = sort(xyzt.(pointgroup(ops)))

    idx = findfirst(pgnum->sort(xyzt.(operations(get_pgops(pgnum, D))))==(xyzt_pgops),  eachindex(POINTGROUP_NUM2LABEL))
    if idx !== nothing
        return POINTGROUP_NUM2LABEL[idx]
    else
        throw("Didn't find a matching point group! Ensure that you are using a canonical axis system.")
    end
end
find_parent_pointgroup(g::AbstractGroup) = find_parent_pointgroup(operations(g))