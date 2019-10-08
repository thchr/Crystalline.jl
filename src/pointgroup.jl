function read_pgops_xyzt(label::String, dim::Integer=2)
    if all(dim .!= [2,3]); throw(DomainError(dim, "dim must be 2 or 3")); end

    filepath = (@__DIR__)*"/../data/pgops/"*string(dim)*"d/"*label*".json"
    symops_str = open(filepath) do io
        JSON2.read(io)
    end
    return symops_str
end

function get_pgops(label::String, dim::Integer=2)
    dim ≠ 2 && throw(DomainError("Only two-dimensional point groups are implemented at this point"))
    pgnum = pointgroup_label2num(label)
    sgops_str = read_pgops_xyzt(label, dim)
    
    return PointGroup(pgnum, label, SymOperation.(sgops_str), dim)
end

function get_pgops(pgnum::Integer, dim::Integer=2)
    dim ≠ 2 && throw(DomainError("Only two-dimensional point groups are implemented at this point"))
    label = pointgroup_num2label(pgnum)
    sgops_str = read_pgops_xyzt(label, dim)

    return PointGroup(pgnum, label, SymOperation.(sgops_str), dim)
end

const POINTGROUP_NUM2LABEL = ("1", "2", "m", "mm2", "4", "4mm", "3", "3m1", "31m", "6", "6mm")
pointgroup_num2label(num::Int64) = POINTGROUP_NUM2LABEL[num]
function pointgroup_label2num(label::String)
    if     label === "1";   return 1
    elseif label === "2";   return 2
    elseif label === "m";   return 3
    elseif label === "mm2"; return 4
    elseif label === "4";   return 5
    elseif label === "4mm"; return 6
    elseif label === "3";   return 7
    elseif label === "3m1"; return 8
    elseif label === "31m"; return 9
    elseif label === "6";   return 10
    elseif label === "6mm"; return 11
    else;                   return nothing
    end
end