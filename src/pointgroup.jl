function read_pgops_xyzt(iuclab::String, D::Integer=2)
    if D ∉ (1,2,3); throw(DomainError(D, "dimension D must be 1, 2, or 3")); end
    
    filepath = (@__DIR__)*"/../data/pgops/"*string(D)*"d/"*replace(iuclab, "/"=>"_slash_")*".json"
    symops_str = open(filepath) do io
        JSON2.read(io)
    end
    return symops_str
end

function get_pgops(iuclab::String, D::Integer=2)
    if D ∉ (1,2,3); throw(DomainError(D, "dimension D must be 1, 2, or 3")); end
    pgnum = pointgroup_iuc2num(iuclab, D)  # this is not generally a particularly meaningful or even well-established numbering
    sgops_str = read_pgops_xyzt(iuclab, D)
    
    return PointGroup(pgnum, iuclab, SymOperation.(sgops_str), D)
end

function get_pgops(pgnum::Integer, D::Integer=3)
    if D ∉ (1,2,3); throw(DomainError(D, "dimension D must be 1, 2, or 3")); end
    iuclab = pointgroup_num2iuc(pgnum, D)
    sgops_str = read_pgops_xyzt(iuclab, D)

    return PointGroup(pgnum, iuclab, SymOperation.(sgops_str), D)
end

# We include several axes settings; as a result, there are more than 32 point groups 
# under the 3D case, because some variations are "setting-degenerate" (but are needed
# to properly match all space group settings)
const NUM2IUC_PGS = (
    ("1", "m"),                                                         # 1D
    ("1", "2", "m", "mm2", "4", "4mm", "3", "3m1", "31m", "6", "6mm"),  # 2D 
    ("1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m", # 3D
     "422", "4mm", 
     "-42m", "-4m2",                             # D2d setting variations
     "4/mmm", "3", "-3", 
     "312", "321", "3m1", "31m", "-31m", "-3m1", # hexagonal lattices w/ setting variations
     "6", "-6", "6/m", "622", "6mm", 
     "-62m", "-6m2",                             # D3d setting variations
     "6/mmm", "23", "m-3", "432", 
     "-43m", "m-3m")
)

pointgroup_num2iuc(pgnum::Integer, D::Integer) = NUM2IUC_PGS[D][pgnum]
function pointgroup_iuc2num(iuclab::String, D::Integer)
    # could use Unrolled.jl to make this as fast as a direct lookup 
    # (not sure if it can unroll through `findfirst`)
    @inbounds for (idx, iuclab′) in enumerate(NUM2IUC_PGS[D])
        iuclab′ == iuclab && return idx
    end
    return nothing
end

pointgroup_iuc2schoenflies(iuclab::String) = IUC2SCHOENFLIES_PGS[iuclab]

# The IUC notation for point groups can be mapped to the Schoenflies notation, but the 
# mapping is not one-to-one but rather one-to-many; e.g. 3m1 and 31m maps to C3c but 
# correspond to different axis orientations. 
# When there is a choice of either hexagonal vs. rhombohedral or unique axes b vs unique
# axes a/c we choose hexagonal and unique axes b, respectively.
const IUC2SCHOENFLIES_PGS = ImmutableDict(
    "1"     => "C1",   "-1"    => "Ci",
    "2"     => "C2",   "m"     => "Cs",   "2/m"   => "C2h",  # unique axes b setting
    "222"   => "D2",   "mm2"   => "C2v",  "mmm"   => "D2h",  "4"    => "C4",
    "-4"    => "S4",   "4/m"   => "C4h",  "422"   => "D4",   "4mm"  => "C4v", 
    "-42m"  => "D2d",  "-4m2"  => "D2d",  # D2d setting variations
    "4/mmm" => "D4h",  "3"     => "C3",   "-3"    => "C3i",  
    "312"   => "D3",   "321"   => "D3",   # D3 setting variations  (hexagonal axes)
    "3m1"   => "C3v",  "31m"   => "C3v",  # C3v setting variations (hexagonal axes)
    "-31m"  => "D3d",  "-3m1"  => "D3d",  # D3d setting variations (hexagonal axes)
    "6"     => "C6",   "-6"    => "C3h",  "6/m"   => "C6h",  "622"  => "D6", 
    "6mm"   => "C6v",  
    "-62m"  => "D3h", "-6m2"   => "D3h",  # D3h setting variations
    "6/mmm" => "D6h",  "23"    => "T",
    "m-3"   => "Th",   "432"   => "O",    "-43m"  => "Td",   "m-3m" => "Oh"
)

function find_parent_pointgroup(g::AbstractGroup)
    D = dim(g)
    xyzt_pgops = sort(xyzt.(pointgroup(g)))

    @inbounds for iuclab in NUM2IUC_PGS[D]
        P = get_pgops(iuclab, D)
        if sort(xyzt.(operations(P)))==(xyzt_pgops)
            return P
        end
    end
    return nothing

end