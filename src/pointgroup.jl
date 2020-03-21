# ===== CONSTANTS =====

# We include several axes settings; as a result, there are more than 32 point groups 
# under the 3D case, because some variations are "setting-degenerate" (but are needed
# to properly match all space group settings)
const PGS_NUM2IUC = (
    (["1"], ["m"]),                                           # 1D
    (["1"], ["2"], ["m"], ["mm2"], ["4"], ["4mm"], ["3"],     # 2D
     ["3m1", "31m"],       # C3v setting variations
     ["6"], ["6mm"]),
    (["1"], ["-1"], ["2"], ["m"], ["2/m"], ["222"], ["mm2"],  # 3D
     ["mmm"], ["4"], ["-4"], ["4/m"], ["422"], ["4mm"], 
     ["-42m", "-4m2"],     # D2d setting variations
     ["4/mmm"], ["3"], ["-3"],
     ["312", "321"],       # D3 setting variations  (hexagonal axes)
     ["3m1", "31m"],       # C3v setting variations (hexagonal axes)
     ["-31m", "-3m1"],     # D3d setting variations (hexagonal axes)
     ["6"], ["-6"], ["6/m"], ["622"], ["6mm"], 
     ["-62m", "-6m2"],     # D3h setting variations
     ["6/mmm"], ["23"], ["m-3"], ["432"], ["-43m"], ["m-3m"])
)
# a flat tuple-listing of all the iuc labels in PGS_NUM2IUC; sliced across dimensions
const PGS_IUCs = map(x->tuple(Iterators.flatten(x)...), PGS_NUM2IUC)
# a tuple of ImmutableDicts, giving maps from iuc label to point group number
const PGS_IUC2NUM = tuple([ImmutableDict([lab=>findfirst(x->lab∈x, PGS_NUM2IUC[D])
                           for lab in PGS_IUCs[D]]...) for D = Base.OneTo(3)]...)
# The IUC notation for point groups can be mapped to the Schoenflies notation, but the 
# mapping is not one-to-one but rather one-to-many; e.g. 3m1 and 31m maps to C3v but 
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


# ===== METHODS =====

# --- Notation ---
function pointgroup_iuc2num(iuclab::String, D::Integer)
    return get(PGS_IUC2NUM[D], iuclab, nothing)
end

pointgroup_iuc2schoenflies(iuclab::String) = IUC2SCHOENFLIES_PGS[iuclab]

# --- Point groups & operators ---
unmangle_pgiuclab(iuclab) = replace(iuclab, "/"=>"_slash_")

function read_pgops_xyzt(iuclab::String, ::Val{D}=Val(3)) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)
    iuclab ∉ PGS_IUCs[D] && throw(DomainError(iuclab, "iuc label not found in database"))

    filepath = (@__DIR__)*"/../data/pgops/"*string(D)*"d/"*unmangle_pgiuclab(iuclab)*".json"
    ops_str = open(filepath) do io
        JSON2.read(io)
    end
    return ops_str
end
read_pgops_xyzt(iuclab::String, D::Integer) = read_pgops_xyzt(iuclab, Val(D))

@inline function pointgroup(iuclab::String, Dᵛ::Val{D}=Val(3)) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)
    pgnum = pointgroup_iuc2num(iuclab, D) # this is not generally a particularly well-established numbering
    ops_str = read_pgops_xyzt(iuclab, Dᵛ)
    
    return PointGroup{D}(pgnum, iuclab, SymOperation{D}.(ops_str))
end
@inline pointgroup(iuclab::String, D::Integer) = pointgroup(iuclab, Val(D))

@inline function pointgroup_num2iuc(pgnum::Integer, Dᵛ::Val{D}, setting::Integer) where D
    iucs = PGS_NUM2IUC[D][pgnum]
    length(iucs) < setting && throw(DomainError(setting, "invalid setting request"))
    return iucs[setting]
end
@inline function pointgroup(pgnum::Integer, Dᵛ::Val{D}=Val(3), setting::Int=1) where D
    D ∉ (1,2,3) && _throw_invaliddim(D)
    iuclab = pointgroup_num2iuc(pgnum, Dᵛ, setting)
    ops_str = read_pgops_xyzt(iuclab, Dᵛ)

    return PointGroup{D}(pgnum, iuclab, SymOperation{D}.(ops_str))
end
@inline pointgroup(pgnum::Integer, D::Integer, setting::Integer=1) = pointgroup(pgnum, Val(D), setting)

# --- POINT GROUPS VS SPACE & LITTLE GROUPS ---
function find_parent_pointgroup(g::AbstractGroup)
    D = dim(g)
    xyzt_pgops = sort(xyzt.(pointgroup(g)))

    @inbounds for iuclab in PGS_IUCs[D]
        P = pointgroup(iuclab, D)
        if sort(xyzt.(P))==(xyzt_pgops)
            return P
        end
    end
    return nothing

end

# --- POINT GROUP IRREPS ---
const DATA_PATH_PGIRREPS_3D = (@__DIR__)*"/../data/pgirreps/3d/irreps_data.jld2"

function _load_pgirreps_data(iuclab::String, jldfile::JLD2.JLDFile)
    jldgroup = jldfile[unmangle_pgiuclab(iuclab)] 
    matrices::Vector{Vector{Matrix{ComplexF64}}} = jldgroup["matrices"]
    types::Vector{Int64}                         = jldgroup["types"]
    cdmls::Vector{String}                        = jldgroup["cdmls"]

    return matrices, types, cdmls
end

# 3D
function get_pgirreps(iuclab::String, ::Val{3})
    pg = pointgroup(iuclab, Val(3)) # operations

    matrices, types, cdmls = JLD2.jldopen(DATA_PATH_PGIRREPS_3D, "r") do irs_jldfile
        _load_pgirreps_data(iuclab, irs_jldfile) # irrep matrices, types, & labels
    end
    
    return PGIrrep{3}.(cdmls, Ref(pg), matrices, types)
end
# 2D
function get_pgirreps(iuclab::String, ::Val{2})
    pg = pointgroup(iuclab, Val(2)) # operations

    # Because the operator sorting and setting is identical* between the shared point groups
    # of 2D and 3D, we can just do a whole-sale transfer of shared irreps from 3D to 2D.
    # (*) Actually, "2" and "m" have different settings in 2D and 3D; but they just have two
    #     operators and irreps each, so the setting difference doesn't matter.
    #     That the settings and sorting indeed agree between 2D and 3D is tested in 
    #     scripts/compare_pgops_3dvs2d.jl
    matrices, types, cdmls = JLD2.jldopen(DATA_PATH_PGIRREPS_3D, "r") do irs_jldfile
        _load_pgirreps_data(iuclab, irs_jldfile) # irrep matrices, types, & labels
    end
    
    return PGIrrep{2}.(cdmls, Ref(pg), matrices, types)
end
# 1D
function get_pgirreps(iuclab::String, ::Val{1})
    pg = pointgroup(iuclab, Val(1))
    # Situation in 1D is sufficiently simple that we don't need to bother with loading from 
    # a disk; just branch on one of the two possibilities
    if iuclab == "1"
        matrices = [[fill(one(ComplexF64), 1, 1)]]
        cdmls    = ["Γ₁"]
    elseif iuclab == "m"
        matrices = [[fill(one(ComplexF64), 1, 1),  fill(one(ComplexF64), 1, 1)], # even
                    [fill(one(ComplexF64), 1, 1), -fill(one(ComplexF64), 1, 1)]] # odd
        cdmls    = ["Γ₁", "Γ₂"]
    else
        throw(DomainError(iuclab, "invalid 1D point group IUC label"))
    end
    return PGIrrep{1}.(cdmls, Ref(pg), matrices, Ref(1))
end
get_pgirreps(iuclab::String, ::Val{D}) where D = _throw_invaliddim(D)
get_pgirreps(iuclab::String, D::Integer) = get_pgirreps(iuclab, Val(D))