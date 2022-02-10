using Crystalline

function make_lgirrep(cdml_suffix::String, lg::LittleGroup{D},
                      scalars_or_matrices::Vector, reality_type::Reality=REAL, 
                      translations=nothing) where D
    cdml = klabel(lg)*cdml_suffix
    if eltype(scalars_or_matrices) <: Number
        matrices = [fill(ComplexF64(v), 1,1) for v in scalars_or_matrices]
    else
        matrices = scalars_or_matrices
    end
    if isnothing(translations)
        translations = [zeros(D) for _ in matrices]
    end
    return LGIrrep{D}(cdml, lg, matrices, translations, reality_type, false)
end

# -----------------------------------------------------------------------------------------
# --- MANUAL TABULATION OF IRREPS OF NONSYMMORPHIC PLANE GROUPS (4, 7, 8, & 12) -----------
# In all cases, the irreps (and associated labels) were obtained simply by comparison with
# their parent space group. The settings are identical between these plane groups and their
# parent space group, except for plane group 4, which corresponds to an yz-plane instead of
# an xy-plane of its parent.

LGIRSD_2D = Dict{Int, Dict{String, Vector{LGIrrep{2}}}}()

# --- Plane group 4 -----------------------------------------------------------------------
# The yz-plane of space group 7
LGIRSD_2D[4] = Dict{String, Vector{LGIrrep{2}}}()
# Γ
kv = KVec(0,0)
lg = LittleGroup(4, kv, "Γ",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["Γ"] = [
    make_lgirrep("₁", lg, [1, +1], REAL)
    make_lgirrep("₂", lg, [1, -1], REAL)
    ]
# X (from 3D B-point)
kv = KVec(0,1/2)
lg = LittleGroup(4, kv, "X",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["X"] = [
    make_lgirrep("₁", lg, [1, 1im],  COMPLEX)
    make_lgirrep("₂", lg, [1, -1im], COMPLEX)
    ]
# Y (from 3D Z-point)
kv = KVec(1/2,0)
lg = LittleGroup(4, kv, "Y",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["Y"] = [
    make_lgirrep("₁", lg, [1, 1],  REAL)
    make_lgirrep("₂", lg, [1, -1], REAL)
    ]
# S (from 3D D-point)
kv = KVec(1/2,1/2)
lg = LittleGroup(4, kv, "S",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["S"] = [
    make_lgirrep("₁", lg, [1, 1im],  COMPLEX)
    make_lgirrep("₂", lg, [1, -1im], COMPLEX)
    ]
# Σ (from 3D F-point)
kv = KVec("0,u")
lg = LittleGroup(4, kv, "Σ",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["Σ"] = [
    make_lgirrep("₁", lg, [1, 1],  COMPLEX, [[0.0,0.0], [0.0,0.5]])
    make_lgirrep("₂", lg, [1, -1], COMPLEX, [[0.0,0.0], [0.0,0.5]])
    ]
# ΣA (from Σ)
kv = KVec("0,-u")
lg = LittleGroup(4, kv, "ΣA",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["ΣA"] = [
    make_lgirrep("₁", lg, [1, 1],  COMPLEX, [[0.0,0.0], [0.0,0.5]]) # same as Σ ...
    make_lgirrep("₂", lg, [1, -1], COMPLEX, [[0.0,0.0], [0.0,0.5]])
    ]
# C (from 3D G-point)
kv = KVec("1/2,u")
lg = LittleGroup(4, kv, "C",  [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["C"] = [
    make_lgirrep("₁", lg, [1, 1],  COMPLEX, [[0.0,0.0], [0.0,0.5]]) # same as Σ ...
    make_lgirrep("₂", lg, [1, -1], COMPLEX, [[0.0,0.0], [0.0,0.5]])
    ]
# CA (from C)
kv = KVec("1/2,-u")
lg = LittleGroup(4, kv, "CA", [S"x,y", S"-x,y+1/2"])
LGIRSD_2D[4]["CA"] = [
    make_lgirrep("₁", lg, [1, 1],  COMPLEX, [[0.0,0.0], [0.0,0.5]]) # same as Σ ...
    make_lgirrep("₂", lg, [1, -1], COMPLEX, [[0.0,0.0], [0.0,0.5]])
    ]


# --- Plane group 7 -----------------------------------------------------------------------
# The xz-plane of space group 28
LGIRSD_2D[7] = Dict{String, Vector{LGIrrep{2}}}()
# Γ
kv = KVec(0,0)
lg = LittleGroup(7, kv, "Γ", [S"x,y", S"-x,-y", S"-x+1/2,y", S"x+1/2,-y"])
LGIRSD_2D[7]["Γ"] = [
    make_lgirrep("₁", lg, [1, +1, +1, +1], REAL)
    make_lgirrep("₂", lg, [1, +1, -1, -1], REAL)
    make_lgirrep("₃", lg, [1, -1, +1, -1], REAL)
    make_lgirrep("₄", lg, [1, -1, -1, +1], REAL)
    ]
# Y
kv = KVec(0,1/2)
lg = LittleGroup(7, kv, "Y", [S"x,y", S"-x,-y", S"-x+1/2,y", S"x+1/2,-y"])
LGIRSD_2D[7]["X"] = [
    make_lgirrep("₁", lg, [1, +1, +1, +1], REAL)
    make_lgirrep("₂", lg, [1, +1, -1, -1], REAL)
    make_lgirrep("₃", lg, [1, -1, +1, -1], REAL)
    make_lgirrep("₄", lg, [1, -1, -1, +1], REAL)
    ]
# X
kv = KVec(1/2,0)
lg = LittleGroup(7, kv, "X", [S"x,y", S"-x,-y", S"-x+1/2,y", S"x+1/2,-y"])
LGIRSD_2D[7]["X"] = [
    make_lgirrep("₁", lg, [[1 0; 0 1], [0 -1; -1 0], [1 0; 0 -1], [0 -1; 1 0]], REAL)
    ]
# S
kv = KVec(1/2,1/2)
lg = LittleGroup(7, kv, "S", [S"x,y", S"-x,-y", S"-x+1/2,y", S"x+1/2,-y"])
LGIRSD_2D[7]["S"] = [
    make_lgirrep("₁", lg, [[1 0; 0 1], [0 -1; -1 0], [1 0; 0 -1], [0 -1; 1 0]], REAL)
    ]
# Σ
kv = KVec("u,0")
lg = LittleGroup(7, kv, "Σ", [S"x,y", S"x+1/2,-y"])
LGIRSD_2D[7]["Σ"] = [
    make_lgirrep("₁", lg, [1, +1], REAL, [[0.0,0.0], [0.5,0.0]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.0]])
    ]
# C
kv = KVec("u,1/2")
lg = LittleGroup(7, kv, "C", [S"x,y", S"x+1/2,-y"])
LGIRSD_2D[7]["C"] = [
    make_lgirrep("₁", lg, [1, +1], REAL, [[0.0,0.0], [0.5,0.0]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.0]])
    ]
# Δ
kv = KVec("0,u")
lg = LittleGroup(7, kv, "Δ", [S"x,y", S"-x+1/2,y"])
LGIRSD_2D[7]["Δ"] = [
    make_lgirrep("₁", lg, [1, -1], REAL)
    make_lgirrep("₂", lg, [1, +1], REAL)
    ]
# D
kv = KVec("1/2,u")
lg = LittleGroup(7, kv, "D", [S"x,y", S"-x+1/2,y"])
LGIRSD_2D[7]["D"] = [
    make_lgirrep("₁", lg, [1, +1], COMPLEX)
    make_lgirrep("₂", lg, [1, -1], COMPLEX)
    ]


# --- Plane group 8 -----------------------------------------------------------------------
# The xz-plane of space group 32
LGIRSD_2D[8] = Dict{String, Vector{LGIrrep{2}}}()
# Γ
kv = KVec(0,0)
lg = LittleGroup(8, kv, "Γ", [S"x,y", S"-x,-y", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["Γ"] = [
    make_lgirrep("₁", lg, [1, +1, +1, +1], REAL)
    make_lgirrep("₂", lg, [1, +1, -1, -1], REAL)
    make_lgirrep("₃", lg, [1, -1, +1, -1], REAL)
    make_lgirrep("₄", lg, [1, -1, -1, +1], REAL)
    ]
# Y
kv = KVec(0,1/2)
lg = LittleGroup(8, kv, "Y", [S"x,y", S"-x,-y", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["Y"] = [
    make_lgirrep("₁", lg, [[1 0; 0 1], [0 -1; -1 0], [0 -1; 1 0], [1 0; 0 -1]], REAL)
    ]
# X
kv = KVec(1/2,0)
lg = LittleGroup(8, kv, "X", [S"x,y", S"-x,-y", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["X"] = [
    make_lgirrep("₁", lg, [[1 0; 0 1], [0 -1; -1 0], [1 0; 0 -1], [0 -1; 1 0]], REAL)
    ]
# S
kv = KVec(1/2,1/2)
lg = LittleGroup(8, kv, "S", [S"x,y", S"-x,-y", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["S"] = [
    make_lgirrep("₁", lg, [1, +1, +1im, +1im], COMPLEX)
    make_lgirrep("₂", lg, [1, +1, -1im, -1im], COMPLEX)
    make_lgirrep("₃", lg, [1, -1, +1im, -1im], COMPLEX)
    make_lgirrep("₄", lg, [1, -1, -1im, +1im], COMPLEX)
    ]
# Σ
kv = KVec("u,0")
lg = LittleGroup(8, kv, "Σ", [S"x,y", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["Σ"] = [
    make_lgirrep("₁", lg, [1, +1], REAL, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.5]])
    ]
# C
kv = KVec("u,1/2")
lg = LittleGroup(8, kv, "C", [S"x,y", S"x+1/2,-y+1/2"])
LGIRSD_2D[8]["C"] = [
    make_lgirrep("₁", lg, [1, -1im], COMPLEX, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, +1im], COMPLEX, [[0.0,0.0], [0.5,0.5]])
    ]
# Δ
kv = KVec("0,u")
lg = LittleGroup(8, kv, "Δ", [S"x,y", S"-x+1/2,y+1/2"])
LGIRSD_2D[8]["Δ"] = [
    make_lgirrep("₁", lg, [1, +1], REAL, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.5]])
    ]
# D
kv = KVec("1/2,u")
lg = LittleGroup(8, kv, "D", [S"x,y", S"-x+1/2,y+1/2"])
LGIRSD_2D[8]["D"] = [
    make_lgirrep("₁", lg, [1, -1im], COMPLEX, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, +1im], COMPLEX, [[0.0,0.0], [0.5,0.5]])
    ]


# --- Plane group 12 ----------------------------------------------------------------------
# The xy-plane of space group 100
LGIRSD_2D[12] = Dict{String, Vector{LGIrrep{2}}}()
# Γ
kv = KVec(0,0)
lg = LittleGroup(12, kv, "Γ", [S"x,y", S"-x,-y", S"-y,x", S"y,-x", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2", S"-y+1/2,-x+1/2", S"y+1/2,x+1/2"])
LGIRSD_2D[12]["Γ"] = [
    make_lgirrep("₁", lg, [1, +1, +1, +1, +1, +1, +1, +1], REAL)
    make_lgirrep("₂", lg, [1, +1, -1, -1, +1, +1, -1, -1], REAL)
    make_lgirrep("₃", lg, [1, +1, -1, -1, -1, -1, +1, +1], REAL)
    make_lgirrep("₄", lg, [1, +1, +1, +1, -1, -1, -1, -1], REAL)
    make_lgirrep("₅", lg, [[1 0; 0 1], [-1 0; 0 -1], [0 -1; 1 0], [0 1; -1 0], [-1 0; 0 1], [1 0; 0 -1], [0 -1; -1 0], [0 1; 1 0]], REAL)
    ]
# M
kv = KVec(1/2,1/2)
lg = LittleGroup(12, kv, "M", [S"x,y", S"-x,-y", S"-y,x", S"y,-x", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2", S"-y+1/2,-x+1/2", S"y+1/2,x+1/2"])
LGIRSD_2D[12]["M"] = [
    make_lgirrep("₁", lg, [1, -1, 1im, -1im, 1im, -1im, 1, -1], COMPLEX)
    make_lgirrep("₂", lg, [1, -1, -1im, 1im, 1im, -1im, -1, 1], COMPLEX)
    make_lgirrep("₃", lg, [1, -1, -1im, 1im, -1im, 1im, 1, -1], COMPLEX)
    make_lgirrep("₄", lg, [1, -1, 1im, -1im, -1im, 1im, -1, 1], COMPLEX)
    make_lgirrep("₅", lg, [[1 0; 0 1], [1 0; 0 1], [-1 0; 0 1], [-1 0; 0 1], [0 -1; 1 0], [0 -1; 1 0], [0 -1; -1 0], [0 -1; -1 0]])
    ]
# X
kv = KVec(0,1/2)
lg = LittleGroup(12, kv, "X", [S"x,y", S"-x,-y", S"-x+1/2,y+1/2", S"x+1/2,-y+1/2"])
LGIRSD_2D[12]["X"] = [
    make_lgirrep("₁", lg, [[1 0; 0 1], [1 0; 0 -1], [0 -1; 1 0], [0 1; 1 0]])
    ]
# Δ
kv = KVec("0,v")
lg = LittleGroup(12, kv, "Δ", [S"x,y", S"-x+1/2,y+1/2"])
LGIRSD_2D[12]["Δ"] = [
    make_lgirrep("₁", lg, [1, 1],  REAL, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.5]])
    ]
# Y
kv = KVec("u,1/2")
lg = LittleGroup(12, kv, "Y", [S"x,y", S"x+1/2,-y+1/2"])
LGIRSD_2D[12]["Y"] = [
    make_lgirrep("₁", lg, [1, -1im], COMPLEX, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, 1im],  COMPLEX, [[0.0,0.0], [0.5,0.5]])
    ]
# Σ 
kv = KVec("u,u")
lg = LittleGroup(12, kv, "Σ", [S"x,y", S"y+1/2,x+1/2"])
LGIRSD_2D[12]["Σ"] = [
    make_lgirrep("₁", lg, [1, 1],  REAL, [[0.0,0.0], [0.5,0.5]])
    make_lgirrep("₂", lg, [1, -1], REAL, [[0.0,0.0], [0.5,0.5]])
    ]


# -----------------------------------------------------------------------------------------
# ADD TRIVIAL Ω IRREP TO EACH LISTING
for (sgnum, lgirsd) in LGIRSD_2D
    local lg = LittleGroup(sgnum, KVec("u,v"), "Ω", [S"x,y"])
    lgirsd["Ω"] = [make_lgirrep("₁", lg, [1], sgnum == 4 ? COMPLEX : REAL)]
    # (reality is COMPLEX at Ω for plane group 4, and REAL for plane groups 7, 8, and 12)
end

# -----------------------------------------------------------------------------------------
# TEST FOR SUBSET OF POSSIBLE TYPOS
for (sgnum, lgirsd) in LGIRSD_2D
    sg = reduce_ops(spacegroup(sgnum, Val(2)), true)
    for (klab,lgirs) in lgirsd
        g = group(first(lgirs))
        if klabel(g) != klab
            error("Tabulation error: mismatched k-label $(klabel(g)) and dict-index $(klab)")
        end
        if reality.(lgirs) != calc_reality.(lgirs, Ref(sg), Ref(Crystalline.TEST_αβγs[2]))
            error("Tabulation error: reality type")
        end
    end
end