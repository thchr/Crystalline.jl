using Pkg
filedir = joinpath(@__DIR__)
(dirname(Pkg.project().path) == filedir) || Pkg.activate(filedir)

using Crystalline
using BandGraphs
using BandGraphs: Connection, LabeledKVec, SubductionTable
using JLD2

# ---------------------------------------------------------------------------------------- #

const C2 = Connection{2}
const LK2 = LabeledKVec{2}
LK2(label::Symbol, kstr::String) = LK2(label, KVec{2}(kstr))
const DIR = joinpath((@__DIR__), "..", "data", "connections", "2d")

# ---------------------------------------------------------------------------------------- #

# Below, "manual addition" = "manual addition relative to parent space group" indicates
# connections that have been added for the plane group, but which are not included in the 
# parent space group; these deviations are made in the hope that the human analysis is eased

# No distinctions are presently made between TR-invariant and TR-broken connections in 2D

# ---------------------------------------------------------------------------------------- #

CONNECTIONSD_2D = Dict{Int, Vector{Connection{2}}}()

# Plane group 1 (p1)
CONNECTIONSD_2D[1] = C2[]

# Plane group 2 (p2)
CONNECTIONSD_2D[2] = C2[
    C2(LK2(:Γ, "0, 0"),      LK2(:Ω, "α, β")), # Γ = [0, 0]      ↓ Ω = [α, β] (manual addition)
    C2(LK2(:B, "1/2, 0"),    LK2(:Ω, "α, β")), # B = [1/2, 0]    ↓ Ω = [α, β] (manual addition)
    C2(LK2(:Y, "0, 1/2"),    LK2(:Ω, "α, β")), # Y = [0, 1/2]    ↓ Ω = [α, β] (manual addition)
    C2(LK2(:A, "1/2, -1/2"), LK2(:Ω, "α, β")), # A = [1/2, -1/2] ↓ Ω = [α, β] (manual addition)
]

# Plane group 3 (p1m1)
CONNECTIONSD_2D[3] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "0, α")),   # Γ = [0, 0]     ↓ Σ = [0, α]
    C2(LK2(:X, "0, 1/2"),   LK2(:Σ, "0, α")),   # X = [0, 1/2]   ↓ Σ = [0, α]
    C2(LK2(:S, "1/2, 1/2"), LK2(:C, "1/2, α")), # S = [1/2, 1/2] ↓ C = [1/2, α]
    C2(LK2(:Y, "1/2, 0"),   LK2(:C, "1/2, α")), # Y = [1/2, 0]   ↓ C = [1/2, α]
]

# Plane group 4 (p1g1)
CONNECTIONSD_2D[4] = C2[
    C2(LK2(:Γ, "0, 0"),      LK2(:Σ, "0, α")),   # Γ = [0, 0]      ↓ Σ = [0, α]
    C2(LK2(:X, "0, 1/2"),    LK2(:Σ, "0, α")),   # X = [0, 1/2]    ↓ Σ = [0, α]
    C2(LK2(:S, "1/2, 1/2"),  LK2(:C, "1/2, α")), # S = [1/2, 1/2]  ↓ C = [1/2, α]
    C2(LK2(:Y, "1/2, 0"),    LK2(:C, "1/2, α")), # Y = [1/2, 0]    ↓ C = [1/2, α]
    C2(LK2(:X′, "0, 3/2"),   LK2(:Σ, "0, α")),   # X′ = [0, 3/2]   ↓ Σ = [0, α]
    C2(LK2(:Γ′, "0, 1"),     LK2(:Σ, "0, α")),   # Γ′ = [0, 1]     ↓ Σ = [0, α]
    C2(LK2(:S′, "1/2, 3/2"), LK2(:C, "1/2, α")), # S′ = [1/2, 3/2] ↓ C = [1/2, α]
    C2(LK2(:Y′, "1/2, 1"),   LK2(:C, "1/2, α")), # Y′ = [1/2, 1]   ↓ C = [1/2, α]
]

# Plane group 5 (c1m1)
CONNECTIONSD_2D[5] = C2[
    C2(LK2(:Γ, "0, 0"), LK2(:Σ, "0, 2α")), # Γ = [0, 0] ↓ Σ = [0, 2α]
    C2(LK2(:Y, "0, 1"), LK2(:Σ, "0, 2α")), # Y = [0, 1] ↓ Σ = [0, 2α]
]

# Plane group 6 (p2mm)
CONNECTIONSD_2D[6] = C2[
    C2(LK2(:S, "1/2, 1/2"), LK2(:C, "α, 1/2")), # S = [1/2, 1/2] ↓ C = [α, 1/2]
    C2(LK2(:Y, "0, 1/2"),   LK2(:C, "α, 1/2")), # Y = [0, 1/2]   ↓ C = [α, 1/2]
    C2(LK2(:S, "1/2, 1/2"), LK2(:D, "1/2, α")), # S = [1/2, 1/2] ↓ D = [1/2, α]
    C2(LK2(:X, "1/2, 0"),   LK2(:D, "1/2, α")), # X = [1/2, 0]   ↓ D = [1/2, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Δ, "0, α")),   # Γ = [0, 0]     ↓ Δ = [0, α]
    C2(LK2(:Y, "0, 1/2"),   LK2(:Δ, "0, α")),   # Y = [0, 1/2]   ↓ Δ = [0, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, 0")),   # Γ = [0, 0]     ↓ Σ = [α, 0]
    C2(LK2(:X, "1/2, 0"),   LK2(:Σ, "α, 0")),   # X = [1/2, 0]   ↓ Σ = [α, 0]
]

# Plane group 7 (p2mg)
CONNECTIONSD_2D[7] = C2[
    C2(LK2(:S, "1/2, 1/2"),  LK2(:C, "α, 1/2")), # S = [1/2, 1/2]  ↓ C = [α, 1/2]
    C2(LK2(:Y, "0, 1/2"),    LK2(:C, "α, 1/2")), # Y = [0, 1/2]    ↓ C = [α, 1/2]
    C2(LK2(:S, "1/2, 1/2"),  LK2(:D, "1/2, α")), # S = [1/2, 1/2]  ↓ D = [1/2, α]
    C2(LK2(:X, "1/2, 0"),    LK2(:D, "1/2, α")), # X = [1/2, 0]    ↓ D = [1/2, α]
    C2(LK2(:Γ, "0, 0"),      LK2(:Δ, "0, α")),   # Γ = [0, 0]      ↓ Δ = [0, α]
    C2(LK2(:Y, "0, 1/2"),    LK2(:Δ, "0, α")),   # Y = [0, 1/2]    ↓ Δ = [0, α]
    C2(LK2(:Γ, "0, 0"),      LK2(:Σ, "α, 0")),   # Γ = [0, 0]      ↓ Σ = [α, 0]
    C2(LK2(:X, "1/2, 0"),    LK2(:Σ, "α, 0")),   # X = [1/2, 0]    ↓ Σ = [α, 0]
    C2(LK2(:S′, "3/2, 1/2"), LK2(:C, "α, 1/2")), # S′ = [3/2, 1/2] ↓ C = [α, 1/2]
    C2(LK2(:Y′, "1, 1/2"),   LK2(:C, "α, 1/2")), # Y′ = [1, 1/2]   ↓ C = [α, 1/2]
    C2(LK2(:Γ′, "1, 0"),     LK2(:Σ, "α, 0")),   # Γ′ = [1, 0]     ↓ Σ = [α, 0]
    C2(LK2(:X′, "3/2, 0"),   LK2(:Σ, "α, 0")),   # X′ = [3/2, 0]   ↓ Σ = [α, 0]
]

# Plane group 8 (p2gg)
CONNECTIONSD_2D[8] = C2[
    C2(LK2(:S, "1/2, 1/2"), LK2(:C, "α, 1/2")), # S = [1/2, 1/2] ↓ C = [α, 1/2]
    C2(LK2(:Y, "0, 1/2"),   LK2(:C, "α, 1/2")), # Y = [0, 1/2]   ↓ C = [α, 1/2]
    C2(LK2(:S, "1/2, 1/2"), LK2(:D, "1/2, α")), # S = [1/2, 1/2] ↓ D = [1/2, α]
    C2(LK2(:X, "1/2, 0"),   LK2(:D, "1/2, α")), # X = [1/2, 0]   ↓ D = [1/2, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Δ, "0, α")),   # Γ = [0, 0]     ↓ Δ = [0, α]
    C2(LK2(:Y, "0, 1/2"),   LK2(:Δ, "0, α")),   # Y = [0, 1/2]   ↓ Δ = [0, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, 0")),   # Γ = [0, 0]     ↓ Σ = [α, 0]
    C2(LK2(:X, "1/2, 0"),   LK2(:Σ, "α, 0")),   # X = [1/2, 0]   ↓ Σ = [α, 0]
    C2(LK2(:Γ′, "[0, 1"),   LK2(:Δ, "0, α")),   # Γ′ = [0, 1]    ↓ Δ = [0, α]
    C2(LK2(:Y′, "[0, 3/2"), LK2(:Δ, "0, α")),   # Y′ = [0, 3/2]  ↓ Δ = [0, α]
    C2(LK2(:Γ′, "[1, 0"),   LK2(:Σ, "α, 0")),   # Γ′ = [1, 0]    ↓ Σ = [α, 0]
    C2(LK2(:X′, "[3/2, 0"), LK2(:Σ, "α, 0")),   # X′ = [3/2, 0]  ↓ Σ = [α, 0]
]

# Plane group 9 (c2mm)
CONNECTIONSD_2D[9] = C2[
    C2(LK2(:Γ, "0, 0"), LK2(:Δ, "0, 2α")), # Γ = [0, 0] ↓ Δ = [0, 2α]
    C2(LK2(:Y, "1, 0"), LK2(:Δ, "0, 2α")), # Y = [1, 0] ↓ Δ = [0, 2α]
    C2(LK2(:Γ, "0, 0"), LK2(:Σ, "2α, 0")), # Γ = [0, 0] ↓ Σ = [2α, 0]
    C2(LK2(:Y, "1, 0"), LK2(:Σ, "2α, 0")), # Y = [1, 0] ↓ Σ = [2α, 0]
    C2(LK2(:S, "1, 2"), LK2(:Ω, "α, β")),  # S = [1, 2] ↓ Ω = [α, β] (manual addition)
    C2(LK2(:Γ, "0, 0"), LK2(:Ω, "α, β")),  # Γ = [0, 0] ↓ Ω = [α, β] (manual addition)
]

# Plane group 10 (p4)
CONNECTIONSD_2D[10] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Ω, "α, β")), # Γ = [0, 0]     ↓ Ω = [α, β] (manual addition)
    C2(LK2(:X, "0, 1/2"),   LK2(:Ω, "α, β")), # X = [0, 1/2]   ↓ Ω = [α, β] (manual addition)
    C2(LK2(:M, "1/2, 1/2"), LK2(:Ω, "α, β")), # M = [1/2, 1/2] ↓ Ω = [α, β] (manual addition)
]

# Plane group 11 (p4mm)
CONNECTIONSD_2D[11] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Δ, "0, α")),   # Γ = [0, 0]     ↓ Δ = [0, α]
    C2(LK2(:X, "0, 1/2"),   LK2(:Δ, "0, α")),   # X = [0, 1/2]   ↓ Δ = [0, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, α")),   # Γ = [0, 0]     ↓ Σ = [α, α]
    C2(LK2(:M, "1/2, 1/2"), LK2(:Σ, "α, α")),   # M = [1/2, 1/2] ↓ Σ = [α, α]
    C2(LK2(:M, "1/2, 1/2"), LK2(:Y, "α, 1/2")), # M = [1/2, 1/2] ↓ Y = [α, 1/2]
    C2(LK2(:X, "0, 1/2"),   LK2(:Y, "α, 1/2")), # X = [0, 1/2]   ↓ Y = [α, 1/2]
]

# Plane group 12 (p4gm)
CONNECTIONSD_2D[12] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Δ, "0, α")),   # Γ = [0, 0]     ↓ Δ = ↓ [0, α]
    C2(LK2(:X, "0, 1/2"),   LK2(:Δ, "0, α")),   # X = [0, 1/2]   ↓ Δ = ↓ [0, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, α")),   # Γ = [0, 0]     ↓ Σ = ↓ [α, α]
    C2(LK2(:M, "1/2, 1/2"), LK2(:Σ, "α, α")),   # M = [1/2, 1/2] ↓ Σ = ↓ [α, α]
    C2(LK2(:M, "1/2, 1/2"), LK2(:Y, "α, 1/2")), # M = [1/2, 1/2] ↓ Y = ↓ [α, 1/2]
    C2(LK2(:X, "0, 1/2"),   LK2(:Y, "α, 1/2")), # X = [0, 1/2]   ↓ Y = ↓ [α, 1/2]
    C2(LK2(:Γ′, "0, 1"),    LK2(:Δ, "0, α")),   # Γ′ = [0, 1]    ↓ Δ = ↓ [0, α]
    C2(LK2(:X′, "0, 3/2"),  LK2(:Δ, "0, α")),   # X′ = [0, 3/2]  ↓ Δ = ↓ [0, α]
]

# Plane group 13 (p3)
CONNECTIONSD_2D[13] = C2[
    C2(LK2(:Γ,  "0, 0"),      LK2(:Ω, "α, β")), # Γ  = [0, 0]      ↓ Ω = [α, β] (manual addition)
    C2(LK2(:K,  "1/3, 1/3"),  LK2(:Ω, "α, β")), # K  = [1/3, 1/3]  ↓ Ω = [α, β] (manual addition)
    C2(LK2(:KA, "2/3, -1/3"), LK2(:Ω, "α, β")), # KA = [2/3, -1/3] ↓ Ω = [α, β] (manual addition)
]

# Plane group 14 (p3m1)
CONNECTIONSD_2D[14] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, 0")), # Γ = [0, 0]     ↓ Σ = [α, 0]
    C2(LK2(:M, "1/2, 0"),   LK2(:Σ, "α, 0")), # M = [1/2, 0]   ↓ Σ = [α, 0]
    C2(LK2(:K, "1/3, 1/3"), LK2(:Ω, "α, β")), # K = [1/3, 1/3] ↓ Ω = [α, β] (manual addition)
    C2(LK2(:Γ, "0, 0"),     LK2(:Ω, "α, β")), # Γ = [0, 0]     ↓ Ω = [α, β] (manual addition)
]

# Plane group 15 (p31m)
CONNECTIONSD_2D[15] = C2[
    C2(LK2(:Γ, "0, 0"),        LK2(:Λ, "-2α, α")), # Γ = [0, 0]        ↓ Λ = [-2α, α]
    C2(LK2(:K, "1/3, 1/3"),    LK2(:Λ, "-2α, α")), # K = [1/3, 1/3]    ↓ Λ = [-2α, α]
    C2(LK2(:Γ, "0, 0"),        LK2(:Λ, "-2α, α")), # Γ = [0, 0]        ↓ Λ = [-2α, α]
    C2(LK2(:KA, "-1/3, -1/3"), LK2(:Λ, "-2α, α")), # KA = [-1/3, -1/3] ↓ Λ = [-2α, α]
    C2(LK2(:Γ, "0, 0"),        LK2(:Λ, "α, -2α")), # Γ = [0, 0]        ↓ Λ = [u, -2α]
    C2(LK2(:M, "1/2, 0"),      LK2(:Λ, "α, -2α")), # M = [1/2, 0]      ↓ Λ = [u, -2α]
]

# Plane group 16 (p6)
CONNECTIONSD_2D[16] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Ω, "α, β")), # Γ = [0, 0]     ↓ Ω = [α, β] (manual addition)
    C2(LK2(:K, "1/3, 1/3"), LK2(:Ω, "α, β")), # K = [1/3, 1/3] ↓ Ω = [α, β] (manual addition)
    C2(LK2(:M, "1/2, 0"),   LK2(:Ω, "α, β")), # M = [1/2, 0]   ↓ Ω = [α, β] (manual addition)
]

# Plane group 17 (p6mm)
CONNECTIONSD_2D[17] = C2[
    C2(LK2(:Γ, "0, 0"),     LK2(:Λ, "α, α")),      # Γ = [0, 0]     ↓ Λ = [α, α]
    C2(LK2(:K, "1/3, 1/3"), LK2(:Λ, "α, α")),      # K = [1/3, 1/3] ↓ Λ = [α, α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Λ, "α, -2α")),    # Γ = [0, 0]     ↓ Λ = [α, -2α]
    C2(LK2(:M, "1/2, 0"),   LK2(:Λ, "α, -2α")),    # M = [1/2, 0]   ↓ Λ = [α, -2α]
    C2(LK2(:Γ, "0, 0"),     LK2(:Σ, "α, 0")),      # Γ = [0, 0]     ↓ Σ = [α, 0]
    C2(LK2(:M, "1/2, 0"),   LK2(:Σ, "α, 0")),      # M = [1/2, 0]   ↓ Σ = [α, 0]
    C2(LK2(:Γ, "1/2, 0"),   LK2(:T, "1/2-α, 2α")), # Γ = [1/2, 0]   ↓ T = [1/2-α, 2α] (manual addition)
    C2(LK2(:M, "1/2, 0"),   LK2(:T, "1/2-α, 2α")), # M = [1/2, 0]   ↓ T = [1/2-α, 2α] (manual addition)
    C2(LK2(:K, "1/3, 1/3"), LK2(:T, "1/2-α, 2α")), # K = [1/3, 1/3] ↓ T = [1/2-α, 2α] (manual addition)
]

# save connections data
for timereversal in (false, true)
    # we do not currently make any distinction between TR-invariant and -broken connections
    jldsave(joinpath(DIR, "connections$(timereversal ? "-tr" : "").jld2"); 
            connectionsd=CONNECTIONSD_2D)
end

## --------------------------------------------------------------------------------------- #
# Compute subduction tables for each connection above

SUBDUCTIONSD_2D_TR = Dict{Int, Vector{SubductionTable{2}}}()
SUBDUCTIONSD_2D    = Dict{Int, Vector{SubductionTable{2}}}()
for sgnum in 1:MAX_SGNUM[2]
    cs = CONNECTIONSD_2D[sgnum]
    lgirsd = lgirreps(sgnum, Val(2))
    sg = spacegroup(sgnum, Val(2))
    for timereversal in (false, true)
        timereversal && (lgirsd = realify(lgirsd))
        subts = SubductionTable{2}[]
        for c in cs
            push!(subts, SubductionTable(c, sg, lgirsd))
        end

        if timereversal
            SUBDUCTIONSD_2D_TR[sgnum] = subts
        else
            SUBDUCTIONSD_2D[sgnum] = subts
        end
    end
end

# save subduction-table data
jldsave(joinpath(DIR, "subductions.jld2");    subductionsd=SUBDUCTIONSD_2D)
jldsave(joinpath(DIR, "subductions-tr.jld2"); subductionsd=SUBDUCTIONSD_2D_TR)