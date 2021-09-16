# The goal of this script is to tabulate and construct the irreps in Φ-Ω (basic domain, Ω; 
# representation domain, Φ) that are missing from ISOTROPY but still feature in the bandreps
# from Bilbao. By manual comparison, we found that there are 145 such irreps, across 20
# space groups, namely:
# ┌───────┬─────────────────────────────────────────────────────────────────────────────┬───────────────┬────────────────────────────────┬──────────────┐
# │       │  Missing                                                                    │  Missing      │  [NS≡nonsymmorph; S≡symmorph]  │  Has fragile │
# │  SGs  │  LGIrrep labels                                                             │  KVec labels  │  Match method                  │  phases?     │
# │───────┼─────────────────────────────────────────────────────────────────────────────┼───────────────┼────────────────────────────────┼──────────────│
# │  23   │  WA₁, WA₂, WA₃, WA₄                                                         │  WA           │  S:  MonoOrthTetraCubic ┐      │  ÷           │
# │  24   │  WA₁                                                                        │  WA           │  NS: Inherits from      └ 23   │  ÷           │
# │  82   │  PA₁, PA₂, PA₃, PA₄                                                         │  PA           │  S:  MonoOrthTetraCubic        │              │
# │  121  │  PA₁, PA₂, PA₃, PA₄, PA₅                                                    │  PA           │  S:  MonoOrthTetraCubic ┐      │              │
# │  122  │  PA₁, PA₂                                                                   │  PA           │  NS: Inherits from      └ 121  │              │
# │  143  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex                    │              │
# │  144  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┐         │  ÷           │
# │  145  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┘         │  ÷           │
# │  150  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex                    │              │
# │  152  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┐         │  ÷           │
# │  154  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┘         │              │
# │  157  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex ───────┐           │              │
# │  159  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Inherits from └ 157       │              │
# │  174  │  HA₁, HA₂, HA₃, HA₄, HA₅, HA₆, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃  │  HA, KA, PA   │  S:  TriHex                    │              │
# │  189  │  HA₁, HA₂, HA₃, HA₄, HA₅, HA₆, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃  │  HA, KA, PA   │  S:  TriHex ───────┐           │              │
# │  190  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃                 │  HA, KA, PA   │  NS: Inherits from └ 189       │              │
# │  197  │  PA₁, PA₂, PA₃, PA₄                                                         │  PA           │  S:  MonoOrthTetraCubic ┐      │  ÷           │
# │  199  │  PA₁, PA₂, PA₃                                                              │  PA           │  NS: Inherits from      └ 197  │  ÷           │
# │  217  │  PA₁, PA₂, PA₃, PA₄, PA₅                                                    │  PA           │  S:  MonoOrthTetraCubic ┐      │              │
# │  220  │  PA₁, PA₂, PA₃                                                              │  PA           │  NS: Inherits from      └ 217  │              │
# └───────┴─────────────────────────────────────────────────────────────────────────────┴───────────────┴────────────────────────────────┴──────────────┘
# In principle, these irreps _could_ be constructed from the existing irreps in ISOTROPY by
# suitable transformations, as described by B&C (e.g. near p. 414). Unfortunately, we have
# thus far not managed to implement that scheme correctly. Instead, we here just go the
# brute-force path and tabulate the things we need manually.
# NB: Of the above space groups, only SG 82 has nontrivial symmetry indicator for bosons
# with TR. Most of the space groups, however, can support symmetry-indicated fragile
# topology (÷ means 'cannot').

# ======================================================================================== #
# ================================== UTILITY FUNCTIONS =================================== #
# ======================================================================================== #

using Crystalline

# manually converting
#    https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_out.pl
# to our own data format for the missing k-points listed in 
#   src/special_representation_domain_kpoints.jl
# TODO: This needs to be merged into the special-points branch before master

function build_lgirrep_with_type(cdml, lg, Psτs, sgops)
    # input handling
    if Psτs isa Vector                    # assume zero-τ factors
        Ps = complex.(float.(Psτs))
        τs = nothing
    elseif Psτs isa Tuple{<:Any, <:Any}   # nonzero τ-factors; input as 2-tuple
        Ps = complex.(float.(Psτs[1]))
        τs = Psτs[2]
    else
        throw("Unexpected input format of Psτs")
    end
    # convert scalar irreps (numbers) to 1×1 matrices
    if eltype(Ps) <: Number   
        Ps = fill.(Ps, 1, 1)
    end

    lgir  = LGIrrep{3}(cdml, lg, Ps, τs, REAL)    # place-holder `REAL` reality type
    reality = calc_reality(lgir, sgops)
    lgir  = LGIrrep{3}(cdml, lg, Ps, τs, reality) # update reality type
end

function prepare_lg_and_sgops(sgnum, kv, klab, ops)
    lg    = LittleGroup{3}(sgnum, kv, klab, ops)
    sgops = reduce_ops(spacegroup(sgnum, Val(3)), centering(sgnum, 3), true) # for herrring
    return lg, sgops
end

function assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
    lg, sgops = prepare_lg_and_sgops(sgnum, kv, klab, lgops)
    cdmls     = Ref(klab) .* string.(1:length(Psτs))
    lgirs     = build_lgirrep_with_type.(cdmls, Ref(lg), Psτs, Ref(sgops))
end

if !isdefined(Main, :cispi)
    cispi(x) = cis(π*x)
end


# ======================================================================================== #
# =========================== MANUALLY COPIED IRREP DATA BELOW =========================== #
# ======================================================================================== #

# "preallocate" a storage dict
sgnums = [23, 24, 82, 121, 122, 143, 144, 145, 150, 152, 154, 157, 159, 174, 189, 190, 197, 
          199, 217, 220]
LGIRS_add = Dict(sgnum=>Dict{String, Vector{LGIrrep{3}}}() for sgnum in sgnums)

# ========= 23 =========
sgnum = 23
# WA₁, WA₂, WA₃, WA₄
klab  = "WA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "x,-y,-z", "-x,y,-z", "-x,-y,z"]) # 1, 2₁₀₀, 2₀₁₀, 2₀₀₁

# listed first across irrep (ascending, e.g. WA1, WA2, WA3, WA4 here) then across lgops
Psτs = [[1, 1, 1, 1],    # Ps = matrices
                         # nonzero τs (translations) can be specified by entering a 2-tuple of Ps and τs instead of a vector of Ps
        [1, -1, -1, 1],
        [1, 1, -1, -1],
        [1, -1, 1, -1],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 24 =========
sgnum = 24
# WA₁
klab  = "WA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "x,-y,-z+1/2", "-x+1/2,y,-z", "-x,-y+1/2,z"]) # 1, {2₁₀₀|00½}, {2₀₁₀|½00}, {2₀₀₁|0½0}
Psτs = [[[1 0; 0 1], [1 0; 0 -1], [0 -im; im 0], [0 1; 1 0]],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 82  =========
sgnum = 82
# PA₁, PA₂, PA₃, PA₄
klab  = "PA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z"]) # 1, 2₀₀₁, -4⁺₀₀₁, -4⁻₀₀₁
Psτs = [[1, 1, 1, 1],
        [1, 1, -1, -1],
        [1, -1, -im, im],
        [1, -1, im, -im],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 121 =========
sgnum = 121
# PA₁, PA₂, PA₃, PA₄, PA₅
klab  = "PA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z", "-x,y,-z", "x,-y,-z", "-y,-x,z", "y,x,z"]) 
                          # 1, 2₀₀₁, -4⁺₀₀₁, -4⁻₀₀₁, 2₀₁₀, 2₁₀₀, m₁₁₀, m₁₋₁₀
Psτs = [[1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, -1, -1, 1, 1, -1, -1],
        [1, 1, -1, -1, -1, -1, 1, 1],
        [1, 1, 1, 1, -1, -1, -1, -1],
        [[1 0; 0 1], [-1 0; 0 -1], [0 -1; 1 0], [0 1; -1 0], [0 1; 1 0], [0 -1; -1 0], [1 0; 0 -1], [-1 0; 0 1]],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 122 =========
sgnum = 122
# PA₁, PA₂
klab  = "PA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z", "-x,y+1/2,-z+1/4", "x,-y+1/2,-z+1/4", "-y,-x+1/2,z+1/4", "y,x+1/2,z+1/4"]) 
                          # 1, 2₀₀₁, -4⁺₀₀₁, -4⁻₀₀₁, {2₀₁₀|0,1/2,1/4}, {2₁₀₀|0,1/2,1/4}, {m₁₁₀|0,1/2,1/4}, {m₁₋₁₀|0,1/2,1/4}

Psτs = [[[1 0; 0 1], [1 0; 0 -1], [1 0; 0 im], [1 0; 0 -im], [0 -1; 1 0], [0 1; 1 0], [0 -im; 1 0], [0 im; 1 0]],
        [[1 0; 0 1], [1 0; 0 -1], [-1 0; 0 -im], [-1 0; 0 im], [0 1; -1 0], [0 -1; -1 0], [0 -im; 1 0], [0 im; 1 0]],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 143 =========
sgnum = 143
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z"]) # 1, 3⁺₀₀₁, 3⁻₀₀₁
# HA₁, HA₂, HA₃
klab  = "HA"
kv    = KVec(-1/3,-1/3,-1/2)
Psτs = [[1, 1, 1],
        [1, cispi(-2/3), cispi(2/3)],
        [1, cispi(2/3), cispi(-2/3)],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops & irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
# ... same lgops & irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 144 =========
sgnum = 144
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]) # 1, {3⁺₀₀₁|0,0,1/3}, {3⁻₀₀₁|0,0,2/3}
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
Psτs = [[1, cispi(-1/3), cispi(-2/3)],
        [1, -1, 1],
        [1, cispi(1/3), cispi(2/3)],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
Psτs = [[1, 1, 1],
        [1, cispi(-2/3), cispi(2/3)],
        [1, cispi(2/3), cispi(-2/3)],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
Psτs = [([1, 1, 1],                    [[0.0,0,0], [0,0,1/3], [0,0,2/3]]), # nonzero τs
        ([1, cispi(-2/3), cispi(2/3)], [[0.0,0,0], [0,0,1/3], [0,0,2/3]]),
        ([1, cispi(2/3), cispi(-2/3)], [[0.0,0,0], [0,0,1/3], [0,0,2/3]]),]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 145 =========
sgnum = 145
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3"]) # 1, {3⁺₀₀₁|0,0,2/3}, {3⁻₀₀₁|0,0,1/3}
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
Psτs = [[1, 1, -1],
        [1, cispi(-2/3), cispi(-1/3)],
        [1, cispi(2/3), cispi(1/3)],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [[1, 1, 1],
        [1, cispi(-2/3), cispi(2/3)],
        [1, cispi(2/3), cispi(-2/3)],]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [([1, 1, 1],                    [[0.0,0,0], [0,0,2/3], [0,0,1/3]]), # nonzero τs
        ([1, cispi(-2/3), cispi(2/3)], [[0.0,0,0], [0,0,2/3], [0,0,1/3]]),
        ([1, cispi(2/3), cispi(-2/3)], [[0.0,0,0], [0,0,2/3], [0,0,1/3]]),]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 150 =========
sgnum = 150
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, 2₁₁₀, 2₁₀₀, 2₀₁₀
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops & irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻
Psτs = [[1, 1, 1],
        [1, cispi(-2/3), cispi(2/3)],
        [1, cispi(2/3), cispi(-2/3)], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 152 =========
sgnum = 152
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3", "y,x,-z", "x-y,-y,-z+2/3", "-x,-x+y,-z+1/3"]) # 1, {3₀₀₁⁺|0,0,⅓}, {3₀₀₁⁻|0,0,⅔}, 2₁₁₀, {2₁₀₀|0,0,⅔}, {2₀₁₀|0,0,⅓}
Psτs = [[1, -1, 1, 1, 1, -1],
        [1, -1, 1, -1, -1, 1],
        [[1 0; 0 1], [cispi(-1/3) 0; 0 cispi(1/3)], [cispi(-2/3) 0; 0 cispi(2/3)], [0 1; 1 0], [0 cispi(-2/3); cispi(2/3) 0], [0 cispi(-1/3); cispi(1/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]) # 1, {3₀₀₁⁺|0,0,⅓}, {3₀₀₁⁻|0,0,⅔}
Psτs = [([1, 1, 1],                    [[0.0,0,0], [0,0,1/3], [0,0,2/3]]),
        ([1, cispi(-2/3), cispi(2/3)], [[0.0,0,0], [0,0,1/3], [0,0,2/3]]), 
        ([1, cispi(2/3), cispi(-2/3)], [[0.0,0,0], [0,0,1/3], [0,0,2/3]]), ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 154 =========
sgnum = 154
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3", "y,x,-z", "x-y,-y,-z+1/3", "-x,-x+y,-z+2/3"]) # 1, {3₀₀₁⁺|0,0,⅔}, {3₀₀₁⁻|0,0,⅓}, 2₁₁₀, {2₁₀₀|0,0,⅓}, {2₀₁₀|0,0,⅔}
Psτs = [[1, 1, -1, -1, 1, -1],
        [1, 1, -1, 1, -1, 1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(-1/3) 0; 0 cispi(1/3)], [0 1; 1 0], [0 cispi(-1/3); cispi(1/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3"]) # 1, {3₀₀₁⁺|0,0,⅔}, {3₀₀₁⁻|0,0,⅓}
Psτs = [([1, 1, 1],                    [[0.0,0,0], [0,0,2/3], [0,0,1/3]]),
        ([1, cispi(-2/3), cispi(2/3)], [[0.0,0,0], [0,0,2/3], [0,0,1/3]]), 
        ([1, cispi(2/3), cispi(-2/3)], [[0.0,0,0], [0,0,2/3], [0,0,1/3]]), ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 157 =========
sgnum = 157
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z", "x-y,-y,z", "-x,-x+y,z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, m₋₁₁₀, m₁₂₀, m₂₁₀
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops & irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
# ... same lgops & irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 159 =========
sgnum = 159
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z+1/2", "x-y,-y,z+1/2", "-x,-x+y,z+1/2"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, {m₋₁₁₀|0,0,½}, {m₁₂₀|0,0,½}, {m₂₁₀|0,0,½}
Psτs = [[1, 1, 1, -1im, -1im, -1im],
        [1, 1, 1, 1im, 1im, 1im],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 -1; 1 0], [0 cispi(-1/3); cispi(-2/3) 0], [0 cispi(1/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [([1, 1, 1, 1, 1, 1],
         [[0.0,0,0], [0.0,0,0], [0.0,0,0], [0,0,1/2], [0,0,1/2], [0,0,1/2]]),
        ([1, 1, 1, -1, -1, -1],
         [[0.0,0,0], [0.0,0,0], [0.0,0,0], [0,0,1/2], [0,0,1/2], [0,0,1/2]]),
        ([[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
         [[0.0,0,0], [0.0,0,0], [0.0,0,0], [0,0,1/2], [0,0,1/2], [0,0,1/2]]),
        ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 174 =========
sgnum = 174
# HA₁, HA₂, HA₃, HA₄, HA₅, HA₆
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, m₀₀₁, -6₀₀₁⁻, -6₀₀₁⁺
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [1, cispi(-2/3), cispi(2/3), 1, cispi(-2/3), cispi(2/3)],
        [1, cispi(-2/3), cispi(2/3), -1, cispi(1/3), cispi(-1/3)],
        [1, cispi(2/3), cispi(-2/3), 1, cispi(2/3), cispi(-2/3)],
        [1, cispi(2/3), cispi(-2/3), -1, cispi(-1/3), cispi(1/3)], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops and irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻
Psτs = [[1, 1, 1],
        [1, cispi(-2/3), cispi(2/3)],
        [1, cispi(2/3), cispi(-2/3)], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 189 =========
sgnum = 189
# HA₁, HA₂, HA₃, HA₄, HA₅, HA₆
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "x,y,-z", "-y,x-y,-z", "-x+y,-x,-z", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z", "y,x,z", "x-y,-y,z", "-x,-x+y,z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, m₀₀₁, -6₀₀₁⁻, -6₀₀₁⁺, 2₁₁₀, 2₁₀₀, 2₀₁₀, m₋₁₁₀, m₁₂₀, m₂₁₀
Psτs = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1],
        [1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [1 0; 0 1],   [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 1; 1 0],   [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [-1 0; 0 -1], [cispi(1/3) 0; 0 cispi(-1/3)], [cispi(-1/3) 0; 0 cispi(1/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 -1; -1 0], [0 cispi(-1/3); cispi(1/3) 0], [0 cispi(1/3); cispi(-1/3) 0]],
        ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops and irreps as HA₁, HA₂, HA₃
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z", "x-y,-y,z", "-x,-x+y,z"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, m₋₁₁₀, m₁₂₀, m₂₁₀
Psτs = [[1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
        ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 190 =========
sgnum = 190
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec(-1/3,-1/3,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "x,y,-z+1/2", "-y,x-y,-z+1/2", "-x+y,-x,-z+1/2", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z", "y,x,z+1/2", "x-y,-y,z+1/2", "-x,-x+y,z+1/2"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, {m₀₀₁|0,0,½}, {-6₀₀₁⁻|0,0,½}, {-6₀₀₁⁺|0,0,½}, 2₁₁₀, 2₁₀₀, 2₀₁₀, {m₋₁₁₀|0,0,½}, {m₁₂₀|0,0,½}, {m₂₁₀|0,0,½}
Psτs = [[[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [1 0; 0 -1], [cispi(-2/3) 0; 0 cispi(-1/3)], [cispi(2/3) 0; 0 cispi(1/3)],   [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 1; -1 0], [0 cispi(2/3); cispi(1/3) 0],   [0 cispi(-2/3); cispi(-1/3) 0]], 
        [[1 0; 0 1], [cispi(2/3) 0; 0 cispi(-2/3)], [cispi(-2/3) 0; 0 cispi(2/3)], [1 0; 0 -1], [cispi(2/3) 0; 0 cispi(1/3)],   [cispi(-2/3) 0; 0 cispi(-1/3)], [0 1; 1 0], [0 cispi(-2/3); cispi(2/3) 0], [0 cispi(2/3); cispi(-2/3) 0], [0 1; -1 0], [0 cispi(-2/3); cispi(-1/3) 0], [0 cispi(2/3); cispi(1/3) 0]], 
        [[1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 -1], [1 0; 0 -1], [1 0; 0 -1], [0 1; 1 0], [0 1; 1 0], [0 1; 1 0], [0 1; -1 0], [0 1; -1 0], [0 1; -1 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
klab = "KA"
kv   = KVec(-1/3,-1/3,0)
# ... same lgops as HA₁, HA₂, HA₃
Psτs = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1],
        [1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [1 0; 0 1],   [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 1; 1 0],   [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
        [[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [-1 0; 0 -1], [cispi(1/3) 0; 0 cispi(-1/3)], [cispi(-1/3) 0; 0 cispi(1/3)], [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 -1; -1 0], [0 cispi(-1/3); cispi(1/3) 0], [0 cispi(1/3); cispi(-1/3) 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z", "y,x,z+1/2", "x-y,-y,z+1/2", "-x,-x+y,z+1/2"]) # 1, 3₀₀₁⁺, 3₀₀₁⁻, {m₋₁₁₀|0,0,½}, {m₁₂₀|0,0,½}, {m₂₁₀|0,0,½}
Psτs = [([1, 1, 1, 1, 1, 1],    [[0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2]]),
        ([1, 1, 1, -1, -1, -1], [[0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2]]),
        ([[1 0; 0 1], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(2/3) 0; 0 cispi(-2/3)],
          [0 1; 1 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
         [[0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2], [0,0,1/2]]), ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 197 =========
sgnum = 197
# PA₁, PA₂, PA₃, PA₄
klab = "PA"
kv   = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "z,x,y", "z,-x,-y", # 1, 2₀₀₁, 2₀₁₀, 2₁₀₀, 3₁₁₁⁺, 3₋₁₁₋₁⁺, 3₋₁₁₁⁻, 3₋₁₋₁₁⁺, 3₁₁₁⁻, 3₋₁₁₁⁺, 3₋₁₋₁₁⁻, 3₋₁₁₋₁⁻
                "-z,-x,y", "-z,x,-y", "y,z,x", "-y,z,-x", "y,-z,-x", "-y,-z,x"])
Psτs = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, cispi(-2/3), cispi(-2/3), cispi(-2/3), cispi(-2/3), cispi(2/3), cispi(2/3), cispi(2/3), cispi(2/3)],
        [1, 1, 1, 1, cispi(2/3), cispi(2/3), cispi(2/3), cispi(2/3), cispi(-2/3), cispi(-2/3), cispi(-2/3), cispi(-2/3)],
        [[1 0 0; 0 1 0; 0 0 1], [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1],
         [-1 0 0; 0 1 0; 0 0 -1], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; 1 0 0; 0 -1 0],
         [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; -1 0 0; 0 1 0], [0 1 0; 0 0 1; 1 0 0],
         [0 -1 0; 0 0 -1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; -1 0 0]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 199 =========
sgnum = 199
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y+1/2,z", "-x+1/2,y,-z", "x,-y,-z+1/2", "z,x,y", # 1, {2₀₀₁|0,½,0}, {2₀₁₀|½,0,0}, {2₁₀₀|0,0,½}, 3₁₁₁⁺, {3₋₁₁₋₁⁺|0,0,½}, {3₋₁₁₁⁻|0,½,0}, {3₋₁₋₁₁⁺|½,0,0}, 3₁₁₁⁻, {3₋₁₁₁⁺|½,0,0}, {3₋₁₋₁₁⁻|0,0,½}, {3₋₁₁₋₁⁻|0,½,0}
                "z,-x,-y+1/2", "-z,-x+1/2,y", "-z+1/2,x,-y", "y,z,x", "-y+1/2,z,-x", "y,-z,-x+1/2", "-y,-z+1/2,x"])
Psτs = [[[1 0; 0 1], [1 0; 0 -1], [0 1; 1 0], [0 -1im; 1im 0],
         [cispi(-1/12)/√2 cispi(-1/12)/√2; cispi(5/12)/√2 cispi(-7/12)/√2],
         [cispi(-1/12)/√2 cispi(11/12)/√2; cispi(5/12)/√2 cispi(5/12)/√2],
         [cispi(-1/12)/√2 cispi(-1/12)/√2; cispi(-7/12)/√2 cispi(5/12)/√2],
         [cispi(5/12)/√2 cispi(-7/12)/√2; cispi(-1/12)/√2 cispi(-1/12)/√2],
         [cispi(1/12)/√2 cispi(-5/12)/√2; cispi(1/12)/√2 cispi(7/12)/√2],
         [cispi(1/12)/√2 cispi(7/12)/√2; cispi(1/12)/√2 cispi(-5/12)/√2],
         [cispi(-5/12)/√2 cispi(1/12)/√2; cispi(7/12)/√2 cispi(1/12)/√2],
         [cispi(1/12)/√2 cispi(-5/12)/√2; cispi(-11/12)/√2 cispi(-5/12)/√2]],
        [[1 0; 0 1], [1 0; 0 -1], [0 1; 1 0], [0 -1im; 1im 0], 
         [cispi(-3/4)/√2 cispi(-3/4)/√2; cispi(-1/4)/√2 cispi(3/4)/√2],
         [cispi(-3/4)/√2 cispi(1/4)/√2; cispi(-1/4)/√2 cispi(-1/4)/√2],
         [cispi(-3/4)/√2 cispi(-3/4)/√2; cispi(3/4)/√2 cispi(-1/4)/√2],
         [cispi(-1/4)/√2 cispi(3/4)/√2; cispi(-3/4)/√2 cispi(-3/4)/√2],
         [cispi(3/4)/√2 cispi(1/4)/√2; cispi(3/4)/√2 cispi(-3/4)/√2],
         [cispi(3/4)/√2 cispi(-3/4)/√2; cispi(3/4)/√2 cispi(1/4)/√2],
         [cispi(1/4)/√2 cispi(3/4)/√2; cispi(-3/4)/√2 cispi(3/4)/√2],
         [cispi(3/4)/√2 cispi(1/4)/√2; cispi(-1/4)/√2 cispi(1/4)/√2]],
        [[1 0; 0 1], [1 0; 0 -1], [0 1; 1 0], [0 -1im; 1im 0], 
         [cispi(7/12)/√2 cispi(7/12)/√2; cispi(-11/12)/√2 cispi(1/12)/√2],
         [cispi(7/12)/√2 cispi(-5/12)/√2; cispi(-11/12)/√2 cispi(-11/12)/√2],
         [cispi(7/12)/√2 cispi(7/12)/√2; cispi(1/12)/√2 cispi(-11/12)/√2],
         [cispi(-11/12)/√2 cispi(1/12)/√2; cispi(7/12)/√2 cispi(7/12)/√2],
         [cispi(-7/12)/√2 cispi(11/12)/√2; cispi(-7/12)/√2 cispi(-1/12)/√2],
         [cispi(-7/12)/√2 cispi(-1/12)/√2; cispi(-7/12)/√2 cispi(11/12)/√2],
         [cispi(11/12)/√2 cispi(-7/12)/√2; cispi(-1/12)/√2 cispi(-7/12)/√2],
         [cispi(-7/12)/√2 cispi(11/12)/√2; cispi(5/12)/√2 cispi(11/12)/√2]], ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 217 =========
sgnum = 217
# PA₁, PA₂, PA₃, PA₄, PA₅
klab = "PA"
kv   = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "z,x,y", "z,-x,-y", # 1, 2₀₀₁, 2₀₁₀, 2₁₀₀, 3₁₁₁⁺, 3₋₁₁₋₁⁺, 3₋₁₁₁⁻, 3₋₁₋₁₁⁺, 3₁₁₁⁻, 3₋₁₁₁⁺, 3₋₁₋₁₁⁻, 3₋₁₁₋₁⁻, m₋₁₁₀, m₁₁₀, -4₀₀₁⁺, -4₀₀₁⁻, m₀₋₁₁, -4₁₀₀⁺, -4₁₀₀⁻, m₀₁₁, m₋₁₀₁, -4₀₁₀⁻, m₁₀₁, -4₀₁₀⁺
                "-z,-x,y", "-z,x,-y", "y,z,x", "-y,z,-x", "y,-z,-x", "-y,-z,x", "y,x,z",
                "-y,-x,z", "y,-x,-z", "-y,x,-z", "x,z,y", "-x,z,-y", "-x,-z,y", "x,-z,-y",
                "z,y,x", "z,-y,-x", "-z,y,-x", "-z,-y,x"])
Psτs = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [[1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1],
         [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(-2/3) 0; 0 cispi(2/3)], [cispi(-2/3) 0; 0 cispi(2/3)],
         [cispi(2/3) 0; 0 cispi(-2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [cispi(2/3) 0; 0 cispi(-2/3)], [cispi(2/3) 0; 0 cispi(-2/3)],
         [0 1; 1 0], [0 1; 1 0], [0 1; 1 0], [0 1; 1 0],
         [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(2/3); cispi(-2/3) 0], [0 cispi(2/3); cispi(-2/3) 0],
         [0 cispi(-2/3); cispi(2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 cispi(-2/3); cispi(2/3) 0], [0 cispi(-2/3); cispi(2/3) 0]],
        [[1 0 0; 0 1 0; 0 0 1], [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; 1 0 0; 0 -1 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; -1 0 0; 0 1 0], [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 -1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; -1 0 0], 
         [1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 0 -1; 0 -1 0], [-1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 1 0], [0 0 1; 0 1 0; 1 0 0], [0 0 -1; 0 -1 0; 1 0 0], [0 0 1; 0 -1 0; -1 0 0], [0 0 -1; 0 1 0; -1 0 0], [0 1 0; 1 0 0; 0 0 1], [0 -1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 1], [0 1 0; -1 0 0; 0 0 -1]],
        [[1 0 0; 0 1 0; 0 0 1], [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; 1 0 0; 0 -1 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; -1 0 0; 0 1 0], [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 -1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; -1 0 0],              # same as first "half" of PA₄
         -[1 0 0; 0 0 1; 0 1 0], -[1 0 0; 0 0 -1; 0 -1 0], -[-1 0 0; 0 0 1; 0 -1 0], -[-1 0 0; 0 0 -1; 0 1 0], -[0 0 1; 0 1 0; 1 0 0], -[0 0 -1; 0 -1 0; 1 0 0], -[0 0 1; 0 -1 0; -1 0 0], -[0 0 -1; 0 1 0; -1 0 0], -[0 1 0; 1 0 0; 0 0 1], -[0 -1 0; 1 0 0; 0 0 -1], -[0 -1 0; -1 0 0; 0 0 1], -[0 1 0; -1 0 0; 0 0 -1]], # negative of second "half" of PA₄
        ]
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)

# ========= 220 =========
sgnum = 220
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y+1/2,z", "-x+1/2,y,-z", "x,-y,-z+1/2", "z,x,y", # 1, {2₀₀₁|0,½,0}, {2₀₁₀|½,0,0}, {2₁₀₀|0,0,½}, 3₁₁₁⁺, {3₋₁₁₋₁⁺|0,0,½}, {3₋₁₁₁⁻|0,½,0}, {3₋₁₋₁₁⁺|½,0,0}, 3₁₁₁⁻, {3₋₁₁₁⁺|½,0,0}, {3₋₁₋₁₁⁻|0,0,½}, {3₋₁₁₋₁⁻|0,½,0}, {m₋₁₁₀|¼,¼,¼}, {m₁₁₀|¾,¼,¼}, {-4₀₀₁⁺|¼,¾,¼}, {-4₀₀₁⁻|¼,¼,¾}, {m₀₋₁₁|¼,¼,¼}, {-4₁₀₀⁺|¼,¼,¾}, {-4₁₀₀⁻|¾,¼,¼}, {m₀₁₁|¼,¾,¼}, {m₋₁₀₁|¼,¼,¼}, {-4₀₁₀⁻|¼,¾,¼}, {m₁₀₁|¼,¼,¾}, {-4₀₁₀⁺|¾,¼,¼}
                "z,-x,-y+1/2", "-z,-x+1/2,y", "-z+1/2,x,-y", "y,z,x", "-y+1/2,z,-x",
                "y,-z,-x+1/2", "-y,-z+1/2,x", "y+1/4,x+1/4,z+1/4", "-y+3/4,-x+1/4,z+1/4",
                "y+1/4,-x+3/4,-z+1/4", "-y+1/4,x+1/4,-z+3/4", "x+1/4,z+1/4,y+1/4",
                "-x+1/4,z+1/4,-y+3/4", "-x+3/4,-z+1/4,y+1/4", "x+1/4,-z+3/4,-y+1/4",
                "z+1/4,y+1/4,x+1/4", "z+1/4,-y+3/4,-x+1/4", "-z+1/4,y+1/4,-x+3/4",
                "-z+3/4,-y+1/4,x+1/4"])
Psτs =      # PA₁
    [[ [1 0; 0 1], [1 0; 0 -1], [0 1; 1 0], [0 -1im; 1im 0],
       [cispi(-3/4)/√2 cispi(-3/4)/√2; cispi(-1/4)/√2 cispi(3/4)/√2 ],
       [cispi(-3/4)/√2 cispi(1/4)/√2;  cispi(-1/4)/√2 cispi(-1/4)/√2],
       [cispi(-3/4)/√2 cispi(-3/4)/√2; cispi(3/4)/√2  cispi(-1/4)/√2],
       [cispi(-1/4)/√2 cispi(3/4)/√2;  cispi(-3/4)/√2 cispi(-3/4)/√2],
       [cispi(3/4)/√2  cispi(1/4)/√2;  cispi(3/4)/√2  cispi(-3/4)/√2],
       [cispi(3/4)/√2  cispi(-3/4)/√2; cispi(3/4)/√2  cispi(1/4)/√2 ], 
       [cispi(1/4)/√2  cispi(3/4)/√2;  cispi(-3/4)/√2 cispi(3/4)/√2 ],
       [cispi(3/4)/√2  cispi(1/4)/√2;  cispi(-1/4)/√2 cispi(1/4)/√2 ],
       [0 -1im; -1 0], [0 1im; -1 0], [-1im 0; 0 -1], [1 0; 0 1im],
       [cispi(-3/4)/√2 cispi(1/4)/√2;  cispi(1/4)/√2  cispi(1/4)/√2 ],
       [cispi(-3/4)/√2 cispi(-3/4)/√2; cispi(1/4)/√2  cispi(-3/4)/√2],
       [cispi(1/4)/√2  cispi(-3/4)/√2; cispi(1/4)/√2  cispi(1/4)/√2 ],
       [cispi(3/4)/√2  cispi(3/4)/√2;  cispi(3/4)/√2  cispi(-1/4)/√2],
       [cispi(1/4)/√2  cispi(3/4)/√2;  cispi(-1/4)/√2 cispi(-3/4)/√2],
       [cispi(1/4)/√2  cispi(-1/4)/√2; cispi(-1/4)/√2 cispi(1/4)/√2 ],
       [cispi(3/4)/√2  cispi(1/4)/√2;  cispi(-3/4)/√2 cispi(-1/4)/√2],
       [cispi(-3/4)/√2 cispi(-1/4)/√2; cispi(-1/4)/√2 cispi(-3/4)/√2], ]]
push!(Psτs, # PA₂ (first 12 irreps same as PA₁, next 12 have sign flipped)
        Psτs[1] .* vcat(ones(Int, 12), fill(-1, 12))
     )
push!(Psτs, # PA₃ (4-dimensional irrep)
     [ [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],        # ---
       [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1],      # ---
       [0 1 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0],   # ---
       [0 -1im 0 0; 1im 0 0 0; 0 0 0 -1; 0 0 -1 0], # ---
       [cispi(7/12)/√2 cispi(7/12)/√2 0 0;          # ---
        cispi(-11/12)/√2 cispi(1/12)/√2 0 0;
        0 0 cispi(-7/12)/√2 cispi(11/12)/√2;
        0 0 cispi(-7/12)/√2 cispi(-1/12)/√2],
       [cispi(7/12)/√2 cispi(-5/12)/√2 0 0;         # ---
        cispi(-11/12)/√2 cispi(-11/12)/√2 0 0;
        0 0 cispi(5/12)/√2 cispi(11/12)/√2;
        0 0 cispi(5/12)/√2 cispi(-1/12)/√2],
       [cispi(7/12)/√2 cispi(7/12)/√2 0 0;          # ---
        cispi(1/12)/√2 cispi(-11/12)/√2 0 0;
        0 0 cispi(5/12)/√2 cispi(-1/12)/√2;
        0 0 cispi(-7/12)/√2 cispi(-1/12)/√2],
       [cispi(-11/12)/√2 cispi(1/12)/√2 0 0;        # ---
        cispi(7/12)/√2 cispi(7/12)/√2 0 0;
        0 0 cispi(-1/12)/√2 cispi(5/12)/√2;
        0 0 cispi(11/12)/√2 cispi(5/12)/√2],
       [cispi(-7/12)/√2 cispi(11/12)/√2 0 0;        # ---
        cispi(-7/12)/√2 cispi(-1/12)/√2 0 0;
        0 0 cispi(7/12)/√2 cispi(7/12)/√2;
        0 0 cispi(-11/12)/√2 cispi(1/12)/√2],
       [cispi(-7/12)/√2 cispi(-1/12)/√2 0 0;        # ---
        cispi(-7/12)/√2 cispi(11/12)/√2 0 0;
        0 0 cispi(-5/12)/√2 cispi(7/12)/√2;
        0 0 cispi(1/12)/√2 cispi(1/12)/√2],
       [cispi(11/12)/√2 cispi(-7/12)/√2 0 0;        # ---
        cispi(-1/12)/√2 cispi(-7/12)/√2 0 0;
        0 0 cispi(1/12)/√2 cispi(-11/12)/√2;
        0 0 cispi(-5/12)/√2 cispi(-5/12)/√2],
       [cispi(-7/12)/√2 cispi(11/12)/√2 0 0;        # ---
        cispi(5/12)/√2 cispi(11/12)/√2 0 0;
        0 0 cispi(-5/12)/√2 cispi(-5/12)/√2;
        0 0 cispi(-11/12)/√2 cispi(1/12)/√2],
       [0 0 1im 0; 0 0 0 1im; 1 0 0 0; 0 1 0 0],    # ---
       [0 0 -1im 0; 0 0 0 1im; 1 0 0 0; 0 -1 0 0],  # ---
       [0 0 0 -1; 0 0 1 0; 0 1 0 0; 1 0 0 0],       # ---
       [0 0 0 -1im; 0 0 -1im 0; 0 -1im 0 0; 1im 0 0 0], # ---
       [0 0 cispi(-1/12)/√2 cispi(-7/12)/√2;        # ---
        0 0 cispi(-1/12)/√2 cispi(5/12)/√2;
        cispi(7/12)/√2 cispi(7/12)/√2 0 0;
        cispi(-11/12)/√2 cispi(1/12)/√2 0 0],
       [0 0 cispi(11/12)/√2 cispi(-7/12)/√2;        # ---
        0 0 cispi(11/12)/√2 cispi(5/12)/√2;
        cispi(7/12)/√2 cispi(-5/12)/√2 0 0;
        cispi(-11/12)/√2 cispi(-11/12)/√2 0 0],
       [0 0 cispi(11/12)/√2 cispi(5/12)/√2;         # ---
        0 0 cispi(-1/12)/√2 cispi(5/12)/√2;
        cispi(7/12)/√2 cispi(7/12)/√2 0 0;
        cispi(1/12)/√2 cispi(-11/12)/√2 0 0],
       [0 0 cispi(5/12)/√2 cispi(11/12)/√2;         # ---
        0 0 cispi(-7/12)/√2 cispi(11/12)/√2;
        cispi(-11/12)/√2 cispi(1/12)/√2 0 0;
        cispi(7/12)/√2 cispi(7/12)/√2 0 0],
       [0 0 cispi(-11/12)/√2 cispi(-11/12)/√2;      # ---
        0 0 cispi(-5/12)/√2 cispi(7/12)/√2;
        cispi(-7/12)/√2 cispi(11/12)/√2 0 0;
        cispi(-7/12)/√2 cispi(-1/12)/√2 0 0],
       [0 0 cispi(1/12)/√2 cispi(-11/12)/√2;        # ---
        0 0 cispi(7/12)/√2 cispi(7/12)/√2;
        cispi(-7/12)/√2 cispi(-1/12)/√2 0 0;
        cispi(-7/12)/√2 cispi(11/12)/√2 0 0],
       [0 0 cispi(7/12)/√2 cispi(-5/12)/√2;         # ---
        0 0 cispi(1/12)/√2 cispi(1/12)/√2;
        cispi(11/12)/√2 cispi(-7/12)/√2 0 0;
        cispi(-1/12)/√2 cispi(-7/12)/√2 0 0],
       [0 0 cispi(1/12)/√2 cispi(1/12)/√2;          # ---
        0 0 cispi(-5/12)/√2 cispi(7/12)/√2;
        cispi(-7/12)/√2 cispi(11/12)/√2 0 0;
        cispi(5/12)/√2 cispi(11/12)/√2 0 0],
     ]
  )
LGIRS_add[sgnum][klab] = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)