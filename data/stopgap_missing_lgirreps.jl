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

    lgir  = LGIrrep{3}(cdml, lg, Ps, τs, 0) # sentinel =0 for reality type
    typeᴴ = herring(lgir, sgops)
    type  = typeᴴ == -1 ? 2 : typeᴴ == 0 ? 3 : 1 # {1,-1,0} ⇒ {1,2,3} (Herring ⇒ ISOTROPY)
    lgir  = LGIrrep{3}(cdml, lg, Ps, τs, type) # update reality type
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

# "preallocate" a storage dict
sgnums = [23, 24, 82, 121, 122, 143, 144, 145, 150, 152, 154, 157, 159, 174, 189, 190, 197, 
          199, 217, 220]
lgirs_dict = Dict(sgnum=>Vector{Vector{LGIrrep{3}}}() for sgnum in sgnums)

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
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))


# ========= 24 =========
sgnum = 24
# WA₁
klab  = "WA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "x,-y,-z+1/2", "-x+1/2,y,-z", "-x,-y+1/2,z"]) # 1, {2₁₀₀|00½}, {2₀₁₀|½00}, {2₀₀₁|0½0}

Psτs = [[[1 0; 0 1], [1 0; 0 -1], [0 -im; im 0], [0 1; 1 0]],
       ]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))


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
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))


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
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# ========= 122 =========
sgnum = 122
# PA₁, PA₂
klab  = "PA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z", "-x,y+1/2,-z+1/4", "x,-y+1/2,-z+1/4", "-y,-x+1/2,z+1/4", "y,x+1/2,z+1/4"]) 
                          # 1, 2₀₀₁, -4⁺₀₀₁, -4⁻₀₀₁, {2₀₁₀|0,1/2,1/4}, {2₁₀₀|0,1/2,1/4}, {m₁₁₀|0,1/2,1/4}, {m₁₋₁₀|0,1/2,1/4}

Psτs = [[[1 0; 0 1], [1 0; 0 -1], [1 0; 0 im], [1 0; 0 -im], [0 -1; 1 0], [0 1; 1 0], [0 -im; 1 0], [0 im; 1 0]],
        [[1 0; 0 1], [1 0; 0 -1], [-1 0; 0 -im], [-1 0; 0 im], [0 1; -1 0], [0 -1; -1 0], [0 -im; 1 0], [0 im; 1 0]],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# ========= 143 =========
sgnum = 143
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z", "-x+y,-x,z"]) # 1, 3⁺₀₀₁, 3⁻₀₀₁
# HA₁, HA₂, HA₃
klab  = "HA"
kv    = KVec("-1/3,-1/3,-1/2")

Psτs = [[1, 1, 1],
        [1, cis(-2π/3), cis(2π/3)],
        [1, cis(2π/3), cis(-2π/3)],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# ... same lgops & irreps as HA₁, HA₂, HA₃
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
# ... same lgops & irreps as HA₁, HA₂, HA₃
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# ========= 144 =========
sgnum = 144
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]) # 1, {3⁺₀₀₁|0,0,1/3}, {3⁻₀₀₁|0,0,2/3}
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")

Psτs = [[1, cis(-π/3), cis(-2π/3)],
        [1, -1, 1],
        [1, cis(π/3), cis(-2π/3)],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
Psτs = [[1, 1, 1],
        [1, cis(-2π/3), cis(2π/3)],
        [1, cis(2π/3), cis(-2π/3)],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
Psτs = [([1, 1, 1],                  [[0.0, 0, 0], [0,0,1/3], [0,0,2/3]]), # nonzero τs
        ([1, cis(-2π/3), cis(2π/3)], [[0.0, 0, 0], [0,0,1/3], [0,0,2/3]]),
        ([1, cis(2π/3), cis(-2π/3)], [[0.0, 0, 0], [0,0,1/3], [0,0,2/3]]),]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# ========= 145 =========
sgnum = 145
lgops = SymOperation{3}.(["x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3"]) # 1, {3⁺₀₀₁|0,0,2/3}, {3⁻₀₀₁|0,0,1/3}
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
Psτs = [[1, 1, -1],
        [1, cis(-2π/3), cis(-π/3)],
        [1, cis(2π/3), cis(π/3)],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
Psτs = [[1, 1, 1],
        [1, cis(-2π/3), cis(2π/3)],
        [1, cis(2π/3), cis(-2π/3)],]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")
Psτs = [([1, 1, 1],                  [[0.0, 0, 0], [0,0,2/3], [0,0,1/3]]), # nonzero τs
        ([1, cis(-2π/3), cis(2π/3)], [[0.0, 0, 0], [0,0,2/3], [0,0,1/3]]),
        ([1, cis(2π/3), cis(-2π/3)], [[0.0, 0, 0], [0,0,2/3], [0,0,1/3]]),]
push!(lgirs_dict[sgnum], assemble_lgirreps(sgnum, kv, klab, lgops, Psτs))

# ========= 150 =========
sgnum = 150
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")

# ========= 152 =========
sgnum = 152
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")

# ========= 154 =========
sgnum = 154
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")

# ========= 157 =========
sgnum = 157
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")

# ========= 159 =========
sgnum = 159
# HA₁, HA₂, HA₃
klab = "HA"
kv   = KVec("-1/3,-1/3,-1/2")
# KA₁, KA₂, KA₃
klab = "KA"
kv   = KVec("-1/3,-1/3,0")
# PA₁, PA₂, PA₃
klab = "PA"
kv   = KVec("-1/3,-1/3,-w")

# ========= 174 =========
sgnum = 174
# HA₁, HA₂, HA₃, HA₄, HA₅, HA₆
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
# PA₁, PA₂, PA₃

# ========= 189 =========
sgnum = 189
# HA₁, HA₂, HA₃, HA₄, HA₅, HA₆
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
# PA₁, PA₂, PA₃

# ========= 190 =========
sgnum = 190
# HA₁, HA₂, HA₃
# KA₁, KA₂, KA₃, KA₄, KA₅, KA₆
# PA₁, PA₂, PA₃

# ========= 197 =========
sgnum = 197
# PA₁, PA₂, PA₃, PA₄

# ========= 199 =========
sgnum = 199
# PA₁, PA₂, PA₃

# ========= 217 =========
sgnum = 217
# PA₁, PA₂, PA₃, PA₄, PA₅

# ========= 220 =========
sgnum = 220
# PA₁, PA₂, PA₃
