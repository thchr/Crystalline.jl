using SGOps

# manually converting
#    https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_out.pl
# to our own data format for the missing k-points listed in 
#   src/special_representation_domain_kpoints.jl

function build_lgirrep_with_type(cdml, lg, Psτs, sgops)
    lgir  = LGIrrep{3}(cdml, lg, Psτs[1], Psτs[2], 0) # sentinel =0 for reality type
    typeᴴ = herring(lgir, sgops)
    type  = typeᴴ == -1 ? 2 : typeᴴ == 0 ? 3 : 1 # {1,-1,0} ⇒ {1,2,3} (Herring ⇒ ISOTROPY)
    lgir  = LGIrrep{3}(cdml, lg, Psτs[1], Psτs[2], type) # update reality type
end

function prepare_lg_and_sgops(sgnum, kv, klab, ops)
    lg    = LittleGroup{3}(sgnum, kv, klab, ops)
    sgops = reduce_ops(spacegroup(sgnum, Val(3)), centering(sgnum, 3), true) # for herrring
    return lg, sgops
end

function allocate_containers(Nirr)
    return Vector{Any}(undef, Nirr),     # matrices (Ps)
           Vector{Any}(undef, Nirr),     # translations (τs)
           0                             # index
end

function assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)
    lg, sgops = prepare_lg_and_sgops(sgnum, kv, klab, lgops)
    cdmls     = Ref(klab) .* string.(1:length(Psτs))
    lgirs     = build_lgirrep_with_type.(cdmls, Ref(lg), Psτs, Ref(sgops))
end

const 𝗶 = fill(ComplexF64(1,0),1,1)
const 𝗼 = zeros(Float64, 3)
const C64 = ComplexF64


# ========= 23 =========
sgnum = 23
klab  = "WA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "x,-y,-z", "-x,y,-z", "-x,-y,z"]) # 1, 2₁₀₀, 2₀₁₀, 2₀₀₁

# sorted in ascending irrep order (e.g. WA1, WA2, WA3, WA4 here)
# ----------------------
Psτs = [([𝗶, 𝗶, 𝗶, 𝗶],      # Ps = matrices     (for WA1, across lgops)
         [𝗼, 𝗼, 𝗼, 𝗼]),   # τs = translations
# ----------------------
        ([𝗶, -𝗶, -𝗶, 𝗶],
         [𝗼, 𝗼, 𝗼, 𝗼]),
# ----------------------
        ([𝗶, 𝗶, -𝗶, -𝗶],
        [𝗼, 𝗼, 𝗼, 𝗼]),
# ----------------------
        ([𝗶, -𝗶, 𝗶, -𝗶],
        [𝗼, 𝗼, 𝗼, 𝗼])]
# ----------------------
lgirs23 = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 24 =========
sgnum = 24
klab  = "WA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "x,-y,-z+1/2", "-x+1/2,y,-z", "-x,-y+1/2,z"]) # 1, {2₁₀₀|00½}, {2₀₁₀|½00}, {2₀₀₁|0½0}

# ----------------------
Psτs = [([C64.([1 0; 0 1]), C64.([1 0; 0 -1]), C64.([0 -im; im 0]), C64.([0 1; 1 0])],
         [𝗼, 𝗼, 𝗼, 𝗼]),]
# ----------------------
lgirs24 = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 82  =========
sgnum = 82
klab  = "PA"
kv    = KVec(-1/2,-1/2,-1/2)
lgops = SymOperation{3}.(["x,y,z", "-x,-y,z", "y,-x,-z", "-y,x,-z"]) # 1, 2₀₀₁, -4⁺₀₀₁, -4⁻₀₀₁

Ps, τs, i = allocate_containers(Nirr)
# ----------------------
Psτs = [([𝗶, 𝗶, 𝗶, 𝗶],
         [𝗼, 𝗼, 𝗼, 𝗼]),
# ----------------------
        ([𝗶, 𝗶, -𝗶, -𝗶],
         [𝗼, 𝗼, 𝗼, 𝗼]),
# ----------------------
        ([𝗶, -𝗶, -𝗶, 𝗶],
         [𝗼, 𝗼, 𝗼, 𝗼]),
# ----------------------
        ([𝗶, -𝗶, 𝗶, -𝗶],
         [𝗼, 𝗼, 𝗼, 𝗼])]
# ----------------------
lgirs82 = assemble_lgirreps(sgnum, kv, klab, lgops, Psτs)


# ========= 121 =========
# ========= 122 =========
# ========= 143 =========
# ========= 144 =========
# ========= 145 =========
# ========= 150 =========
# ========= 152 =========
# ========= 154 =========
# ========= 157 =========
# ========= 159 =========
# ========= 174 =========
# ========= 189 =========
# ========= 190 =========
# ========= 197 =========
# ========= 199 =========
# ========= 217 =========
# ========= 220 =========