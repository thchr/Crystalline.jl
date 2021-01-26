using Crystalline

# The listings below include all unique, labelled k-points in the representation domain
# whose little group is not trivial. In addition, Γ=(0,0) and Ω=(u,v) are included for 
# every group. The notation follows of layer groups as given by LKVEC on the Bilbao server
# (https://www.cryst.ehu.es/subperiodic/get_layer_kvec.html).
PLANE2KVEC = Dict(
    1  => ([],               KVec{2}.([])),
    2  => (["Y", "B", "A"],  KVec{2}.(["0,1/2",   "1/2,0",   "1/2,-1/2"])),
    3  => (["X", "Y", "S", 
            "Σ", "ΣA", "C",
            "CA"],           KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "0,u",   "0,-u",  "1/2,u",    "1/2,-u"])),
    4  => (["X", "Y", "S",   #= non-symmorphic =#
            "Σ", "ΣA", "C",
            "CA"],           KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "0,u",   "0,-u",  "1/2,u",    "1/2,-u"])),
    5  => (["Y", "Σ", "ΣA",
            "C", "CA"],      KVec{2}.(["0,1",     "0,2u",    "0,-2u",   "1,2u",  "1,-2u"])),
    6  => (["X", "S", "Y",
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["1/2,0",   "1/2,1/2", "0,1/2",   "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    7  => (["Y", "X", "S",   #= non-symmorphic =#
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    8  => (["Y", "X", "S",   #= non-symmorphic =#
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    9  => (["Y", "S", "Σ",
            "Δ", "F", "C"],  KVec{2}.(["1,0",     "1/2,1/2", "2u,0",    "0,2u",  "1,2u",  "2u,1"])),
    10 => (["X", "M"],       KVec{2}.(["0,1/2",   "1/2,1/2"])),
    11 => (["X", "M", "Σ",
            "Δ", "Y"],       KVec{2}.(["0,1/2",   "1/2,1/2", "u,u",     "0,u",   "u,1/2"])),
    12 => (["M", "X", "Δ",   #= non-symmorphic =#
            "Y", "Σ"],       KVec{2}.(["1/2,1/2", "0,1/2",   "0,u",     "u,1/2", "u,u"])),
    13 => (["K", "KA"],      KVec{2}.(["1/3,1/3", "2/3,-1/3"])),
    14 => (["K", "M", "Σ",
            "ΣA"],           KVec{2}.(["1/3,1/3", "1/2,0",   "u,0",   "0,u"])),
    15 => (["K", "KA", "M",
            "Λ", "ΛA", "T",
            "TA"],           KVec{2}.(["1/3,1/3", "2/3,-1/3", "1/2,0",  "u,u",   "2u,-u", "1/2-u,2u", "1/2+u,-2u"])),
    16 => (["K", "M"],       KVec{2}.(["1/3,1/3", "1/2,0"])),
    17 => (["K", "M", "Σ",
            "Λ", "T"],       KVec{2}.(["1/3,1/3", "1/2,0",   "u,0",   "u,u", "1/2-u,2u"]))
)
# push Γ and Ω k-points, which are in all plane groups, to front & end of the above vectors
foreach(values(PLANE2KVEC)) do d
    pushfirst!(d[1], "Γ");         push!(d[1], "Ω")
    pushfirst!(d[2], KVec("0,0")); push!(d[2], KVec("α,β"))
end

# should be called with a reduced set of operations (i.e. no centering translation copies)
function _find_isomorphic_parent_pointgroup(G)
    D = dim(G)
    ctᴳ = MultTable(G).table
    @inbounds for iuclab in PGS_IUCs[D]
        P = pointgroup(iuclab, D)
        ctᴾ = MultTable(P).table
        if ctᴳ == ctᴾ # bit sloppy; would miss ismorphisms that are "concealed" by row/column swaps
            return P
        end
    end
    return nothing # in case we didn't find any isomorphic parent
end

# build little group irreps of tabulated k-points by matching up to the assoc. point group
LGIRSD_2D = Dict{Int, Dict{String, Vector{LGIrrep{2}}}}()
for (sgnum, (klabs, kvs)) in PLANE2KVEC
    issymmorph(sgnum, 2) || continue  # treat only symmorphic groups
    
    sg   = spacegroup(sgnum, Val(2))
    lgs  = littlegroup.(Ref(sg), kvs) # little groups at each tabulated k-point
    ops′ = reduce_ops.(lgs)           # reduce to operations without centering "copies"
    lgs .= LittleGroup.(sgnum, kvs, klabs, ops′)
    pgs  = _find_isomorphic_parent_pointgroup.(lgs) # find the parent point group 

    pglabs    = label.(pgs) # labels of the parent point groups
    pgirs_vec = get_pgirreps.(pglabs, Val(2))
    
    LGIRSD_2D[sgnum] = Dict{String, Vector{LGIrrep{2}}}()
    for (klab, kvec, lg, pgirs) in zip(klabs, kvs, lgs, pgirs_vec)
        # regarding cdml labels: in 3D some CDML labels can have a ± superscript; in 2D, 
        # however, the labelling scheme is just numbers; taking `last(label(pgir))` picks
        # out this number ('₁', '₂', etc.)
        cdmls        = Ref(klab).*last.(label.(pgirs))
        matrices     = getfield.(pgirs, Ref(:matrices))
        translations = [zeros(2) for _ in lg]
        pgirs_realities = reality.(pgirs)

        lgirs = LGIrrep{2}.(cdmls, Ref(lg), matrices, Ref(translations), pgirs_realities, false)
        LGIRSD_2D[sgnum][klab] = lgirs
    end
end

# correct the reality by calling out to `calc_reality` explicitly (the thing to guard
# against here is that the Herring criterion and the Frobenius-Schur criterion need not 
# agree: i.e. LGIrreps do not necessarily inherit the reality of a parent PGIrrep)
LGIRSD_2D′ = Dict{Int, Dict{String, Vector{LGIrrep{2}}}}()
for (sgnum, LGIRSD_2D) in LGIRSD_2D
    sg = reduce_ops(spacegroup(sgnum, Val(2)))

    LGIRSD_2D′[sgnum] = Dict{String, Vector{LGIrrep{2}}}()
    for (klab, lgirs) in LGIRSD_2D
        lgirs′ = LGIrrep{2}[]
        for lgir in lgirs
            reality_type  = calc_reality(lgir, sg)
            lgir′ = LGIrrep{2}(label(lgir), group(lgir), lgir.matrices, lgir.translations,
                               reality_type , false)
            push!(lgirs′, lgir′)        
        end
        LGIRSD_2D′[sgnum][klab]= lgirs′
    end
end

# -----------------------------------------------------------------------------------------
# ADD LGIRREPS OR THE NON-SYMMORPHIC PLANE GROUPS (FROM HAND-TABULATION)
let
    # "load" the `LGIrrep`s of plane groups 4, 7, 8, and 12 into `LGIRSD_2D` (using `let`
    # statement to avoid namespace clash w/ existing `LGIRSD_2D`)
    include("setup_2d_littlegroup_irreps_nonsymmorph.jl") 
    merge!(LGIRSD_2D′, LGIRSD_2D)
end

# -----------------------------------------------------------------------------------------
# CREATE A SORTED VECTOR OF `Vector{LGIrrep}` WHOSE INDICES MATCH THE SGNUM
LGIRS_2D′ = Vector{valtype(LGIRSD_2D′)}(undef, 17)
foreach(LGIRSD_2D′) do (sgnum, lgirsd)
    LGIRS_2D′[sgnum] = lgirsd
end