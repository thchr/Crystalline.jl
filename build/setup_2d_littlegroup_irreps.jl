using Crystalline

# NOTE/TODO: Currently only deals with symmorphic groups
PLANE2KVEC = Dict(
    1  => ([],               KVec{2}.([])),
    2  => (["Y", "B", "A"],  KVec{2}.(["0,1/2",    "1/2,0",   "1/2,-1/2"])),
    3  => (["Y", "X", "S"],  KVec{2}.(["0,1/2",    "1/2,0",   "1/2,1/2" ])),
    5  => (["Y", "S"],       KVec{2}.(["1/2,1/2",  "0,1/2"  ])),
    6  => (["X", "S", "Y"],  KVec{2}.(["1/2,0",    "1/2,1/2", "0,1/2"])),
    9  => (["Y", "S"],       KVec{2}.(["1,0",      "1/2,1/2"])),
    10 => (["X", "M"],       KVec{2}.(["0,1/2",    "1/2,1/2"])),
    11 => (["X", "M" ],      KVec{2}.(["0,1/2",    "1/2,1/2"])),
    13 => (["KA", "K"],      KVec{2}.(["2/3,-1/3", "1/3,1/3"])),
    14 => (["K", "M"],       KVec{2}.(["1/3,1/3",  "1/2,0"  ])),
    15 => (["KA", "K", "M"], KVec{2}.(["2/3,-1/3", "1/3,1/3", "1/2,0"])),
    16 => (["K", "M"],       KVec{2}.(["1/3,1/3",  "1/2,0"  ])),
    17 => (["K", "M"],       KVec{2}.(["1/3,1/3",  "1/2,0"  ]))
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
        type         = -1 # sentinel value for undetermined reality type; fixed below
        pgirs_types  = getfield.(pgirs, Ref(:type))

        lgirs = LGIrrep{2}.(cdmls, Ref(lg), matrices, Ref(translations), pgirs_types, false)
        LGIRSD_2D[sgnum][klab] = lgirs
    end
end

# correct the reality type by calling out to herring explicitly
LGIRSD_2D′ = Dict{Int, Dict{String, Vector{LGIrrep{2}}}}()
for (sgnum, LGIRSD_2D) in LGIRSD_2D
    sg    = reduce_ops(spacegroup(sgnum, Val(2)))
    LGIRSD_2D′[sgnum] = Dict{String, Vector{LGIrrep{2}}}()

    for (klab, lgirs) in LGIRSD_2D
        lgirs′ = LGIrrep{2}[]
        for lgir in lgirs
            lgir′ = LGIrrep{2}(label(lgir), group(lgir), lgir.matrices, lgir.translations,
                               herring(lgir, sg), false)
            push!(lgirs′, lgir′)        
        end
        LGIRSD_2D′[sgnum][klab]= lgirs
    end
end
LGIRS_2D′ = [LGIRSD_2D for (key, LGIRSD_2D) in LGIRSD_2D′]