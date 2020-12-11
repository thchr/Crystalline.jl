using Crystalline

# NOTE/TODO: Currently only deals with symmorphic groups
PLANE2KVEC = Dict{String, Any}(
    "p1"   => (["Γ", "Ω"],            KVec.(["0,0", "α,β"])),
    "p2"   => (["Γ", "Y", "B", "A"],  KVec.(["0,0", "0,1/2",    "1/2,0",    "1/2,-1/2"])),
    "p1m1" => (["Γ", "Y", "X", "S"],  KVec.(["0,0", "0,1/2",    "1/2,0",    "1/2,1/2" ])),
    "c1m1" => (["Γ", "Y", "S"],       KVec.(["0,0", "1/2,1/2",  "0,1/2"  ])),
    "c2mm" => (["Γ", "Y", "S"],       KVec.(["0,0", "1,0",      "1/2,1/2"])),
    "p2mm" => (["Γ", "X", "S", "Y"],  KVec.(["0,0", "1/2,0",    "1/2,1/2",   "0,1/2"])),
    "p4"   => (["Γ", "X", "M"],       KVec.(["0,0", "0,1/2",    "1/2,1/2"])),
    "p4mm" => (["Γ", "X", "M" ],      KVec.(["0,0", "0,1/2",    "1/2,1/2"])),
    "p3"   => (["Γ", "KA", "K"],      KVec.(["0,0", "2/3,-1/3", "1/3,1/3"])),
    "p3m1" => (["Γ", "K", "M"],       KVec.(["0,0", "1/3,1/3",  "1/2,0"  ])),
    "p31m" => (["Γ", "KA", "K", "M"], KVec.(["0,0", "2/3,-1/3", "1/3,1/3",   "1/2,0"])),
    "p6"   => (["Γ", "K", "M"],       KVec.(["0,0", "1/3,1/3",  "1/2,0"  ])),
    "p6mm" => (["Γ", "K", "M"],       KVec.(["0,0", "1/3,1/3",  "1/2,0"  ]))
)

PLANEGROUP_IUC = Dict("p1"   => 1,  "p2" => 2,  "p1m1" => 3,  "c1m1" => 5,  "p2mm" => 6,
                      "c2mm" => 9,  "p4" => 10, "p4mm" => 11, "p3"   => 13, "p3m1" => 14,
                      "p31m" => 15, "p6" => 16, "p6mm" => 17)

function plane2littlegroup(planegroup::String)
    littlegroup.(Ref(spacegroup(PLANEGROUP_IUC[planegroup], Val(2))), 
                 PLANE2KVEC[planegroup][2])
end
function _find_isomorphic_parent_pointgroup(G)
    D = dim(G)
    ctᴳ = MultTable(reduce_ops(operations(G), 'p')).table
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
LGIRSD_2D = Dict{String, Dict{String, Vector{LGIrrep{2}}}}()
for (iuclab, (klabs, kvs)) in PLANE2KVEC
    lgs  = plane2littlegroup(iuclab) # the associated little group at the k kpoints
    ops′ = reduce_ops.(lgs)
    lgs .= LittleGroup.(num.(lgs), kvec.(lgs), klabel.(lgs), ops′)
    pgs  = _find_isomorphic_parent_pointgroup.(lgs) # find the parent point group 

    pglabs    = label.(pgs) #Labels of the parent point groups
    pgirreps  = get_pgirreps.(pglabs, Val(2)) #Use the labels of the point groups to find irreps    

    LGIRSD_2D[iuclab] = Dict{String, Vector{LGIrrep{2}}}()
    for (klab, kvec, lg, pgirs) in zip(klabs, kvs, lgs, pgirreps)
        cdmls        = Ref(klab).*last.(label.(pgirs))
        matrices     = getfield.(pgirs, Ref(:matrices))
        translations = [zeros(2) for _ in lg]
        type         = -1 # sentinel value for undetermined reality type; fixed below
        pgirs_types  = getfield.(pgirs, Ref(:type))
        lgirs = LGIrrep{2}.(cdmls, Ref(lg), matrices, Ref(translations), pgirs_types, false)
        LGIRSD_2D[iuclab][klab] = lgirs
    end
end

# correct the reality type by calling out to herring explicitly
LGIRSD_2D′ = Dict{String, Dict{String, Vector{LGIrrep{2}}}}()
for (iuclab, LGIRSD_2D) in LGIRSD_2D
    sgnum = num(first(LGIRSD_2D["Γ"]))
    sg    = reduce_ops(spacegroup(sgnum, Val(2)))
    LGIRSD_2D′[iuclab] = Dict{String, Vector{LGIrrep{2}}}()

    for (klab, lgirs) in LGIRSD_2D
        lgirs′ = LGIrrep{2}[]
        for lgir in lgirs
            lgir′ = LGIrrep{2}(lgir.cdml, lgir.lg, lgir.matrices, lgir.translations,
                               herring(lgir, sg), false)
            push!(lgirs′, lgir′)        
        end
        LGIRSD_2D′[iuclab][klab]= lgirs
    end
end
LGIRS_2D′ = [LGIRSD_2D for (key, LGIRSD_2D) in LGIRSD_2D′]