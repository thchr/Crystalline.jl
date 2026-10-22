"""
    LinePlaneLittleGroups

Little groups of the plane (2D) and line (1D) groups: their tabulated k-manifolds, with
labels, and their operations. Written to disk by `write_littlegroups`.
"""
module LinePlaneLittleGroups

using Crystalline, JLD2
using Crystalline: reduce_ops

export write_littlegroups, PLANE2KVEC, LINE2KVEC

# The listings below include all unique, labelled k-points in the representation domain
# whose little group is not trivial, and every time-reversal invariant momentum (TRIM), even
# if its little group is trivial (as in p1, c1m1, & p3): with time-reversal, such a TRIM is
# still special, e.g., hosting Kramers pairs for spinful particles. In addition, Γ=(0,0)
# and Ω=(u,v) are included for every group. The notation follows Bilbao's LKVEC tool [^1]
# entirely. Note that this implies some occasionally odd-looking choices, e.g. (see (*)
# marks below):
#   - *Plane groups 3 & 4*: the X and Y labels are defined counter-intuitively, such that
#           X = (0,1/2) and Y = (1/2,0). This also differs from the definitions in the
#           other primitive rectangular (op) Bravais lattices (i.e. 6, 7, & 8), which
#           instead have X = (1/2,0) etc. The difference also exists for k-line labels.
#   - *Plane group 5 vs. 9*: although both are centred rectangular (oc) Bravais lattices,
#           their k-vector labels do not agree, but differ by coordinate permutations.
# The origin for these differences lies in different settings for the plane group
# operations. Because we use Bilbao's conventions for plane group operations, we also stick
# with their conventions for k-vector labels here, even though they are somewhat unsightly
# occasionally.
# Litvin & Wike's 'Character tables and compatibility relations of the eighty layer groups
# and seventeen plane groups' (1991), give labels that do not have these deficiencies (and
# otherwise agree with LKVEC) in table 24 and figure 7. However, this implies a different
# setting for the symmetry operations, so we cannot really adopt it.
# [^1]: [Bilbao's LKVEC](https://cryst.ehu.es/subperiodic/get_layer_kvec.html)
const PLANE2KVEC = [
    #=1 =# (["Y", "B", "A"],  KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2"])),
    #=2 =# (["Y", "B", "A"],  KVec{2}.(["0,1/2",   "1/2,0",   "1/2,-1/2"])),
    #=3 =# (["X", "Y", "S",                         #= (*) =#
            "Σ", "ΣA", "C",
            "CA"],           KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "0,u",   "0,-u",  "1/2,u",    "1/2,-u"])),
    #=4 =# (["X", "Y", "S",   #= non-symmorphic =#  #= (*) =#
            "Σ", "ΣA", "C",
            "CA"],           KVec{2}.(["0,1/2",   "1/2,0",   "1/2,1/2", "0,u",   "0,-u",  "1/2,u",    "1/2,-u"])),
    #=5 =# (["Y", "S", "Σ",                         #= (*) =#
            "ΣA", "C", "CA"],KVec{2}.(["0,1",     "1/2,1/2", "0,2u",    "0,-2u", "1,2u",  "1,-2u"])),
    #=6 =# (["X", "Y", "S",
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["1/2,0",   "0,1/2",   "1/2,1/2", "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    #=7 =# (["X", "Y", "S",   #= non-symmorphic =#
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["1/2,0",   "0,1/2",   "1/2,1/2", "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    #=8 =# (["X", "Y", "S",   #= non-symmorphic =#
            "Σ", "C", "Δ",
            "D"],            KVec{2}.(["1/2,0",   "0,1/2",   "1/2,1/2", "u,0",   "u,1/2", "0,u",      "1/2,u"])),
    #=9 =# (["Y", "S", "Σ",
            "Δ", "F", "C"],  KVec{2}.(["1,0",     "1/2,1/2", "2u,0",    "0,2u",  "1,2u",  "2u,1"])),
    #=10=# (["X", "M"],       KVec{2}.(["0,1/2",   "1/2,1/2"])),
    #=11=# (["X", "M", "Σ",
            "Δ", "Y"],       KVec{2}.(["0,1/2",   "1/2,1/2", "u,u",     "0,u",   "u,1/2"])),
    #=12=# (["X", "M", "Σ",   #= non-symmorphic =#
            "Δ", "Y"],       KVec{2}.(["0,1/2",   "1/2,1/2", "u,u",     "0,u",   "u,1/2"])),
    #=13=# (["K", "KA", "M"], KVec{2}.(["1/3,1/3", "2/3,-1/3", "1/2,0"])),
    #=14=# (["K", "M", "Σ",
            "ΣA"],           KVec{2}.(["1/3,1/3", "1/2,0",   "u,0",   "0,u"])),
    #=15=# (["K", "KA", "M",
            "Λ", "ΛA", "T",
            "TA"],           KVec{2}.(["1/3,1/3", "2/3,-1/3", "1/2,0",  "u,u",   "2u,-u", "1/2-u,2u", "1/2+u,-2u"])),
    #=16=# (["K", "M"],       KVec{2}.(["1/3,1/3", "1/2,0"])),
    #=17=# (["K", "M", "Σ",
            "Λ", "T"],       KVec{2}.(["1/3,1/3", "1/2,0",   "u,0",   "u,u", "1/2-u,2u"]))
]
# Line groups: the only k-point with a nontrivial little group, besides Γ, is X=(1/2)
const LINE2KVEC = [(["X"], KVec{1}.(["1/2"])) for n in 1:MAX_SGNUM[1]]

# push Γ & Ω k-points, which are in all line + plane groups, to front & end of above vectors
for (D, klabs2kvs) in ((1, LINE2KVEC), (2, PLANE2KVEC))
    for (klabs, kvs) in klabs2kvs
        pushfirst!(klabs, "Γ"); pushfirst!(kvs, KVec{D}(D == 1 ? "0" : "0,0"))
        push!(klabs, "Ω");      push!(kvs, KVec{D}(D == 1 ? "α" : "α,β"))
    end
end

"""
    write_littlegroups(::Val{D}) --> String

Write the little groups of all D-dimensional space groups, at the k-manifolds of
`PLANE2KVEC` (`D = 2`) or `LINE2KVEC` (`D = 1`), to
`data/irreps/lgs/{D}d/littlegroups_data.jld2`, in the layout of the 3D data (see
`write_littlegroup_irreps.jl`). The operations are those of the space group, without
centring copies.
"""
function write_littlegroups(::Val{D}) where D
    D in (1, 2) || throw(ArgumentError("D must be 1 or 2"))
    outpath = joinpath(
        dirname(@__DIR__), "data", "irreps", "lgs", "$(D)d", "littlegroups_data.jld2"
    )
    # Crystalline keeps the data files open for reading, which blocks overwriting them
    jldfile = Crystalline.LGS_JLDFILES[D]
    isassigned(jldfile) && close(jldfile[])

    JLD2.jldopen(outpath, "w") do f
        klabs2kvs = D == 1 ? LINE2KVEC : PLANE2KVEC
        for (sgnum, (klabs, kvs)) in enumerate(klabs2kvs)
            sg = spacegroup(sgnum, Val(D))
            sgops = reduce_ops(sg)
            # indices into `sgops` of the operations of each little group
            opsidxs = [first(littlegroup(sgops, kv, centering(sg))) for kv in kvs]
            f["$(sgnum)/sgops"]       = xyzt.(sgops)
            f["$(sgnum)/klab_list"]   = klabs
            f["$(sgnum)/kstr_list"]   = [filter(!isspace, chop(string(kv); head=1, tail=1)) for kv in kvs]
            f["$(sgnum)/opsidx_list"] = [Int16.(idxs) for idxs in opsidxs]
        end
    end
    jldfile[] = JLD2.jldopen(outpath, "r") # reopen, for subsequent use of `littlegroups`
    println("wrote $(outpath)  ($(round(filesize(outpath)/2^10; digits=1)) KiB)")
    return outpath
end

end # module LinePlaneLittleGroups
