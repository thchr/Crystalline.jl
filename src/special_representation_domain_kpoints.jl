# The basic goal of this endeavor is only to find the irreps that 
# are missing in ISOTROPY but feature in the bandreps from Bilbao. 
# There are 145 such irreps, across 20 space groups, they are:                                (... I manually pieced this table together)
# ┌───────┬─────────────────────────────────────────────────────────────────────────────┬───────────────┬────────────────────────────────┐
# │       │  Missing                                                                    │  Missing      │  [NS≡nonsymmorph; S≡symmorph]  │
# │  SGs  │  LGIrrep labels                                                             │  KVec labels  │  Match method                  │
# │───────┼─────────────────────────────────────────────────────────────────────────────┼───────────────┼────────────────────────────────│
# │  23   │  WA₁, WA₂, WA₃, WA₄                                                         │  WA           │  S:  MonoOrthTetraCubic ┐      │
# │  24   │  WA₁                                                                        │  WA           │  NS: Inherits from      └ 23   │
# │  82   │  PA₁, PA₂, PA₃, PA₄                                                         │  PA           │  S:  MonoOrthTetraCubic        │
# │  121  │  PA₁, PA₂, PA₃, PA₄, PA₅                                                    │  PA           │  S:  MonoOrthTetraCubic ┐      │
# │  122  │  PA₁, PA₂                                                                   │  PA           │  NS: Inherits from      └ 121  │
# │  143  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex                    │
# │  144  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┐         │
# │  145  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┘         │
# │  150  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex                    │
# │  152  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┐         │
# │  154  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Orphan (type b) ┘         │
# │  157  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  S:  TriHex ───────┐           │
# │  159  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, PA₁, PA₂, PA₃                                │  HA, KA, PA   │  NS: Inherits from └ 157       │
# │  174  │  HA₁, HA₂, HA₃, HA₄, HA₅, HA₆, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃  │  HA, KA, PA   │  S:  TriHex                    │
# │  189  │  HA₁, HA₂, HA₃, HA₄, HA₅, HA₆, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃  │  HA, KA, PA   │  S:  TriHex ───────┐           │
# │  190  │  HA₁, HA₂, HA₃, KA₁, KA₂, KA₃, KA₄, KA₅, KA₆, PA₁, PA₂, PA₃                 │  HA, KA, PA   │  NS: Inherits from └ 189       │
# │  197  │  PA₁, PA₂, PA₃, PA₄                                                         │  PA           │  S:  MonoOrthTetraCubic ┐      │
# │  199  │  PA₁, PA₂, PA₃                                                              │  PA           │  NS: Inherits from      └ 197  │
# │  217  │  PA₁, PA₂, PA₃, PA₄, PA₅                                                    │  PA           │  S:  MonoOrthTetraCubic ┐      │
# │  220  │  PA₁, PA₂, PA₃                                                              │  PA           │  NS: Inherits from      └ 217  │
# └───────┴─────────────────────────────────────────────────────────────────────────────┴───────────────┴────────────────────────────────┘
# Of course, in principle, we might consider trying to create ALL the missing points from Φ-Ω,
# but for what we need, this is sufficient.

# --- STRUCTS ---

struct KVecMapping
    kᴬlab::String
    kᴮlab::String
    op::SymOperation
end
function show(io::IO, ::MIME"text/plain", kvmap::KVecMapping)
    println(io, kvmap.kᴬlab, " => ", kvmap.kᴮlab, " via R = ", xyzt(kvmap.op))
end


# --- HARDCODED CONSTANTS ---
# Cached results of `Tuple(tuple(_find_holosymmetric_sgnums(D)...) for D = 1:3)`
HOLOSYMMETRIC_SGNUMS = (
    (2,),                                                            # 1D
    (2, 6, 7, 8, 9, 11, 12, 17),                                     # 2D
    (2, 10, 11, 12, 13, 14, 15, 47, 48, 49, 50, 51, 52, 53, 54, 55,  # 3D
    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 
    133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 166, 167, 191,
    192, 193, 194, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230)
)

# Holosymmetric point group labels (get corresponding pointgroups from with get_pgsops(...))
# ordered in ascending point group order (number of operators).
# For 3D, we used Table 2 of J. Appl. Cryst. (2018). 51, 1481–1491 (https://doi.org/10.1107/S1600576718012724)
# For 2D and 1D, the cases can be inferred from case-by-case enumeration.
const HOLOSYMMETRIC_PGS_IUCLAB = (
    ("m",),                                                         # 1D
    ("2", "mm2", "4mm",  "6mm"),                                    # 2D
    ("-1", "2/m", "mmm", "-31m", "-3m1", "4/mmm", "6/mmm", "m-3m")  # 3D
)
# Table 1 of J. Appl. Cryst. (2018). 51, 1481–1491 (https://doi.org/10.1107/S1600576718012724)
const HOLOSYMMETRIC_PG_FOR_BRAVAISTYPE = ImmutableDict(
    "aP"                       =>     ["-1"], # (not ideal to be using Vector here, but meh...)
    (("mP", "mC")             .=> Ref(["2/m"]))...,
    (("oP", "oI", "oF", "oC") .=> Ref(["mmm"]))...,
    "hR"                       =>     ["-31m","-3m1"], # special case: two possible settings (-31m and -3m1)
    "hP"                       =>     ["6/mmm"],
    (("tP", "tI")             .=> Ref(["4/mmm"]))...,
    (("cP", "cI", "cF")       .=> Ref(["m-3m"]))...
)

# Mnemonized/cached data from calling 
# Tuple(tuple(getindex.(get_arith_crystalclass_partner.(1:MAX_SGNUM[D], D), 2)...) for D in 1:3)
const ISOGONAL_PARENT_GROUPS = (
    (1,2),                                                                                          # 1D
    (1,2,3,3,5,6,6,6,9,10,11,11,13,14,15,16,17),                                                    # 2D
    (1,2,3,3,5,6,6,8,8,10,10,12,10,10,12,16,16,16,16,21,21,22,23,23,25,25,25,25,25,25,25,25,25,25,  # 3D
     35,35,35,38,38,38,38,42,42,44,44,44,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,65,65,65,
     65,65,65,69,69,71,71,71,71,75,75,75,75,79,79,81,82,83,83,83,83,87,87,89,89,89,89,89,89,89,89,
     97,97,99,99,99,99,99,99,99,99,107,107,107,107,111,111,111,111,115,115,115,115,119,119,121,121,
     123,123,123,123,123,123,123,123,123,123,123,123,123,123,123,123,139,139,139,139,143,143,143,146,
     147,148,149,150,149,150,149,150,155,156,157,156,157,160,160,162,162,164,164,166,166,168,168,168,
     168,168,168,174,175,175,177,177,177,177,177,177,183,183,183,183,187,187,189,189,191,191,191,191,
     195,196,197,195,197,200,200,202,202,204,200,204,207,207,209,209,211,207,207,211,215,216,217,215,
     216,217,221,221,221,221,225,225,225,225,229,229)
)

# Orphan space group numbers in (3D only)
const ORPHAN_SGNUMS = (
    (198,),                                                     # type (a) in B&C p. 414
    (76, 78, 144, 145, 151, 152, 153, 154, 169, 170, 171, 172), # type (b)
    (91, 95, 92, 96, 178, 179, 180, 181, 212, 213),             # type (c); no new k-vecs though
    (205,)                                                      # type (d)
)

# The group to supergroups (G => G′) map from CDML Table 4.1; necessary to 
# construct the irrep mapping for orphans of type (a) and (b)
const ORPHAN_AB_SUPERPARENT_SGNUMS = ImmutableDict(
    76=>91, 78=>95, 144=>178, 145=>179, 151=>181, 152=>178,
    153=>180, 154=>179, 169=>178, 170=>179, 171=>180, 172=>181,
    198=>212)

# Dict of group (G) => supergroup (G₀) relations, along with their transformation operators,
# for tricky corner cases that cannot be treated by a naïve subgroup check which 
# doesn't account for the changes in crystal setting between space groups (centering,
# orientation etc.); this is rather tedious. We manually read off the appropriate
# supergroups - and picked one from the list of holosymmetric sgs - and then subsequently
# manually verified that that choice makes G a normal/invariant of G₀, i.e. that G◁G₀.
# This extraction was done using Bilbao's MINSUP program (www.cryst.ehu.es/cryst/minsup.html),
# so the setting is already consistent with ITA.
# Abbreviations below: min-sup-sg ≡ minimal (normal and holosymmetric) supergroup
const CORNERCASES_SUBSUPER_NORMAL_SGS = SGOps.ImmutableDict(
     # group sgnum => (supergroup sgnum, transformation rotation P, transformation translation p)
     17 => (51,  copy.(unpack(SymOperation("z,x,y")))...),            # min-sup-sg
     26 => (51,  copy.(unpack(SymOperation("z,x,y")))...),            # min-sup-sg
     28 => (51,  copy.(unpack(SymOperation("x,-z,y")))...),           # min-sup-sg
     29 => (54,  copy.(unpack(SymOperation("-z,y,x+1/4")))...),       # min-sup-sg
     30 => (52,  copy.(unpack(SymOperation("z,x+1/4,y+1/4")))...),    # min-sup-sg
     33 => (52,  copy.(unpack(SymOperation("x+1/4,z,-y+1/4")))...),   # min-sup-sg
     38 => (63,  copy.(unpack(SymOperation("y,z,x+1/4")))...),        # min-sup-sg
     39 => (64,  copy.(unpack(SymOperation("y+1/4,z,x+1/4")))...),    # min-sup-sg
     40 => (63,  copy.(unpack(SymOperation("-z,y,x")))...),           # min-sup-sg
     41 => (64,  copy.(unpack(SymOperation("-z,y,x")))...),           # min-sup-sg
     43 => (70,  copy.(unpack(SymOperation("z,x+3/8,y+3/8")))...),    # min-sup-sg
     46 => (72,  copy.(unpack(SymOperation("-z,y,x+1/4")))...),       # min-sup-sg
     80 => (141, copy.(unpack(SymOperation("x+1/2,y+3/4,z")))...),    # NOT min-sup-sg: instead, a cycle through 88=>141 which ensures normality/invariance 
     86 => (133, copy.(unpack(SymOperation("x,y+1/2,z")))...),        # min-sup-sg
     88 => (141, copy.(unpack(SymOperation("x,y+1/2,z")))...),        # min-sup-sg
     90 => (127, copy.(unpack(SymOperation("x,y+1/2,z")))...),        # min-sup-sg
     98 => (141, copy.(unpack(SymOperation("x,y+1/4,z+3/8")))...),    # min-sup-sg
    109 => (141, copy.(unpack(SymOperation("x,y+1/4,z")))...),        # min-sup-sg
    110 => (142, copy.(unpack(SymOperation("x,y+1/4,z")))...),        # min-sup-sg
    122 => (141, copy.(unpack(SymOperation("x,y+1/4,z+3/8")))...),    # min-sup-sg
    210 => (227, copy.(unpack(SymOperation("x+3/8,y+3/8,z+3/8")))...) # min-sup-sg
)

# Transformation matrices from CDML to ITA settings
include((@__DIR__)*"/../data/transformation_matrices_CDML2ITA.jl") # ⇒ defines TRANSFORMS_CDML2ITA::ImmutableDict 


# --- FUNCTIONS ---

function is_orphan_sg(sgnum::Integer, D::Integer=3)
    D ≠ 3 && throw_2d_not_yet_implemented(D)  # 2D not considered in CDML
    for orphantypeidx in eachindex(ORPHAN_SGNUMS)
        sgnum ∈ ORPHAN_SGNUMS[orphantypeidx] && return orphantypeidx
    end
    return 0 # ⇒ not an orphan
end

"""
    _find_holosymmetric_sgnums(D::Integer)

We compute the list of holosymmetric space group numbers by first finding the "maximum"
isogonal point group of each Bravais type (looping through all the space groups
in that Bravais type); then we subsequently compare the isogonal point groups of 
each space group to this maximal (Bravais-type-specific) point group; if they agree
the space group is holosymmetric.

See `is_holosymmetric` for description of holosymmetric space groups and of 
their connection to the representation and basic domains Φ and Ω.
"""
function _find_holosymmetric_sgnums(D::Integer)
    bravaistypes = bravaistype.(Base.OneTo(MAX_SGNUM[D]),D)
    uniquebravaistypes = unique(bravaistypes) # Bravais types (1, 5, & 14 in 1D, 2D, & 3D)

    # find maximum point groups for each bravais type
    maxpointgroups = Dict(ubt=>Vector{SymOperation}() for ubt in uniquebravaistypes)
    for (sgnum,bt) in enumerate(bravaistypes)
        pg = sort(pointgroup(get_sgops(sgnum,D)), by=xyzt)
        if length(pg) > length(maxpointgroups[bt])
            maxpointgroups[bt] = pg;
        end
    end

    # determine whether each space group is a holosymmetric space group
    # then accumulate the `sgnum`s of the holosymmetric space groups
    holosymmetric_sgnums = Vector{Int64}()
    for (sgnum,bt) in enumerate(bravaistypes)
        pg = sort(pointgroup(get_sgops(sgnum,D)), by=xyzt)
        if length(pg) == length(maxpointgroups[bt])
            push!(holosymmetric_sgnums, sgnum)
        end
    end

    return holosymmetric_sgnums
end

"""
    is_holosymmetric(sgnum::Integer, D::Integer) --> Bool

Return a Boolean answer for whether the representation domain Φ equals
the basic domain Ω, i.e. whether the space group is holosymmetric (see
CDML p. 31 and 56). Φ and Ω are defined such that the Brillouin zone BZ
can be generated from Φ through the point group-parts of the space group
operations g∈P (the "isogonal point group" is CDML) and from Ω through
the crystal system's point group operations g∈P₀ (the "holosymmetric
point group") i.e. BZ≡∑_(g∈Ḡ)gΦ and BZ≡∑_(r∈Gᶜʳʸˢ)rΦ. If Φ=Ω, we say 
that the space group is holosymmetric; otherwise, Φ is an integer multiple
of Ω and we say that the space group is non-holosymmetric.

In practice, rather than compute explicitly every time, we use a cache of
holosymmetric space group numbers obtained from `_find_holosymmetric_sgnums`
(from the `const` `HOLOSYMMERIC_SGNUMS`).
"""
is_holosymmetric(sgnum::Integer, D::Integer=3) = (sgnum ∈ HOLOSYMMETRIC_SGNUMS[D])
is_holosymmetric(sg::SpaceGroup) = is_holosymmetric(num(sg), dim(sg))


"""
    find_holosymmetric_supergroup(G::SpaceGroup)
    find_holosymmetric_supergroup(sgnum::Integer, D::Integer)
                --> PointGroup

Finds the minimal holosymmetric super point group `P` of a space group `G`
(alternatively specified by its number `sgnum` and dimension `D`), such
that the isogonal point group of `G`, denoted `F`, is a subgroup of `P`,
i.e. such that `F`≤`P` with `P` restricted to the lattice type of `G` (see
`HOLOSYMMETRIC_PG_FOR_BRAVAISTYPE`). For holosymmetric space groups `F=P`.
"""
function find_holosymmetric_superpointgroup(G::SpaceGroup)
    D = dim(G)
    F = pointgroup(G) # isogonal point group of G (::Vector{SymOperation})
    if D == 3
        # In 3D there are cases (162, 163, 164, 165) where the distinctions
        # between hP and hR need to be accounted for explicitly, so there we
        # have to explicitly ensure that we only compare with pointgroups 
        # that indeed match the lattice type (hR and hP have different holohedries)
        bt = bravaistype(num(G), D)
        for pglab in HOLOSYMMETRIC_PG_FOR_BRAVAISTYPE[bt]
            # this check is actually redunant for everything except the hR groups; we do it anyway
            P = get_pgops(pglab, D) # holosymmetric point group (::PointGroup)
            if issubgroup(operations(P), F)
                return P
            end
        end

    else
        # No tricky cornercasees in 2D or 1D: we can iterate through 
        # the holosymmetric point groups to find the minimal holosymmetric
        # supergroup, because we already sorted HOLOSYMMETRIC_PGS_IUCLAB 
        # by the point group order      
        for pglab in HOLOSYMMETRIC_PGS_IUCLAB[D]
            P = get_pgops(pglab, D) # holosymmetric point group (::PointGroup)
            if issubgroup(operations(P), F)
                return P
            end
        end
    end
    throw("Did not find a holosymmetric super point group: if the setting is conventional, this should never happen")
end
find_holosymmetric_superpointgroup(sgnum::Integer, D::Integer=3) = find_holosymmetric_superpointgroup(get_sgops(sgnum, D))


"""
    find_holosymmetric_parent(sgnum::Integer, D::Integer)
                --> Union{SpaceGroup, Nothing}

Find a holosymmetric space group G₀ such that the space group G with number `sgnum`
is an invariant subgroup of G₀, i.e. G◁G₀: a "holosymmetric parent spacegroup" of G, 
in the notation of B&C and CDML.
This identification is not necessarily unique; several distinct parents may exist.

The meaning of invariant subgroup (also called a "normal subgroup") is that G is
both a subgroup of G₀, i.e. G<G₀, and G is invariant under conjugation by any
element g₀ of G₀.
Thus, to find if G is an invariant subgroup of G₀ we first check if G is a subgroup
of a G₀; then check that for every g₀∈G₀ and g∈G we have g₀gg₀⁻¹∈G. If we find 
these conditions to be fulfilled, we return the space group G₀. 

The search is complicated by the fact that the appropriate G₀ may be in a different
setting than G is. At the moment a poor-man's fix is implemented that manually 
circumvents this case by explicit tabulation of cornercases; though this is only 
tested in 3D.

If we find no match, we return `nothing`: the space groups for which this is true
should agree with those in `ORPHAN_SGNUMS`.
"""
function find_holosymmetric_parent(sgnum::Integer, D::Integer=3)
    D ≠ 3 && throw_2d_not_yet_implemented()
    G = get_sgops(sgnum, D)

    if !is_holosymmetric(sgnum, D) # nontrivial case: find invariant subgroup
        cntr = centering(sgnum, D)

        if !haskey(CORNERCASES_SUBSUPER_NORMAL_SGS, sgnum)
            # Naïve attempt to find a holosymmetric supergroup; this fails in 21 cases:
            # see data/cornercases_normal_subsupergroup_pairs_with_transformations.jl.
            # Basically, this naïve approach assumes that the crystal setting between 
            # G and its supergroup G⁰ agrees: it does not in general.
            for sgnum₀ in HOLOSYMMETRIC_SGNUMS[D]
                G₀ = get_sgops(sgnum₀, D)
                cntr₀ = centering(sgnum₀, D)

                # === check whether G is a subgroup of G₀, G<G₀ ===
                # won't be a naïve subgroup if centerings are different (effectively, a cheap short-circuit check)
                cntr ≠ cntr₀ && continue
                # now check if G<G₀ operator-element-wise
                issubgroup(G₀, G)[1] || continue

                # check if G is an _invariant_ subgroup of G₀, i.e. if g₀gg₀⁻¹∈G
                # for every g₀∈G₀ and g∈G
                isnormal(G₀, G) && return G₀
            end

        else    # G must be a cornercase
            # The above naïve search does not work generally because different sgs are in
            # different settings; even when the Bravais types match, different choices can
            # be made for the coordinate axis (rotation, orientation, origo, handedness).
            # In general, that means that two sets of symmetry operations can be considered
            # to be a group-subgroup pair if there exists a single transformation operator 
            # such that their operators agree. 
            # Put more abstractly, in the above we check for a normal supergroup; but what 
            # we really want is to check for a normal supergroup TYPE (the notion of "space 
            # group type" differentiates from a concrete "space group", signifying a specific
            # choice of setting: usually, the distinction is ignored.).
            # We manually found off the appropriate normal and holosymmetric supergroup by
            # using Bilbao's MINSUP program, along with the appropriate transformation matrix;
            # see data/cornercases_normal_subsupergroup_pairs_with_transformations.jl and 
            # the constant Tuple `CORNERCASES_SUBSUPER_NORMAL_SGS`.
            # So far only checked for 3D.
            sgnum₀, P₀, p₀ = CORNERCASES_SUBSUPER_NORMAL_SGS[sgnum]

            G₀ = get_sgops(sgnum₀, 3)
            # find the supergroup G₀ in its transformed setting using P₀ and p₀
            opsG₀′ = transform.(operations(G₀), Ref(P₀), Ref(p₀))
            G₀′ = SpaceGroup(sgnum₀, opsG₀′, D)

            # verify that with G₀′ in this setting, G is a subgroup of G₀′, G is normal in G₀′
            # and that G₀′ is a holosymmetric space group
            @assert issubgroup(G₀′, G)
            @assert isnormal(G₀′, G)
            @assert is_holosymmetric(G₀′)

            return G₀′
        end

        # G doesn't have a holosymmetric parent space group ⇒ one of the "orphans" from CDML/B&C
        return nothing 
    else # trivial case; a holosymmetric sg is its own parent
        return G
    end
end

"""
    find_map_from_Ω_to_ΦnotΩ(G::SpaceGroup)
    find_map_from_Ω_to_ΦnotΩ(sgnum::Integer, D::Integer)

Find the point group operations `Rs` that map the basic domain Ω to the 
"missing" domains in Φ-Ω for the space group `G`; this is akin to a coset
decomposition of the isogonal point group `F` of `G` into the holosymmetric 
super point group `P` of `G`
"""
function find_map_from_Ω_to_ΦnotΩ(G::SpaceGroup)
    if is_holosymmetric(num(G), dim(G))
        return nothing
    else # G is not a holosymmetric sg, so Φ-Ω is finite
        P = operations(find_holosymmetric_superpointgroup(G)) # holosymmetric supergroup of G (::Vector{SymOperation})
        F = pointgroup(G)                                     # isogonal point group of G     (::Vector{SymOperation})

        index::Int64 = length(P)/length(F) # index of F in P = 2^"number of needed maps"
        if index == 1; println(num(G)); return nothing; end
        N::Int64 = log2(index)
        Rs = Vector{SymOperation}(undef, N)
        matchidx = 0
        # Find N coset representatives
        for i in 1:N
            matchidx = findnext(op->op∉F, P, matchidx+1) # find a coset representative for F in P
            Rs[i] = P[matchidx]
            if N > 1 && i < N
                # a bit of extra work to make sure we find an operator which
                # is not already in the coset P[matchidx]∘F by growing F with
                # the "newly added" elements that P[matchidx] produce
                F = append!(F, Ref(Rs[i]).∘F)
            end
        end
        return Rs
    end
end
find_map_from_Ω_to_ΦnotΩ(sgnum::Integer, D::Integer=3) = find_map_from_Ω_to_ΦnotΩ(get_sgops(sgnum, D))
# A sanity check for the above implementation, is to compare the number of distinct maps
# obtained by first defining 
#   Nmaps = [length(ops) for ops in SGOps.find_map_from_Ω_to_ΦnotΩ.(keys(SGOps.ΦNOTΩ_KVECS_AND_MAPS),3)]
# against the tabulated maps
#   maxfreq(x) = maximum(values(StatsBase.countmap(x)))
#   Nmapsᶜᵈᵐˡ = [maxfreq([x.kᴬlab for x in mapset]) for mapset in values(SGOps.ΦNOTΩ_KVECS_AND_MAPS)]
# Then, Nmaps .> Nmapsᶜᵈᵐˡ (not equal, necessarily, because some maps might be trivially related to
# a k-star or to a G-vector). 

function find_new_kvecs(G::SpaceGroup)
    Rs = find_map_from_Ω_to_ΦnotΩ(G)
    Rs === nothing && return nothing # holosymmetric case
    cntr = centering(num(G), dim(G))
    
    # load the KVecs in Ω from the ISOTROPY dataset
    lgs = get_littlegroups(num(G), dim(G))
    # We are only interested in mapping kvs from the basic domain Ω; but ISOTROPY already 
    # includes some of the ZA points that are in Φ-Ω, so we strip these already here (and
    # so find these points anew effectively). Also, there is no point in trying to map the
    # general point (also denoted Ω by us) to a new point, since it can be obtained by 
    # parameter variation; so we filter that out as well.
    filter!(lg-> length(klabel(lg)) == 1 && klabel(lg) != "Ω", lgs)
    kvs = kvec.(lgs)
    klabs = klabel.(lgs)
    Nk = length(lgs)    

    newkvs = [KVec[] for _ in Base.OneTo(Nk)]
    newklabs = [String[] for _ in Base.OneTo(Nk)]
    for (kidx, (kv, klab)) in enumerate(zip(kvs, klabs))
        for (Ridx, R) in enumerate(Rs)
            newkv = R∘kv
            # compare newkv to the stars of kv: if it is equivalent to any member of star{kv},
            # it is not a "new" k-vector. We do not have to compare to _all_ kv′∈Ω (i.e. to all
            # kvs), because they must have inequivalent stars to kv in the first place and also
            # cannot be mapped to them by a point group operation of the isogonal parent lattice
            kvstar = kstar(G, kv)
            if any(kvˢᵗᵃʳ->isapprox(newkv, kvˢᵗᵃʳ, cntr), kvstar)
                continue # jump out of loop if newkv is equivalent to any star{kv′}
            end
    
            # if we are not careful we may now add a "new" k-vector which is equivalent 
            # to a previously added new k-vector (i.e. is equivalent to a k-vector in its 
            # star): this can e.g. happen if R₁ and R₂ maps to the same KVec, which is a real
            # possibility. We check against the k-star of the new k-vectors just added.
            newkv_bool_vs_ΦnotΩ = true
            for newkv′ in newkvs[kidx]
                newkv′star = kstar(G, newkv′) # k-star of previously added new k-vector
                if any(kvˢᵗᵃʳ->isapprox(newkv, kvˢᵗᵃʳ, cntr), newkv′star)
                    newkv_bool_vs_ΦnotΩ = false
                    continue
                end
            end
            !newkv_bool_vs_ΦnotΩ && continue # jump out of loop if newkv is equivalent to any just added new KVec

            # newkv must now be a bonafide "new" KVec which is inequivalent to any original KVecs
            # kv′∈Ω, as well as any elements in their stars, and similarly inequivalent to any 
            # previously added new KVecs newkv′∈Φ-Ω      ⇒      add it to the list of new KVecs!
            push!(newkvs[kidx],   newkv)
            push!(newklabs[kidx], klab*('@'+Ridx)*'′') # new labels, e.g. if klab = "Γ", then newklab = "ΓA′", "ΓB′", "ΓC′", ...

            # TODO: It seems we would need to rule out additional k-vectors by some rule
            #       that I don't quite know yet. One option may be to consider k-vectors
            #       that can be obtained by parameter-variation equivalent; I think this 
            #       is what is called the "uni-arm" description in 
            #           https://www.cryst.ehu.es/cryst/help/definitions_kvec.html#uniarm
            #       Another concern could be that the symmetry of the reciprocal lattice
            #       is not same as that of the direct lattice; that is also discussed in 
            #       link. In general, we should be finding the same k-vectors as those 
            #       generated by https://www.cryst.ehu.es/cryst/get_kvec.html; right now
            #       we occassionally find too many (148 is a key example; 75 another).
            #       See also their companion paper to that utility.
            #       
            #       Easiest ways to inspect our results is as:
            #           sgnum = 75
            #           v=SGOps.find_new_kvecs(sgnum, 3)
            #           [v[3] string.(v[1]) first.(v[4]) string.(first.(v[2]))] # (only shows "A" generated KVecs though...)
            #
            #       We also have some debugging code for this issue in test/holosymmetric.jl,
            #       but at present it is not a priority to fix this, since we have the nominally
            #       correct additions from CDML.
        end
    end

    # only retain the cases where something "new" was added
    mappedkvs   = [kvs[kidx] for (kidx, newkvsₖᵥ) in enumerate(newkvs) if !isempty(newkvsₖᵥ)]
    mappedklabs = [klabs[kidx] for (kidx, newkvsₖᵥ) in enumerate(newkvs) if !isempty(newkvsₖᵥ)]
    filter!(!isempty, newkvs)
    filter!(!isempty, newklabs)
    # return arrays of {original KVec in Ω}, {new KVec(s) in Φ-Ω}, {original KVec label in Ω}, {new KVec label(s) in Φ-Ω},
    return mappedkvs, newkvs, mappedklabs, newklabs
end
find_new_kvecs(sgnum::Integer, D::Integer=3) = find_new_kvecs(get_sgops(sgnum, D))

"""
    get_arith_crystalclass_partner(sg::SpaceGroup)
    get_arith_crystalclass_partner(sgnum::Integer, D::Integer=3)

Find the symmorphic partner space group `sg′` of a space group `sg`, i.e.
find the space group `sg′` whose group embodies the arithmetic crystal class
of `sg`. In practice, this means that they have the same point group operations
and that they share the same Bravais lattice (i.e. centering) and Brillouin zone. 
For symmorphic space groups, this is simply the space group itself.
Another way of putting this is to find the isogonal parent group of `sg`.

The answer is returned as a pair of space group numbers: `num(sg)=>num(sg′)`.
"""
function get_arith_crystalclass_partner(sg::SpaceGroup)
    D = dim(sg)

    sgnum = num(sg)
    cntr = centering(sgnum, D)
    if sgnum ∈ SYMMORPH_SGNUMS[D]
        return sgnum=>sgnum
        
    else
        pgops_sorted = sort(pointgroup(sg), by=xyzt)
        N_pgops = length(pgops_sorted)
        N_arith_crys_class = length(SYMMORPH_SGNUMS[D])
        # exploit that the isogonal partner typically is close to the index
        # of the space group, but is (at least for 3D space groups) always smaller
        idx₀ = findfirst(x->sgnum<x, SYMMORPH_SGNUMS[D])
        if idx₀ === nothing # no point group of smaller index than the provided sg
            idx₀ = N_arith_crys_class
        end 
        for idx′ in Iterators.flatten((idx₀-1:-1:1, idx₀:N_arith_crys_class))
            sgnum′ = SYMMORPH_SGNUMS[D][idx′]
            # We have to take the `pointgroup(...)` operation here even though `sgops′` 
            # is symmorphic, because the `SymOperation`s are in a conventional basis 
            # and may contain contain trivial primitive translations if `cntr′ ≠ 'P'`
            sgops′ = sort(pointgroup(get_sgops(sgnum′, D)), by=xyzt)
            if (length(sgops′) == N_pgops &&
                all(sgops′ .== pgops_sorted) && # ... this assumes the coordinate setting
                                                # is the same between sg and sg′ (not
                                                # generally valid; but seems to work here,
                                                # probably because the setting is shared
                                                # by definition between these "partners"?)
                centering(sgnum′, D) == cntr) # the centering (actually, Bravais type, 
                                              # but doesn't seem to make a difference) 
                                              # must also agree!
                return sgnum=>sgnum′
            end
        end
    end

    # if we reached this point, no match was found ⇒ implies an error..!
    throw(ErrorException("No isogonal partner group was found: this should"*
                         "never happen for well-defined groups."))
end

function get_arith_crystalclass_partner(sgnum::Integer, D::Integer=3)
    get_arith_crystalclass_partner(get_sgops(sgnum, D))
end


"""
    _ΦnotΩ_kvecs_and_maps_imdict(;verbose::Bool=false)

Load the special k-points kᴮ∈Φ-Ω that live in the representation domain Φ, but not 
the basic domain Ω, as well parent space group operation R that link them to a related
k-point in the basic domain kᴬ∈Ω. Returns an `ImmutableDict` over `KVecMapping` `Tuple`s, 
each a containing {kᴬ, kᴮ, R} triad, with kᴬ and kᴮ represented by their labels.
The keys of the `ImmutableDict` give the associated space group number.
The data used here is obtained from Table 3.11 of CDML. Because they split symmetry 
operations up into two groups, {monoclinic, orthorhombic, tetragonal, and cubic} vs. 
{trigonal, hexagonal}, we do the same when we load it; but the end-result concatenates
both types into a single `ImmutableDict`.

For convenience we later write the results of this function to constants
`ΦNOTΩ_KVECS_AND_MAPS`.
"""
function _ΦnotΩ_kvecs_and_maps_imdict(;verbose::Bool=false)
    # prepare for loading csv files into ImmutableDict
    baseloadpath = (@__DIR__)*"/../data/CDML_RepresentationDomainSpecialKPoints_"
    kwargs = (header=["kᴮ_num", "kᴮ", "R", "kᴬ"], comment="#", delim=", ",
              ignoreemptylines=true, types=[Int64, String, Int64, String])

    d = Base.ImmutableDict{Int64, NTuple{N, KVecMapping} where N}()
    for loadtype in ["MonoOrthTetraCubic", "TriHex"]
        # load crystal-specific variables
        loadpath = baseloadpath*loadtype*".csv"
        if loadtype == "MonoOrthTetraCubic"
            # Monoclinic, cubic, tetragonal, orthorhombic
            num2ops = Dict( 2=>SymOperation("x,-y,-z"), 3=>SymOperation("-x,y,-z"), 
                            4=>SymOperation("-x,-y,z"), 16=>SymOperation("y,x,-z"), 
                            26=>SymOperation("-x,y,z"), 28=>SymOperation("x,y,-z"),
                            37=>SymOperation("y,x,z") )
        elseif loadtype == "TriHex"
            # Trigonal & hexagonal
            num2ops = Dict( 2=>SymOperation("x-y,x,z"),    6=>SymOperation("y,-x+y,z"), 
                            8=>SymOperation("x,x-y,-z"),   9=>SymOperation("y,x,-z"), 
                            10=>SymOperation("-x+y,y,-z"), 16=>SymOperation("x,y,-z"),
                            23=>SymOperation("x,x-y,z"),   24=>SymOperation("y,x,z") )
        end

        # read data from csv file
        data = CSV.read(loadpath; kwargs...)

        # convert csv file to ImmutableDict of KVecMapping NTuples
        cur_sgnum = 0
        tempvec = Vector{KVecMapping}()
        for (idx, row) in enumerate(eachrow(data))
            if iszero(row[:kᴮ_num])
                if cur_sgnum ≠ 0
                    d = Base.ImmutableDict(d, cur_sgnum=>tuple(tempvec...))
                end
                cur_sgnum = row[:R]
                tempvec = Vector{KVecMapping}()
            else
                R = num2ops[row[:R]]
                # CDML setting is not the same as ITA settings for all cases; for
                # the cur_sgnums where the two differ, we store the transformations
                # between the two settings in TRANSFORMS_CDML2ITA and apply it here 
                if haskey(TRANSFORMS_CDML2ITA, cur_sgnum) 
                    P, p = TRANSFORMS_CDML2ITA[cur_sgnum]
                    if verbose
                        Pstr = xyzt(SymOperation(hcat(P, p)))
                        println("SG $(cur_sgnum): transforming from CDML to ITA setting via P = $(Pstr):\n",
                                "    From R = $(xyzt(R))\t= $(seitz(R))")
                    end
                    R = transform(R, P, p)

                    verbose && println("     to R′ = $(xyzt(R))\t= $(seitz(R))")
                end
                push!(tempvec, KVecMapping(row[:kᴬ], row[:kᴮ], R))
            end
            if idx == size(data, 1)
                d = Base.ImmutableDict(d, cur_sgnum=>tuple(tempvec...))
            end
        end
    end
    return d
end



""" 
    ΦnotΩ_kvecs(sgnum::Integer, D::Integer)

For a given space group number `sgnum` and dimensionality `D`, find 
the parent holosymmetric group - or, if an orphan, its parent supergroup, 
and return the special representation domain Φ k-points kᴮ∈Φ-Ω that do
not reside in the basic domain Ω. Return a tuple of kᴮs, their associated
basic-domain kᴬs, and the transformation from one to the other R, such 
that kᴮ = Rkᴬ.

Returns `nothing` if there are no new special k-points in Φ-Ω.

Presently only works in 3D, since CDML does not tabulate the special 
k-points in Φ-Ω in 2D.
"""
function ΦnotΩ_kvecs(sgnum::Integer, D::Integer=3)
    D ≠ 3 && throw_2d_not_yet_implemented(D)

    # No new k-points for holosymmetric sgs where Φ = Ω
    is_holosymmetric(sgnum, D) && return nothing 

    # Check if sg is an "orphan" type (if 0 ⇒ not an orphan) in the sense discussed in CDML p. 70-71.
    orphantype = is_orphan_sg(sgnum, D)
    if orphantype == 3 # no new k-vectors after all (CDML p. 71, bullet (i))   [type (c)]
        return nothing

    elseif orphantype == 4 # sgnum = 205 needs special treatment               [type (d)]
        # TODO: Return ZA, AA, and BA irreps from CDML p. C491
        #       ZA is already included in ISOTROPY, and we don't really need AA or BA?
        @warn "not yet implemented"
        return missing
        
    elseif orphantype ∈ (1,2) # supergroup case                                [types (a,b)]
        # A supergroup G′ exists which contains some (but not all) necessary operations 
        # that otherwise exist in the holosymmetric group; see table 4.1 of CDML; by exploiting
        # their connection to an isomorphic partner, we can get all the necessary operations
        # (see B&C).
        parent_sgnum = ORPHAN_AB_SUPERPARENT_SGNUM[sgnum]
    end

    # If we reached here, orphantype is either 0 - in which case we can follow
    # the holosymmetric parent space group (G₀) approach - or orphantype is 
    # 1 (a) or 2 (b); in that case, we should use the supergroup G′ to identify
    # the irrep mapping vectors {R|v}; but we can still use the point group operation
    # R inferred from ISOGONAL_PARENT_GROUPS and ΦNOTΩ_KVECS_AND_MAPS (if featuring)

    # We need the isogonal parent group to find the correct key in ΦNOTΩ_KVECS_AND_MAPS
    # (since it only stores symmorphic space groups)
    isogonal_parent_sgnum = ISOGONAL_PARENT_GROUPS[D][sgnum]

    # There are some space groups (i.e. some isogonal super groups) that just
    # do not get truly new k-points in Φ-Ω since they either already feature in 
    # star{𝐤} or are equivalent (either via a primitive reciprocal G-vector or 
    # by trivial parameter variation) to an existing KVec in Ω. In that case,
    # there is no entry for isogonal_parent_sgnum in ΦNOTΩ_KVECS_AND_MAPS and
    # there is nothing to add. This is the case for sgnums 1, 16, 22, 89, 148,
    # 177, 207, & 211). Terminate early in that case.
    !haskey(ΦNOTΩ_KVECS_AND_MAPS, isogonal_parent_sgnum) && return nothing

    # Now we know there is a matching set of mappings in ΦNOTΩ_KVECS_AND_MAPS 
    # that contains ALL the new k-vecs and R mappings from CDML Table 3.11; get it
    kvmaps_pg = ΦNOTΩ_KVECS_AND_MAPS[isogonal_parent_sgnum] # point group operations R

    # now we need to find the appropriate parent group: if not an orphan, 
    # this is the parent holosymmetric space group G₀; if an orphan of type
    # (a) or (b) it is the supergroup G′ (already stored then)
    if iszero(orphantype)
        parent_sg = find_holosymmetric_parent(sgnum, D)
        parent_sg === nothing && throw("Unexpected error")
    else
        # TODO: load the orphan's parent_sg from parent_sgnum (maybe consolidate approaches?)
    end

    # The mapping kvmaps_isogonal refers to an operation {R|0} which may or may not 
    # exist in the parent space group `parent_sg`; but an operation {R|v} _will_ exist:
    # we now find that operation (this is the operation needed when mapping the irreps)
    # TODO: Continue implementation from here?
    kvmaps = Vector{KVecMapping}()
    for kvmap_pg in kvmaps_pg
        pgop = kvmap_pg.op                    # R in CDML notation
        matchidx = findfirst(op′->rotation(pgop)==rotation(op′), operations(parent_sg))
        if matchidx === nothing; throw("Should never happen: R = $(xyzt(pgop)) and parent_sgnum = $(num(parent_sg))"); end
        sgop = operations(parent_sg)[matchidx] # {R|v} in CDML notation

        push!(kvmaps, KVecMapping(kvmap_pg.kᴬlab, kvmap_pg.kᴮlab, sgop))
    end

    return kvmaps
end



# -- "GENERATED" CONSTANTS ---

# Mnemonized data from calling `_ΦnotΩ_kvecs_and_maps_imdict()` (in 3D only) as an ImmutableDict
const ΦNOTΩ_KVECS_AND_MAPS = _ΦnotΩ_kvecs_and_maps_imdict()