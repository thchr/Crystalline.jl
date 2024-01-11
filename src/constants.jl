# --- SCALARS ---
const DEFAULT_ATOL = 1e-12        # absolute tolerance for approximate equality
const NULL_ATOL    = 1e-11        # absolute tolerance for nullspace

@doc raw"""
    MAX_SGNUM :: Tuple{Int,Int,Int}

Return the number of distinct space group types across dimensions 1, 2, and 3 (indexable by
dimensionality).
"""
const MAX_SGNUM    = (2, 17, 230)

@doc raw"""
    MAX_SUBGNUM :: ImmutableDict

An immutable dictionary with values `v::Int` and keys `k::Tuple{Int,Int}`, where `v` is
the number of distinct subperiodic group types for a given key `k = (D,P)` describing a
subperiodic group of dimensionality `D` and periodicity `P`:

- layer groups: `(D,P) = (3,2)`
- rod groups: `(D,P) = (3,1)`
- frieze groups: `(D,P) = (2,1)`
"""
const MAX_SUBGNUM = ImmutableDict((3,2)=>80, (3,1)=>75, (2,1)=>7)

@doc raw"""
    MAX_MSGNUM :: Tuple{Int,Int,Int}

Analogous to `MAX_SUBGNUM`, but for the number of magnetic space groups.
"""
const MAX_MSGNUM  = (7, 80, 1651) # Figure 1.2.1 of Litvin's book


@doc raw"""
    MAX_MSUBGNUM :: Tuple{Int,Int,Int}

Analogous to `MAX_SUBGNUM`, but for the number of magnetic subperiodic groups.
"""
const MAX_MSUBGNUM = ImmutableDict((3,2)=>528, (3,1)=>394, (2,1)=>31) # Litvin, Fig. 1.2.1

# --- VECTORS ---
# arbitrary test vector for e.g. evaluating KVecs lines/planes/volumes; 3D by default
const TEST_αβγ =  [0.123123123, 0.456456456, 0.789789789]
const TEST_αβγs = (TEST_αβγ[1:1], TEST_αβγ[1:2], TEST_αβγ)

# --- LISTS ---
# The following space group numbers are symmorphic; in 3D they each embody one of 
# the 73 arithmetic crystal classes. Basically, this is the combination of the 14  
# Bravais lattices with the 32 point groups/geometric crystal classes. Obtained from e.g.
#   > tuple((tuple([i for i in 1:MAX_SGNUM[D] if issymmorph(i, D)]...) for D = 1:3)...)
const SYMMORPH_SGNUMS = (# 1D (1st index; all elements are symmorph)
                         (1, 2), 
                         # 2D (2nd index)
                         ((1, 2, 3, 5, 6, 9, 10, 11, 13, 14, 15, 16, 17)), 
                         # 3D (3rd index)
                         (1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47,65,
                          69,71,75,79,81,82,83,87,89,97,99,107,111,115,119,121,
                          123,139,143,146,147,148,149,150,155,156,157,160,162,
                          164,166,168,174,175,177,183,187,189,191,195,196,197,
                          200,202,204,207,209,211,215,216,217,221,225,229)
                        )


@doc raw"""
    ENANTIOMORPHIC_PAIRS :: NTuple{11, Pair{Int,Int}}

Return the space group numbers of the 11 enantiomorphic space group pairs in 3D.

The space group types associated with each such pair `(sgnum, sgnum')` are related by a
mirror transformation: i.e. there exists a transformation 
``\mathbb{P} = \{\mathbf{P}|\mathbf{p}\}`` between the two groups ``G = \{g\}`` and
``G' = \{g'\}`` such that ``G' = \mathbb{P}^{-1}G\mathbb{P}`` where ``\mathbf{P}`` is
improper (i.e. ``\mathrm{det}\mathbf{P} < 0``).

We define distinct space group *types* as those that cannot be related by a proper
transformation (i.e. with ``\mathrm{det}\mathbf{P} > 0``). With that view, there are 230
space group types. If the condition is relaxed to allow improper rotations, there are 
``230-11 = 219`` distinct *affine* space group types. See e.g. ITA5 Section 8.2.2.

The enantiomorphic space group types are also chiral space group types in 3D. There are no
enantiomorphic pairs in lower dimensions; in 3D all enantiomorphic pairs involve screw
symmetries, whose direction is inverted between pairs (i.e. have opposite handedness).
"""
const ENANTIOMORPHIC_PAIRS = (76 => 78,   91 => 95,   92 => 96,   144 => 145, 152 => 154,
                              151 => 153, 169 => 170, 171 => 172, 178 => 179, 180 => 181,
                              213 => 212)
# this is just the cached result of
#   pairs_label = [ # cf. ITA5 p. 727
#       "P4₁" => "P4₃", "P4₁22" => "P4₃22", "P4₁2₁2" => "P4₃2₁2", "P3₁" => "P3₂", 
#       "P3₁21" => "P3₂21", "P3₁12" => "P3₂12", "P6₁" => "P6₅", "P6₂" => "P6₄", 
#       "P6₁22" => "P6₅22", "P6₂22" => "P6₄22", "P4₁32" => "P4₃32"
#       ]
#   iucs = iuc.(1:230)
#   pairs = map(((x,y),) -> findfirst(==(x), iucs) => findfirst(==(y), iucs), ps)