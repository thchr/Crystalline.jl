# By indexing into the following arrays, one obtains the "parent" group number, associated
# with the index' group number. As an example, `PLANE2SPACE_NUM[16] = 168`, meaning that the
# 2D plane group #16 has a parent in the 3D space group #168.

# Manual comparison to Bilbao's listings, and consulting Litvin's book's Table 30 (which
# doesn't fully give the correct 1-to-1 matches, because the conventions changed later on)
const PLANE2LAYER_NUMS = (
    1  #= p1   ⇒ p1   =#, 3  #= p2   ⇒ p112 =#, 11 #= p1m1 ⇒ pm11 =#,
    12 #= p1g1 ⇒ pb11 =#, 13 #= c1m1 ⇒ cm11 =#, 23 #= p2mm ⇒ pmm2 =#,
    24 #= p2mg ⇒ pma2 =#, 25 #= p2gg ⇒ pba2 =#, 26 #= c2mm ⇒ cmm2 =#,
    49 #= p4   ⇒ p4   =#, 55 #= p4mm ⇒ p4mm =#, 56 #= p4gm ⇒ p4bm =#,
    65 #= p3   ⇒ p3   =#, 69 #= p3m1 ⇒ p3m1 =#, 70 #= p31m ⇒ p3m1 =#,
    73 #= p6   ⇒ p6   =#, 77 #= p6mm ⇒ p6mm =#
    )

# Data from Table 1 of the SI of Watanabe, Po, and Vishwanath's 2017 Nature Commun.
const LAYER2SPACE_NUMS = (
    1, 2, 3, 6, 7, 10, 13, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 12, 16, 17,
    18, 21, 25, 28, 32, 35, 25, 26, 26, 27, 28, 31, 29, 30, 38, 39, 47,
    49, 50, 51, 51, 53, 54, 55, 57, 59, 65, 67, 75, 81, 83, 85, 89, 90,
    99, 100, 111, 113, 115, 117, 123, 125, 127, 129, 143, 147, 149, 150,
    156, 157, 162, 164, 168, 174, 175, 177, 183, 187, 189, 191
    )

# this is just `LAYER2SPACE_NUMS[[PLANE2LAYER_NUMS...]]`
const PLANE2SPACE_NUM = (
    1   #= p1   ⇒ P1   =#, 3   #= p2   ⇒ P2   =#, 6   #= p1m1 ⇒ Pm   =#,
    7   #= p1g1 ⇒ Pc   =#, 8   #= c1m1 ⇒ Cm   =#, 25  #= p2mm ⇒ Pmm2 =#,
    28  #= p2mg ⇒ Pma2 =#, 32  #= p2gg ⇒ Pba2 =#, 35  #= c2mm ⇒ Cmm2 =#,
    75  #= p4   ⇒ P4   =#, 99  #= p4mm ⇒ P4mm =#, 100 #= p4gm ⇒ P4bm =#,
    143 #= p3   ⇒ P3   =#, 156 #= p3m1 ⇒ P3m1 =#, 157 #= p31m ⇒ P31m =#,
    168 #= p6   ⇒ P6   =#, 183 #= p6mm ⇒ P6mm =#,
    )

const LINE2FRIEZE_NUM = (1 #= p1 ⇒ p1 =#, 3 #= p1m ⇒ p1m1 =#)

# notation from https://www.cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen
const LAYERGROUP_IUCs = ( # 80 layer groups
    "1",       "-1",     "112",    "11m",   "11a",    "112/m",  "112/a",   "211",
    "p2₁11",   "c211",   "pm11",   "pb11",  "cm11",   "p2/m11", "p2₁/m11", "p2/b11",
    "p2₁/b11", "c2/m11", "p222",   "p2₁22", "p2₁2₁2", "c222",   "pmm2",    "pma2",
    "pba2",    "cmm2",   "pm2m",   "pm2₁b", "pb2₁m",  "pb2b",   "pm2a",    "pm2₁n",
    "pb2₁a",   "pb2n",   "cm2m",   "cm2e",  "pmmm",   "pmaa",   "pban",    "pmam",
    "pmma",    "pman",   "pbaa",   "pbam",  "pbma",   "pmmn",   "cmmm",    "cmme",
    "p4",      "p-4",    "p4/m",   "p4/n",  "p422",   "p42₁2",  "p4mm",    "p4bm",
    "p-42m",   "p-42₁m", "p-4m2",  "p-4b2", "p4/mmm", "p4/nbm", "p4/mbm",  "p4/nmm",
    "p3",      "p-3",    "p312",   "p321",  "p3m1",   "p31m",   "p-31m",   "p-3m1",
    "p6",      "p-6",    "p6/m",   "p622",  "p6mm",   "p-6m2",  "p-62m",   "p6/mmm"
)

#const RODGROUP_IUCs    = () # 75 rod groups #TODO
const FRIEZEGROUP_IUCs = ("p1", "p2", "p1m1", "p11m", "p11g", "p2mm", "p2mg") # 7 frieze groups