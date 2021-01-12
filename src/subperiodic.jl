# -- ASSOCIATIONS BETWEEN LAYER, ROD, & FRIEZE GROUPS VS. SPACE, PLANE, & LINE GROUPS ---
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

# TODO: ROD2SPACE_NUM
# TODO: LINE2ROD_NUM

# -- NOTATION FOR LAYER, ROD, AND FRIEZE GROUPS ---
# notation from https://www.cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen
const LAYERGROUP_IUCs = ( # 80 layer groups
    #=1=#  "p1",      #=2=#  "p-1",     #=3=#  "p112",   #=4=#  "p11m",   #=5=#  "p11a",
    #=6=#  "p112/m",  #=7=#  "p112/a",  #=8=#  "p211",   #=9=#  "p2₁11",  #=10=# "c211",
    #=11=# "pm11",    #=12=# "pb11",    #=13=# "cm11",   #=14=# "p2/m11", #=15=# "p2₁/m11",
    #=16=# "p2/b11",  #=17=# "p2₁/b11", #=18=# "c2/m11", #=19=# "p222",   #=20=# "p2₁22",
    #=21=# "p2₁2₁2",  #=22=# "c222",    #=23=# "pmm2",   #=24=# "pma2",   #=25=# "pba2",
    #=26=# "cmm2",    #=27=# "pm2m",    #=28=# "pm2₁b",  #=29=# "pb2₁m",  #=30=# "pb2b",
    #=31=# "pm2a",    #=32=# "pm2₁n",   #=33=# "pb2₁a",  #=34=# "pb2n",   #=35=# "cm2m",
    #=36=# "cm2e",    #=37=# "pmmm",    #=38=# "pmaa",   #=39=# "pban",   #=40=# "pmam",
    #=41=# "pmma",    #=42=# "pman",    #=43=# "pbaa",   #=44=# "pbam",   #=45=# "pbma",
    #=46=# "pmmn",    #=47=# "cmmm",    #=48=# "cmme",   #=49=# "p4",     #=50=# "p-4",
    #=51=# "p4/m",    #=52=# "p4/n",    #=53=# "p422",   #=54=# "p42₁2",  #=55=# "p4mm", 
    #=56=# "p4bm",    #=57=# "p-42m",   #=58=# "p-42₁m", #=59=# "p-4m2",  #=60=# "p-4b2", 
    #=61=# "p4/mmm",  #=62=# "p4/nbm",  #=63=# "p4/mbm", #=64=# "p4/nmm", #=65=# "p3",
    #=66=# "p-3",     #=67=# "p312",    #=68=# "p321",   #=69=# "p3m1",   #=70=# "p31m",
    #=71=# "p-31m",   #=72=# "p-3m1",   #=73=# "p6",     #=74=# "p-6",    #=75=# "p6/m", 
    #=76=# "p622",    #=77=# "p6mm",    #=78=# "p-6m2",  #=79=# "p-62m",  #=80=# "p6/mmm"
    )

const RODGROUP_IUCs    = ( # 75 rod groups 
                           # (for multiple setting choices, we always pick setting 1)
    #=1=#  "𝑝1",      #=2=#  "𝑝-1",     #=3=#  "𝑝211",   #=4=#  "𝑝m11",   #=5=#  "𝑝c11",
    #=6=#  "𝑝2/m11",  #=7=#  "𝑝2/c1",   #=8=#  "𝑝112",   #=9=#  "𝑝112₁",  #=10=# "𝑝11m",
    #=11=# "𝑝112/m",  #=12=# "𝑝112₁/m", #=13=# "𝑝222",   #=14=# "𝑝222₁",  #=15=# "𝑝mm2",
    #=16=# "𝑝cc2",    #=17=# "𝑝mc2₁",   #=18=# "𝑝2mm",   #=19=# "𝑝2cm",   #=20=# "𝑝mmm",
    #=21=# "𝑝ccm",    #=22=# "𝑝mcm",    #=23=# "𝑝4",     #=24=# "𝑝4₁",    #=25=# "𝑝4₂",
    #=26=# "𝑝4₃",     #=27=# "𝑝-4",     #=28=# "𝑝4/m",   #=29=# "𝑝4₂/m",  #=30=# "𝑝422",
    #=31=# "𝑝4₁22",   #=32=# "𝑝4₂22",   #=33=# "𝑝4₃22",  #=34=# "𝑝4mm",   #=35=# "𝑝4₂cm",
    #=36=# "𝑝4cc",    #=37=# "𝑝-42",    #=38=# "𝑝-42c",  #=39=# "𝑝4/mmm", #=40=# "𝑝4/mcc",
    #=41=# "𝑝4₂/mmc", #=42=# "𝑝3",      #=43=# "𝑝3₁",    #=44=# "𝑝3₂",    #=45=# "𝑝-3",
    #=46=# "𝑝312",    #=47=# "𝑝3₁12",   #=48=# "𝑝3₂12",  #=49=# "𝑝3m1",   #=50=# "𝑝3c1",
    #=51=# "𝑝-31m",   #=52=# "𝑝-31c",   #=53=# "𝑝6",     #=54=# "𝑝6₁",    #=55=# "𝑝6₂",
    #=56=# "𝑝6₃",     #=57=# "𝑝6₄",     #=58=# "𝑝6₅",    #=59=# "𝑝-6",    #=60=# "𝑝6/m",
    #=61=# "𝑝6₃/m",   #=62=# "𝑝622",    #=63=# "𝑝6₁22",  #=64=# "𝑝6₂22",  #=65=# "𝑝6₃22",
    #=66=# "𝑝6₄22",   #=67=# "𝑝6₅22",   #=68=# "𝑝6mm",   #=69=# "𝑝6cc",   #=70=# "𝑝6₃mc",
    #=71=# "𝑝-6m2",   #=72=# "𝑝-6c2",   #=73=# "𝑝6/mmm", #=74=# "𝑝6/mcc", #=75=# "𝑝6/mmc",
    )

const FRIEZEGROUP_IUCs = ( # 7 frieze groups
    #=1=#  "p1",      #=2=#  "p2",      #=3=#  "p1m1",   #=4=#  "p11m",   #=5=#  "p11g",
    #=6=#  "p2mm",    #=7=#  "p2mg"
    )

# Band topology check for layer groups:
#   w/ time-reversal symmetry
#       "Z₂"    = [2, 3, 7, 49, 50, 52, 66, 73]
#       "Z₂×Z₂" = [6, 51, 75]
#   w/o time-reversal symmetry: 
#       "Z₂"    = [2, 3, 7]
#       "Z₃"    = [65]
#       "Z₄"    = [49, 50, 52]
#       "Z₆"    = [66, 73]
#       "Z₂×Z₂" = [6]
#       "Z₃×Z₃" = [74]
#       "Z₄×Z₄" = [51]
#       "Z₆×Z₆" = [75]
# cf. Tables S18 and S20 of https://doi.org/10.1038/s41467-017-00133-2
