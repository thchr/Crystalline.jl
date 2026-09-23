# ---------------------------------------------------------------------------------------- #
# NOTATION FOR LAYER, ROD, AND FRIEZE GROUPS
# notation from https://www.cryst.ehu.es/cgi-bin/subperiodic/programs/nph-sub_gen with a
# few manual corrections from Kopsky & Litvin's 'Nomenclature, Symbols and Classification
# of the Subperiodic Groups' report.

const LAYERGROUP_IUCs = ( # 80 layer groups
    #=1=#  "𝑝1",      #=2=#  "𝑝-1",     #=3=#  "𝑝112",   #=4=#  "𝑝11m",   #=5=#  "𝑝11a",
    #=6=#  "𝑝112/m",  #=7=#  "𝑝112/a",  #=8=#  "𝑝211",   #=9=#  "𝑝2₁11",  #=10=# "𝑐211",
    #=11=# "𝑝m11",    #=12=# "𝑝b11",    #=13=# "𝑐m11",   #=14=# "𝑝2/m11", #=15=# "𝑝2₁/m11",
    #=16=# "𝑝2/b11",  #=17=# "𝑝2₁/b11", #=18=# "𝑐2/m11", #=19=# "𝑝222",   #=20=# "𝑝2₁22",
    #=21=# "𝑝2₁2₁2",  #=22=# "𝑐222",    #=23=# "𝑝mm2",   #=24=# "𝑝ma2",   #=25=# "𝑝ba2",
    #=26=# "𝑐mm2",    #=27=# "𝑝m2m",    #=28=# "𝑝m2₁b",  #=29=# "𝑝b2₁m",  #=30=# "𝑝b2b",
    #=31=# "𝑝m2a",    #=32=# "𝑝m2₁n",   #=33=# "𝑝b2₁a",  #=34=# "𝑝b2n",   #=35=# "𝑐m2m",
    #=36=# "𝑐m2e",    #=37=# "𝑝mmm",    #=38=# "𝑝maa",   #=39=# "𝑝ban",   #=40=# "𝑝mam",
    #=41=# "𝑝mma",    #=42=# "𝑝man",    #=43=# "𝑝baa",   #=44=# "𝑝bam",   #=45=# "𝑝bma",
    #=46=# "𝑝mmn",    #=47=# "𝑐mmm",    #=48=# "𝑐mme",   #=49=# "𝑝4",     #=50=# "𝑝-4",
    #=51=# "𝑝4/m",    #=52=# "𝑝4/n",    #=53=# "𝑝422",   #=54=# "𝑝42₁2",  #=55=# "𝑝4mm", 
    #=56=# "𝑝4bm",    #=57=# "𝑝-42m",   #=58=# "𝑝-42₁m", #=59=# "𝑝-4m2",  #=60=# "𝑝-4b2", 
    #=61=# "𝑝4/mmm",  #=62=# "𝑝4/nbm",  #=63=# "𝑝4/mbm", #=64=# "𝑝4/nmm", #=65=# "𝑝3",
    #=66=# "𝑝-3",     #=67=# "𝑝312",    #=68=# "𝑝321",   #=69=# "𝑝3m1",   #=70=# "𝑝31m",
    #=71=# "𝑝-31m",   #=72=# "𝑝-3m1",   #=73=# "𝑝6",     #=74=# "𝑝-6",    #=75=# "𝑝6/m", 
    #=76=# "𝑝622",    #=77=# "𝑝6mm",    #=78=# "𝑝-6m2",  #=79=# "𝑝-62m",  #=80=# "𝑝6/mmm")

const RODGROUP_IUCs    = ( # 75 rod groups 
                           # (for multiple setting choices, we always pick setting 1)
    #=1=#  "𝓅1",      #=2=#  "𝓅-1",     #=3=#  "𝓅211",   #=4=#  "𝓅m11",   #=5=#  "𝓅c11",
    #=6=#  "𝓅2/m11",  #=7=#  "𝓅2/c1",   #=8=#  "𝓅112",   #=9=#  "𝓅112₁",  #=10=# "𝓅11m",
    #=11=# "𝓅112/m",  #=12=# "𝓅112₁/m", #=13=# "𝓅222",   #=14=# "𝓅222₁",  #=15=# "𝓅mm2",
    #=16=# "𝓅cc2",    #=17=# "𝓅mc2₁",   #=18=# "𝓅2mm",   #=19=# "𝓅2cm",   #=20=# "𝓅mmm",
    #=21=# "𝓅ccm",    #=22=# "𝓅mcm",    #=23=# "𝓅4",     #=24=# "𝓅4₁",    #=25=# "𝓅4₂",
    #=26=# "𝓅4₃",     #=27=# "𝓅-4",     #=28=# "𝓅4/m",   #=29=# "𝓅4₂/m",  #=30=# "𝓅422",
    #=31=# "𝓅4₁22",   #=32=# "𝓅4₂22",   #=33=# "𝓅4₃22",  #=34=# "𝓅4mm",   #=35=# "𝓅4₂cm",
    #=36=# "𝓅4cc",    #=37=# "𝓅-42m",   #=38=# "𝓅-42c",  #=39=# "𝓅4/mmm", #=40=# "𝓅4/mcc",
    #=41=# "𝓅4₂/mmc", #=42=# "𝓅3",      #=43=# "𝓅3₁",    #=44=# "𝓅3₂",    #=45=# "𝓅-3",
    #=46=# "𝓅312",    #=47=# "𝓅3₁12",   #=48=# "𝓅3₂12",  #=49=# "𝓅3m1",   #=50=# "𝓅3c1",
    #=51=# "𝓅-31m",   #=52=# "𝓅-31c",   #=53=# "𝓅6",     #=54=# "𝓅6₁",    #=55=# "𝓅6₂",
    #=56=# "𝓅6₃",     #=57=# "𝓅6₄",     #=58=# "𝓅6₅",    #=59=# "𝓅-6",    #=60=# "𝓅6/m",
    #=61=# "𝓅6₃/m",   #=62=# "𝓅622",    #=63=# "𝓅6₁22",  #=64=# "𝓅6₂22",  #=65=# "𝓅6₃22",
    #=66=# "𝓅6₄22",   #=67=# "𝓅6₅22",   #=68=# "𝓅6mm",   #=69=# "𝓅6cc",   #=70=# "𝓅6₃mc",
    #=71=# "𝓅-6m2",   #=72=# "𝓅-6c2",   #=73=# "𝓅6/mmm", #=74=# "𝓅6/mcc", #=75=# "𝓅6₃/mmc")

const FRIEZEGROUP_IUCs = ( # 7 frieze groups
    #=1=#  "𝓅1",      #=2=#  "𝓅2",      #=3=#  "𝓅1m1",   #=4=#  "𝓅11m",   #=5=#  "𝓅11g",
    #=6=#  "𝓅2mm",    #=7=#  "𝓅2mg")

# ---------------------------------------------------------------------------------------- #
# ASSOCIATIONS BETWEEN LAYER, ROD, & FRIEZE GROUPS VS. SPACE, PLANE, & LINE GROUPS
# By indexing into the following arrays, one obtains the "parent" group number, associated
# with the index' group number. As an example, `PLANE2SPACE_NUM[16] = 168`, meaning that the
# 2D plane group ⋕16 has a parent in the 3D space group ⋕168.

# Manual comparison to Bilbao's listings, and consulting Litvin's book's Table 30 (which
# doesn't fully give the correct 1-to-1 matches, because the conventions changed later on)
const PLANE2LAYER_NUMS = (
    1  #= p1   ⇒ 𝑝1   =#, 3  #= p2   ⇒ 𝑝112 =#, 11 #= p1m1 ⇒ 𝑝m11 =#,
    12 #= p1g1 ⇒ 𝑝b11 =#, 13 #= c1m1 ⇒ 𝑐m11 =#, 23 #= p2mm ⇒ 𝑝mm2 =#,
    24 #= p2mg ⇒ 𝑝ma2 =#, 25 #= p2gg ⇒ 𝑝ba2 =#, 26 #= c2mm ⇒ 𝑐mm2 =#,
    49 #= p4   ⇒ 𝑝4   =#, 55 #= p4mm ⇒ 𝑝4mm =#, 56 #= p4gm ⇒ 𝑝4bm =#,
    65 #= p3   ⇒ 𝑝3   =#, 69 #= p3m1 ⇒ 𝑝3m1 =#, 70 #= p31m ⇒ 𝑝3m1 =#,
    73 #= p6   ⇒ 𝑝6   =#, 77 #= p6mm ⇒ 𝑝6mm =#)

# Data from Table 1 of the SI of Watanabe, Po, and Vishwanath's 2017 Nature Commun.
const LAYER2SPACE_NUMS = (
    1   #= 𝑝1     ⇒ P1 =#,     2   #= 𝑝-1     ⇒ P-1 =#,    3   #= 𝑝112    ⇒ P2 =#, 
    6   #= 𝑝11m   ⇒ Pm =#,     7   #= 𝑝11a    ⇒ Pc =#,     10  #= 𝑝112/m  ⇒ P2/m =#,
    13  #= 𝑝112/a ⇒ P2/c =#,   3   #= 𝑝211    ⇒ P2 =#,     4   #= 𝑝2₁11   ⇒ P2₁ =#,
    5   #= 𝑐211   ⇒ C2 =#,     6   #= 𝑝m11    ⇒ Pm =#,     7   #= 𝑝b11    ⇒ Pc =#, 
    8   #= 𝑐m11   ⇒ Cm =#,     10  #= 𝑝2/m11  ⇒ P2/m =#,   11  #= 𝑝2₁/m11 ⇒ P2₁/m =#, 
    13  #= 𝑝2/b11 ⇒ P2/c =#,   14  #= 𝑝2₁/b11 ⇒ P2₁/c =#,  12  #= 𝑐2/m11  ⇒ C2/m =#,
    16  #= 𝑝222   ⇒ P222 =#,   17  #= 𝑝2₁22   ⇒ P222₁ =#,  18  #= 𝑝2₁2₁2  ⇒ P2₁2₁2 =#,
    21  #= 𝑐222   ⇒ C222 =#,   25  #= 𝑝mm2    ⇒ Pmm2 =#,   28  #= 𝑝ma2    ⇒ Pma2 =#, 
    32  #= 𝑝ba2   ⇒ Pba2 =#,   35  #= 𝑐mm2    ⇒ Cmm2 =#,   25  #= 𝑝m2m    ⇒ Pmm2 =#, 
    26  #= 𝑝m2₁b  ⇒ Pmc2₁ =#,  26  #= 𝑝b2₁m   ⇒ Pmc2₁ =#,  27  #= 𝑝b2b    ⇒ Pcc2 =#,
    28  #= 𝑝m2a   ⇒ Pma2 =#,   31  #= 𝑝m2₁n   ⇒ Pmn2₁ =#,  29  #= 𝑝b2₁a   ⇒ Pca2₁ =#, 
    30  #= 𝑝b2n   ⇒ Pnc2 =#,   38  #= 𝑐m2m    ⇒ Amm2 =#,   39  #= 𝑐m2e    ⇒ Aem2 =#, 
    47  #= 𝑝mmm   ⇒ Pmmm =#,   49  #= 𝑝maa    ⇒ Pccm =#,   50  #= 𝑝ban    ⇒ Pban =#, 
    51  #= 𝑝mam   ⇒ Pmma =#,   51  #= 𝑝mma    ⇒ Pmma =#,   53  #= 𝑝man    ⇒ Pmna =#, 
    54  #= 𝑝baa   ⇒ Pcca =#,   55  #= 𝑝bam    ⇒ Pbam =#,   57  #= 𝑝bma    ⇒ Pbcm =#, 
    59  #= 𝑝mmn   ⇒ Pmmn =#,   65  #= 𝑐mmm    ⇒ Cmmm =#,   67  #= 𝑐mme    ⇒ Cmme =#, 
    75  #= 𝑝4     ⇒ P4 =#,     81  #= 𝑝-4     ⇒ P-4 =#,    83  #= 𝑝4/m    ⇒ P4/m =#, 
    85  #= 𝑝4/n   ⇒ P4/n =#,   89  #= 𝑝422    ⇒ P422 =#,   90  #= 𝑝42₁2   ⇒ P42₁2 =#, 
    99  #= 𝑝4mm   ⇒ P4mm =#,   100 #= 𝑝4bm    ⇒ P4bm =#,   111 #= 𝑝-42m   ⇒ P-42m =#, 
    113 #= 𝑝-42₁m ⇒ P-42₁m =#, 115 #= 𝑝-4m2   ⇒ P-4m2 =#,  117 #= 𝑝-4b2   ⇒ P-4b2 =#,
    123 #= 𝑝4/mmm ⇒ P4/mmm =#, 125 #= 𝑝4/nbm  ⇒ P4/nbm =#, 127 #= 𝑝4/mbm  ⇒ P4/mbm =#,
    129 #= 𝑝4/nmm ⇒ P4/nmm =#, 143 #= 𝑝3      ⇒ P3 =#,     147 #= 𝑝-3     ⇒ P-3 =#,
    149 #= 𝑝312   ⇒ P312 =#,   150 #= 𝑝321    ⇒ P321 =#,   156 #= 𝑝3m1    ⇒ P3m1 =#, 
    157 #= 𝑝31m   ⇒ P31m =#,   162 #= 𝑝-31m   ⇒ P-31m =#,  164 #= 𝑝-3m1   ⇒ P-3m1 =#, 
    168 #= 𝑝6     ⇒ P6 =#,     174 #= 𝑝-6     ⇒ P-6 =#,    175 #= 𝑝6/m    ⇒ P6/m =#, 
    177 #= 𝑝622   ⇒ P622 =#,   183 #= 𝑝6mm    ⇒ P6mm =#,   187 #= 𝑝-6m2   ⇒ P-6m2 =#,
    189 #= 𝑝-62m  ⇒ P-62m =#,  191 #= 𝑝6/mmm  ⇒ P6/mmm =#)

# this is just `LAYER2SPACE_NUMS[[PLANE2LAYER_NUMS...]]`
const PLANE2SPACE_NUMS = (
    1   #= p1   ⇒ P1   =#, 3   #= p2   ⇒ P2   =#, 6   #= p1m1 ⇒ Pm   =#, 
    7   #= p1g1 ⇒ Pc   =#, 8   #= c1m1 ⇒ Cm   =#, 25  #= p2mm ⇒ Pmm2 =#,
    28  #= p2mg ⇒ Pma2 =#, 32  #= p2gg ⇒ Pba2 =#, 35  #= c2mm ⇒ Cmm2 =#,
    75  #= p4   ⇒ P4   =#, 99  #= p4mm ⇒ P4mm =#, 100 #= p4gm ⇒ P4bm =#,
    143 #= p3   ⇒ P3   =#, 156 #= p3m1 ⇒ P3m1 =#, 157 #= p31m ⇒ P31m =#,
    168 #= p6   ⇒ P6   =#, 183 #= p6mm ⇒ P6mm =#)

const FRIEZE2SPACE_NUMS = (
    1  #= 𝓅1   ⇒ P1   =#, 3 #= 𝓅2   ⇒ P2 =#, 6  #= 𝓅1m1 ⇒ Pm   =#,
    6  #= 𝓅11m ⇒ Pm   =#, 7 #= 𝓅11g ⇒ Pc =#, 25 #= 𝓅2mm ⇒ Pmm2 =#,
    28 #= 𝓅2mg ⇒ Pma2 =#)

const FRIEZE2LAYER_NUMS = (
    1  #= 𝓅1   ⇒ 𝑝1   =#, 3 #= 𝓅2   ⇒ 𝑝112 =#, 11 #= 𝓅1m1 ⇒ 𝑝m11 =#,
    4  #= 𝓅11m ⇒ 𝑝11m =#, 5 #= 𝓅11g ⇒ 𝑝11a =#, 23 #= 𝓅2mm ⇒ 𝑝mm2 =#,
    24 #= 𝓅2mg ⇒ 𝑝ma2 =#)

const FRIEZE2ROD_NUMS = (
     1  #= 𝓅1   ⇒ 𝓅1   =#, 3 #= 𝓅2   ⇒ 𝓅211 =#, 10 #= 𝓅1m1 ⇒ 𝓅11m =#, 
     4  #= 𝓅11m ⇒ 𝓅m11 =#, 5 #= 𝓅11g ⇒ 𝓅c11 =#, 18 #= 𝓅2mm ⇒ 𝓅2mm =#,
     19 #= 𝓅2mg ⇒ 𝓅2cm =#)

const FRIEZE2PLANE_NUMS  = (
    1 #= 𝓅1   ⇒ p1   =#, 2 #= 𝓅2   ⇒ p2   =#, 3 #= 𝓅1m1 ⇒ p1m1 =#,
    3 #= 𝓅11m ⇒ p1m1 =#, 4 #= 𝓅11g ⇒ p1g1 =#, 6 #= 𝓅2mm ⇒ p2mm =#,
    7 #= 𝓅2mg ⇒ p2mg =#)

const ROD2SPACE_NUMS = (
    1   #= 𝓅1     ⇒ P1     =#, 2   #= 𝓅-1     ⇒ P-1     =#, 3   #= 𝓅211    ⇒ P2      =#,
    6   #= 𝓅m11   ⇒ Pm     =#, 7   #= 𝓅c11    ⇒ Pc      =#, 10  #= 𝓅2/m11  ⇒ P2/m    =#,
    13  #= 𝓅2/c1  ⇒ P2/c   =#, 3   #= 𝓅112    ⇒ P2      =#, 4   #= 𝓅112₁   ⇒ P2₁     =#,
    6   #= 𝓅11m   ⇒ Pm     =#, 10  #= 𝓅112/m  ⇒ P2/m    =#, 11  #= 𝓅112₁/m ⇒ P2₁/m   =#,
    16  #= 𝓅222   ⇒ P222   =#, 17  #= 𝓅222₁   ⇒ P222₁   =#, 25  #= 𝓅mm2    ⇒ Pmm2    =#,
    27  #= 𝓅cc2   ⇒ Pcc2   =#, 26  #= 𝓅mc2₁   ⇒ Pmc2₁   =#, 25  #= 𝓅2mm    ⇒ Pmm2    =#,
    28  #= 𝓅2cm   ⇒ Pma2   =#, 47  #= 𝓅mmm    ⇒ Pmmm    =#, 49  #= 𝓅ccm    ⇒ Pccm    =#,
    63  #= 𝓅mcm   ⇒ Cmcm   =#, 75  #= 𝓅4      ⇒ P4      =#, 76  #= 𝓅4₁     ⇒ P4₁     =#,
    77  #= 𝓅4₂    ⇒ P4₂    =#, 78  #= 𝓅4₃     ⇒ P4₃     =#, 81  #= 𝓅-4     ⇒ P-4     =#,
    83  #= 𝓅4/m   ⇒  P4/m  =#, 84  #= 𝓅4₂/m   ⇒ P4₂/m   =#, 89  #= 𝓅422    ⇒ P422    =#,
    91  #= 𝓅4₁22  ⇒ P4₁22  =#, 93  #= 𝓅4₂22   ⇒ P4₂22   =#, 95  #= 𝓅4₃22   ⇒ P4₃22   =#,
    99  #= 𝓅4mm   ⇒  P4mm  =#, 101 #= 𝓅4₂cm   ⇒ P4₂cm   =#, 103 #= 𝓅4cc    ⇒ P4cc    =#,
    111 #= 𝓅-42m  ⇒ P-42m  =#, 112 #= 𝓅-42c   ⇒ P-42c   =#, 123 #= 𝓅4/mmm  ⇒ P4/mmm  =#,
    124 #= 𝓅4/mcc ⇒ P4/mcc =#, 131 #= 𝓅4₂/mmc ⇒ P4₂/mmc =#, 143 #= 𝓅3      ⇒ P3      =#,
    144 #= 𝓅3₁    ⇒  P3₁   =#, 145 #= 𝓅3₂     ⇒ P3₂     =#, 147 #= 𝓅-3     ⇒ P-3     =#,
    149 #= 𝓅312   ⇒ P312   =#, 151 #= 𝓅3₁12   ⇒ P3₁12   =#, 153 #= 𝓅3₂12   ⇒ P3₂12   =#,
    156 #= 𝓅3m1   ⇒ P3m1   =#, 158 #= 𝓅3c1    ⇒ P3c1    =#, 162 #= 𝓅-31m   ⇒ P-31m   =#,
    163 #= 𝓅-31c  ⇒ P-31c  =#, 168 #= 𝓅6      ⇒ P6      =#, 169 #= 𝓅6₁     ⇒ P6₁     =#,
    171 #= 𝓅6₂    ⇒ P6₂    =#, 173 #= 𝓅6₃     ⇒ P6₃     =#, 172 #= 𝓅6₄     ⇒ P6₄     =#,
    170 #= 𝓅6₅    ⇒ P6₅    =#, 174 #= 𝓅-6     ⇒ P-6     =#, 175 #= 𝓅6/m    ⇒ P6/m    =#,
    176 #= 𝓅6₃/m  ⇒ P6₃/m  =#, 177 #= 𝓅622    ⇒ P622    =#, 178 #= 𝓅6₁22   ⇒ P6₁22   =#,
    180 #= 𝓅6₂22  ⇒ P6₂22  =#, 182 #= 𝓅6₃22   ⇒ P6₃22   =#, 181 #= 𝓅6₄22   ⇒ P6₄22   =#,
    179 #= 𝓅6₅22  ⇒ P6₅22  =#, 183 #= 𝓅6mm    ⇒ P6mm    =#, 184 #= 𝓅6cc    ⇒ P6cc    =#,
    186 #= 𝓅6₃mc  ⇒ P6₃mc  =#, 187 #= 𝓅-6m2   ⇒ P-6m2   =#, 188 #= 𝓅-6c2   ⇒ P-6c2   =#,
    191 #= 𝓅6/mmm ⇒ P6/mmm =#, 192 #= 𝓅6/mcc  ⇒ P6/mcc  =#, 194 #= 𝓅6/mmc  ⇒ P6₃/mmc =#)

const LINE2ROD_NUMS = (1 #= p1 ⇒ 𝓅1 =#, 10 #= p1m ⇒ 𝓅11m =#)

const LINE2FRIEZE_NUMS = (1 #= p1 ⇒ 𝓅1 =#, 3 #= p1m ⇒ 𝓅1m1 =#)

# ---------------------------------------------------------------------------------------- #
# GROUP ELEMENTS

"""
$(TYPEDEF)$(TYPEDFIELDS)

A subperiodic group of embedding dimension `D` and periodicity dimension `P`. 

Fields: 
- `operations`: the `SymOperation`s of the finite factor group ``G/T``, where ``G`` is the
subperiodic group and ``T`` is the translation group of the associated lattice.
- `num`: the canonical number of the group, following the International Tables for
Crystallography, Volume E.
"""
struct SubperiodicGroup{D,P} <: AbstractSpaceGroup{D, SymOperation{D}}
    num :: Int
    operations :: Vector{SymOperation{D}}
end
const LayerGroup  = SubperiodicGroup{3,2}
const RodGroup    = SubperiodicGroup{3,1}
const FriezeGroup = SubperiodicGroup{2,1}

function _throw_subperiodic_domain(D::Integer, P::Integer)
    throw(DomainError((D, P), "invalid dimension and periodicity for subperiodic group"))
end

@noinline function _throw_subperiodic_num(num::Integer, D::Integer, P::Integer)
    maxnum, sub = (D==3 && P==2) ? (80, "layer") :
                  (D==3 && P==1) ? (75, "rod") :
                  (D==2 && P==1) ? (7, "frieze") :
                  _throw_subperiodic_domain(D, P)

    throw(DomainError(num,
        "group number must be between 1 and $maxnum for $sub groups (D=$D, P=$P)"))
end

@inline function _subperiodic_kind(D, P)
    # TODO: Move to a testing-utils module (only used in /test/)
    if D == 3 && P == 2
        return "layer"
    elseif D == 3 && P == 1
        return "rod"
    elseif D == 2 && P == 1
        return "frieze"
    else
        _throw_subperiodic_domain(D, P)
    end
end

function _check_valid_subperiodic_num_and_dim(num::Integer, D::Integer, P::Integer)
    if D == 3 && P == 2     # layer groups
        num > 80 && _throw_subperiodic_num(num, D, P)
    elseif D == 3 && P == 1 # rod groups
        num > 75  && _throw_subperiodic_num(num, D, P)
    elseif D == 2 && P == 1 # frieze groups
        num > 7   && _throw_subperiodic_num(num, D, P)
    else
        _throw_subperiodic_domain(D,P)
    end
    num < 1 && throw(DomainError(num, "group number must be a positive integer"))
    return nothing
end

label(g::SubperiodicGroup{D,P}) where {D,P} = _subperiodic_label(num(g), D, P)
iuc(g::SubperiodicGroup) = label(g)
@inline function _subperiodic_label(num::Integer, D::Integer, P::Integer)
    @boundscheck _check_valid_subperiodic_num_and_dim(num, D, P)
    if D == 3 && P == 2
        return LAYERGROUP_IUCs[num]
    elseif D == 3 && P == 1
        return RODGROUP_IUCs[num]
    elseif D == 2 && P == 1
        return FRIEZEGROUP_IUCs[num]
    else
        _throw_subperiodic_domain(D, P) # unreachable under boundschecking
    end
end

centering(g::SubperiodicGroup{D,P}) where {D,P} = centering(num(g), D, P)
function centering(num::Integer, D::Integer, P::Integer)
    if D == P
        return centering(num, D)
    else
        lab = _subperiodic_label(num, D, P) # (also checks input validity)
        return first(lab)
    end
    error("unreachable")
end

# ---------------------------------------------------------------------------------------- #

function reduce_ops(subg::SubperiodicGroup{<:Any,P}, conv_or_prim::Bool=true,
                    modw::Bool=true) where P
    return reduce_ops(operations(subg), centering(subg), conv_or_prim, modw, Val(P))
end

function primitivize(subg::SubperiodicGroup{<:Any,P}, modw::Bool=true) where P
    return typeof(subg)(num(subg), reduce_ops(subg, false, modw))
end

# TODO: conventionalize(::SubperiodicGroup)

# ---------------------------------------------------------------------------------------- #

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

# ---------------------------------------------------------------------------------------- #