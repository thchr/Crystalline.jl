# generated from original `xyzt`-based data via:
#   io = stdout
#   println(io, "SG_GENS_CODES_Ds = (")
#   for D in 1:3
#       sgnums = 1:MAX_SGNUM[D]
#       gensv = generators_from_xyzt.(sgnums, SpaceGroup{D})
#       ROTATIONS = D == 3 ? Crystalline.ROTATIONS_3D : 
#                   D == 2 ? Crystalline.ROTATIONS_2D :
#                   D == 1 ? Crystalline.ROTATIONS_1D : error()
#       TRANSLATIONS = D == 3 ? Crystalline.TRANSLATIONS_3D : 
#                      D == 2 ? Crystalline.TRANSLATIONS_2D :
#                      D == 1 ? Crystalline.TRANSLATIONS_1D : error()
#       println(io, "# ", SpaceGroup{D})
#       println(io, "Vector{Tuple{Int8,Int8}}[")
#       for (sgnum, gens) in zip(sgnums, gensv)
#           Nops = length(gens)
#           pg_codes = Vector{Tuple{Int8, Int8}}(undef, Nops)
#           idx = 0
#           for (n,op) in enumerate(gens)
#               r, t = rotation(op), translation(op)
#               r_idx = Int8(@something findfirst(==(r), ROTATIONS) error("could not find r = $r"))
#               t_idx = Int8(@something findfirst(==(t), TRANSLATIONS) error("could not find t = $t"))
#               pg_codes[idx+=1] = (r_idx, t_idx)
#           end  
#           print(io, "#=", sgnum, "=# [")
#           for (idx, pg_code) in enumerate(pg_codes)
#               print(io, "(")
#               join(io, pg_code, ",") 
#               print(io, ")")
#               idx ≠ length(pg_codes) && print(io, ",")
#           end
#           println(io, "],")
#       end
#       println(io, "],")
#   end
#   println(io, ")")

const SG_GENS_CODES_Vs = (
# SpaceGroup{1}
Vector{Tuple{Int8,Int8}}[
#=1=# [(1,1)],
#=2=# [(2,1)],
],
# SpaceGroup{2}
Vector{Tuple{Int8,Int8}}[
#=1=# [(1,1)],
#=2=# [(2,1)],
#=3=# [(3,1)],
#=4=# [(3,2)],
#=5=# [(3,1),(1,3)],
#=6=# [(2,1),(3,1)],
#=7=# [(2,1),(3,4)],
#=8=# [(2,1),(3,3)],
#=9=# [(2,1),(3,1),(1,3)],
#=10=# [(2,1),(5,1)],
#=11=# [(2,1),(5,1),(3,1)],
#=12=# [(2,1),(5,1),(3,3)],
#=13=# [(9,1)],
#=14=# [(9,1),(8,1)],
#=15=# [(9,1),(7,1)],
#=16=# [(9,1),(2,1)],
#=17=# [(9,1),(2,1),(8,1)],
],
# SpaceGroup{3}
Vector{Tuple{Int8,Int8}}[
#=1=# [(1,1)],
#=2=# [(2,1)],
#=3=# [(3,1)],
#=4=# [(3,4)],
#=5=# [(3,1),(1,5)],
#=6=# [(4,1)],
#=7=# [(4,2)],
#=8=# [(4,1),(1,5)],
#=9=# [(4,2),(1,5)],
#=10=# [(3,1),(2,1)],
#=11=# [(3,4),(2,1)],
#=12=# [(3,1),(2,1),(1,5)],
#=13=# [(3,2),(2,1)],
#=14=# [(3,7),(2,1)],
#=15=# [(3,2),(2,1),(1,5)],
#=16=# [(6,1),(3,1)],
#=17=# [(6,2),(3,2)],
#=18=# [(6,1),(3,5)],
#=19=# [(6,6),(3,7)],
#=20=# [(6,2),(3,2),(1,5)],
#=21=# [(6,1),(3,1),(1,5)],
#=22=# [(6,1),(3,1),(1,7),(1,6)],
#=23=# [(6,1),(3,1),(1,8)],
#=24=# [(6,6),(3,7),(1,8)],
#=25=# [(6,1),(4,1)],
#=26=# [(6,2),(4,2)],
#=27=# [(6,1),(4,2)],
#=28=# [(6,1),(4,3)],
#=29=# [(6,2),(4,3)],
#=30=# [(6,1),(4,7)],
#=31=# [(6,6),(4,6)],
#=32=# [(6,1),(4,5)],
#=33=# [(6,2),(4,5)],
#=34=# [(6,1),(4,8)],
#=35=# [(6,1),(4,1),(1,5)],
#=36=# [(6,2),(4,2),(1,5)],
#=37=# [(6,1),(4,2),(1,5)],
#=38=# [(6,1),(4,1),(1,7)],
#=39=# [(6,1),(4,4),(1,7)],
#=40=# [(6,1),(4,3),(1,7)],
#=41=# [(6,1),(4,5),(1,7)],
#=42=# [(6,1),(4,1),(1,7),(1,6)],
#=43=# [(6,1),(4,9),(1,7),(1,6)],
#=44=# [(6,1),(4,1),(1,8)],
#=45=# [(6,1),(4,5),(1,8)],
#=46=# [(6,1),(4,3),(1,8)],
#=47=# [(6,1),(3,1),(2,1)],
#=48=# [(6,5),(3,6),(2,1)],
#=49=# [(6,1),(3,2),(2,1)],
#=50=# [(6,5),(3,3),(2,1)],
#=51=# [(6,3),(3,1),(2,1)],
#=52=# [(6,3),(3,8),(2,1)],
#=53=# [(6,6),(3,6),(2,1)],
#=54=# [(6,3),(3,2),(2,1)],
#=55=# [(6,1),(3,5),(2,1)],
#=56=# [(6,5),(3,7),(2,1)],
#=57=# [(6,2),(3,7),(2,1)],
#=58=# [(6,1),(3,8),(2,1)],
#=59=# [(6,5),(3,4),(2,1)],
#=60=# [(6,8),(3,2),(2,1)],
#=61=# [(6,6),(3,7),(2,1)],
#=62=# [(6,6),(3,4),(2,1)],
#=63=# [(6,2),(3,2),(2,1),(1,5)],
#=64=# [(6,7),(3,7),(2,1),(1,5)],
#=65=# [(6,1),(3,1),(2,1),(1,5)],
#=66=# [(6,1),(3,2),(2,1),(1,5)],
#=67=# [(6,4),(3,4),(2,1),(1,5)],
#=68=# [(6,3),(3,2),(2,1),(1,5)],
#=69=# [(6,1),(3,1),(2,1),(1,7),(1,6)],
#=70=# [(6,13),(3,12),(2,1),(1,7),(1,6)],
#=71=# [(6,1),(3,1),(2,1),(1,8)],
#=72=# [(6,1),(3,5),(2,1),(1,8)],
#=73=# [(6,6),(3,7),(2,1),(1,8)],
#=74=# [(6,4),(3,4),(2,1),(1,8)],
#=75=# [(6,1),(9,1)],
#=76=# [(6,2),(9,23)],
#=77=# [(6,1),(9,2)],
#=78=# [(6,2),(9,24)],
#=79=# [(6,1),(9,1),(1,8)],
#=80=# [(6,8),(9,27),(1,8)],
#=81=# [(6,1),(11,1)],
#=82=# [(6,1),(11,1),(1,8)],
#=83=# [(6,1),(9,1),(2,1)],
#=84=# [(6,1),(9,2),(2,1)],
#=85=# [(6,5),(9,3),(2,1)],
#=86=# [(6,5),(9,7),(2,1)],
#=87=# [(6,1),(9,1),(2,1),(1,8)],
#=88=# [(6,6),(9,29),(2,1),(1,8)],
#=89=# [(6,1),(9,1),(3,1)],
#=90=# [(6,1),(9,5),(3,5)],
#=91=# [(6,2),(9,23),(3,1)],
#=92=# [(6,2),(9,25),(3,25)],
#=93=# [(6,1),(9,2),(3,1)],
#=94=# [(6,1),(9,8),(3,8)],
#=95=# [(6,2),(9,24),(3,1)],
#=96=# [(6,2),(9,26),(3,26)],
#=97=# [(6,1),(9,1),(3,1),(1,8)],
#=98=# [(6,8),(9,27),(3,39),(1,8)],
#=99=# [(6,1),(9,1),(4,1)],
#=100=# [(6,1),(9,1),(4,5)],
#=101=# [(6,1),(9,2),(4,2)],
#=102=# [(6,1),(9,8),(4,8)],
#=103=# [(6,1),(9,1),(4,2)],
#=104=# [(6,1),(9,1),(4,8)],
#=105=# [(6,1),(9,2),(4,1)],
#=106=# [(6,1),(9,2),(4,5)],
#=107=# [(6,1),(9,1),(4,1),(1,8)],
#=108=# [(6,1),(9,1),(4,2),(1,8)],
#=109=# [(6,8),(9,27),(4,1),(1,8)],
#=110=# [(6,8),(9,27),(4,2),(1,8)],
#=111=# [(6,1),(11,1),(3,1)],
#=112=# [(6,1),(11,1),(3,2)],
#=113=# [(6,1),(11,1),(3,5)],
#=114=# [(6,1),(11,1),(3,8)],
#=115=# [(6,1),(11,1),(4,1)],
#=116=# [(6,1),(11,1),(4,2)],
#=117=# [(6,1),(11,1),(4,5)],
#=118=# [(6,1),(11,1),(4,8)],
#=119=# [(6,1),(11,1),(4,1),(1,8)],
#=120=# [(6,1),(11,1),(4,2),(1,8)],
#=121=# [(6,1),(11,1),(3,1),(1,8)],
#=122=# [(6,1),(11,1),(3,39),(1,8)],
#=123=# [(6,1),(9,1),(3,1),(2,1)],
#=124=# [(6,1),(9,1),(3,2),(2,1)],
#=125=# [(6,5),(9,3),(3,3),(2,1)],
#=126=# [(6,5),(9,3),(3,6),(2,1)],
#=127=# [(6,1),(9,1),(3,5),(2,1)],
#=128=# [(6,1),(9,1),(3,8),(2,1)],
#=129=# [(6,5),(9,3),(3,4),(2,1)],
#=130=# [(6,5),(9,3),(3,7),(2,1)],
#=131=# [(6,1),(9,2),(3,1),(2,1)],
#=132=# [(6,1),(9,2),(3,2),(2,1)],
#=133=# [(6,5),(9,6),(3,3),(2,1)],
#=134=# [(6,5),(9,6),(3,6),(2,1)],
#=135=# [(6,1),(9,2),(3,5),(2,1)],
#=136=# [(6,1),(9,8),(3,8),(2,1)],
#=137=# [(6,5),(9,6),(3,4),(2,1)],
#=138=# [(6,5),(9,6),(3,7),(2,1)],
#=139=# [(6,1),(9,1),(3,1),(2,1),(1,8)],
#=140=# [(6,1),(9,1),(3,2),(2,1),(1,8)],
#=141=# [(6,6),(9,31),(3,6),(2,1),(1,8)],
#=142=# [(6,6),(9,31),(3,3),(2,1),(1,8)],
#=143=# [(17,1)],
#=144=# [(17,35)],
#=145=# [(17,36)],
#=146=# [(17,1),(1,53)],
#=147=# [(17,1),(2,1)],
#=148=# [(17,1),(2,1),(1,53)],
#=149=# [(17,1),(14,1)],
#=150=# [(17,1),(13,1)],
#=151=# [(17,35),(14,36)],
#=152=# [(17,35),(13,1)],
#=153=# [(17,36),(14,35)],
#=154=# [(17,36),(13,1)],
#=155=# [(17,1),(13,1),(1,53)],
#=156=# [(17,1),(15,1)],
#=157=# [(17,1),(16,1)],
#=158=# [(17,1),(15,2)],
#=159=# [(17,1),(16,2)],
#=160=# [(17,1),(15,1),(1,53)],
#=161=# [(17,1),(15,2),(1,53)],
#=162=# [(17,1),(14,1),(2,1)],
#=163=# [(17,1),(14,2),(2,1)],
#=164=# [(17,1),(13,1),(2,1)],
#=165=# [(17,1),(13,2),(2,1)],
#=166=# [(17,1),(13,1),(2,1),(1,53)],
#=167=# [(17,1),(13,2),(2,1),(1,53)],
#=168=# [(17,1),(6,1)],
#=169=# [(17,35),(6,2)],
#=170=# [(17,36),(6,2)],
#=171=# [(17,36),(6,1)],
#=172=# [(17,35),(6,1)],
#=173=# [(17,1),(6,2)],
#=174=# [(17,1),(8,1)],
#=175=# [(17,1),(6,1),(2,1)],
#=176=# [(17,1),(6,2),(2,1)],
#=177=# [(17,1),(6,1),(13,1)],
#=178=# [(17,35),(6,2),(13,35)],
#=179=# [(17,36),(6,2),(13,36)],
#=180=# [(17,36),(6,1),(13,36)],
#=181=# [(17,35),(6,1),(13,35)],
#=182=# [(17,1),(6,2),(13,1)],
#=183=# [(17,1),(6,1),(15,1)],
#=184=# [(17,1),(6,1),(15,2)],
#=185=# [(17,1),(6,2),(15,2)],
#=186=# [(17,1),(6,2),(15,1)],
#=187=# [(17,1),(8,1),(15,1)],
#=188=# [(17,1),(8,2),(15,2)],
#=189=# [(17,1),(8,1),(13,1)],
#=190=# [(17,1),(8,2),(13,1)],
#=191=# [(17,1),(6,1),(13,1),(2,1)],
#=192=# [(17,1),(6,1),(13,2),(2,1)],
#=193=# [(17,1),(6,2),(13,2),(2,1)],
#=194=# [(17,1),(6,2),(13,1),(2,1)],
#=195=# [(6,1),(3,1),(33,1)],
#=196=# [(6,1),(3,1),(33,1),(1,7),(1,6)],
#=197=# [(6,1),(3,1),(33,1),(1,8)],
#=198=# [(6,6),(3,7),(33,1)],
#=199=# [(6,6),(3,7),(33,1),(1,8)],
#=200=# [(6,1),(3,1),(33,1),(2,1)],
#=201=# [(6,5),(3,6),(33,1),(2,1)],
#=202=# [(6,1),(3,1),(33,1),(2,1),(1,7),(1,6)],
#=203=# [(6,13),(3,12),(33,1),(2,1),(1,7),(1,6)],
#=204=# [(6,1),(3,1),(33,1),(2,1),(1,8)],
#=205=# [(6,6),(3,7),(33,1),(2,1)],
#=206=# [(6,6),(3,7),(33,1),(2,1),(1,8)],
#=207=# [(6,1),(3,1),(33,1),(13,1)],
#=208=# [(6,1),(3,1),(33,1),(13,8)],
#=209=# [(6,1),(3,1),(33,1),(13,1),(1,7),(1,6)],
#=210=# [(6,7),(3,5),(33,1),(13,33),(1,7),(1,6)],
#=211=# [(6,1),(3,1),(33,1),(13,1),(1,8)],
#=212=# [(6,6),(3,7),(33,1),(13,30)],
#=213=# [(6,6),(3,7),(33,1),(13,29)],
#=214=# [(6,6),(3,7),(33,1),(13,29),(1,8)],
#=215=# [(6,1),(3,1),(33,1),(16,1)],
#=216=# [(6,1),(3,1),(33,1),(16,1),(1,7),(1,6)],
#=217=# [(6,1),(3,1),(33,1),(16,1),(1,8)],
#=218=# [(6,1),(3,1),(33,1),(16,8)],
#=219=# [(6,1),(3,1),(33,1),(16,8),(1,7),(1,6)],
#=220=# [(6,6),(3,7),(33,1),(16,9),(1,8)],
#=221=# [(6,1),(3,1),(33,1),(13,1),(2,1)],
#=222=# [(6,5),(3,6),(33,1),(13,2),(2,1)],
#=223=# [(6,1),(3,1),(33,1),(13,8),(2,1)],
#=224=# [(6,5),(3,6),(33,1),(13,5),(2,1)],
#=225=# [(6,1),(3,1),(33,1),(13,1),(2,1),(1,7),(1,6)],
#=226=# [(6,1),(3,1),(33,1),(13,8),(2,1),(1,7),(1,6)],
#=227=# [(6,41),(3,42),(33,1),(13,41),(2,1),(1,7),(1,6)],
#=228=# [(6,44),(3,45),(33,1),(13,47),(2,1),(1,7),(1,6)],
#=229=# [(6,1),(3,1),(33,1),(13,1),(2,1),(1,8)],
#=230=# [(6,6),(3,7),(33,1),(13,29),(2,1),(1,8)],
],
)