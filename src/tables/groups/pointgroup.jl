# generated from original `xyzt`-based data via:
#    PG_CODES_3D_D = Dict{String,Vector{Int8}}()
#    pglabs = Crystalline.PG_IUCs[3]
#    pgs = pointgroup.(pglabs, Val(3))
#    for pg in pgs
#        pglab = label(pg)
#        Nops = length(pg)
#        pg_codes = Vector{Int8}(undef, Nops-1)
#        idx = 0
#        for (n,op) in enumerate(pg)
#            isone(op) && continue # skip identity (trivially in  group): no need to list
#            r = rotation(op)
#            r_idx = Int8(something(findfirst(==(r), Crystalline.ROTATIONS_3D)))
#            pg_codes[idx+=1] = r_idx
#        end
#        PG_CODES_3D_D[pglab] = pg_codes
#    end
#    for pglab in pglabs
#        println(pglab => PG_CODES_3D_D[pglab], ",")
#    end

const PG_CODES_Ds = (
# PointGroup{1}
Dict{String, Vector{Int8}}(
"1"   => [],
"m"   => [2],
),
# PointGroup{2}
Dict{String, Vector{Int8}}(
"1"   => [],
"2"   => [2],
"m"   => [3],
"mm2" => [2,4,3],
"4"   => [2,5,6],
"4mm" => [2,5,6,4,3,8,7],
"3"   => [9,10],
"3m1" => [9,10,8,11,12],
"31m" => [9,10,7,13,14],
"6"   => [9,10,2,15,16],
"6mm" => [9,10,2,15,16,8,11,12,7,13,14],
),
# PointGroup{3}
Dict{String, Vector{Int8}}(
"1"     => [],
"-1"    => [2],
"2"     => [3],
"m"     => [4],
"2/m"   => [3,2,4],
"222"   => [6,3,5],
"mm2"   => [6,4,7],
"mmm"   => [6,3,5,2,8,4,7],
"4"     => [6,9,10],
"-4"    => [6,11,12],
"4/m"   => [6,9,10,2,8,11,12],
"422"   => [6,9,10,3,5,13,14],
"4mm"   => [6,9,10,4,7,15,16],
"-42m"  => [6,11,12,3,5,15,16],
"-4m2"  => [6,11,12,4,7,13,14],
"4/mmm" => [6,9,10,3,5,13,14,2,8,11,12,4,7,15,16],
"3"     => [17,18],
"-3"    => [17,18,2,19,20],
"312"   => [17,18,14,22,21],
"321"   => [17,18,13,23,24],
"3m1"   => [17,18,15,25,26],
"31m"   => [17,18,16,28,27],
"-31m"  => [17,18,14,22,21,2,19,20,16,28,27],
"-3m1"  => [17,18,13,23,24,2,19,20,15,25,26],
"6"     => [17,18,6,30,29],
"-6"    => [17,18,8,32,31],
"6/m"   => [17,18,6,30,29,2,19,20,8,32,31],
"622"   => [17,18,6,30,29,13,23,24,14,22,21],
"6mm"   => [17,18,6,30,29,15,25,26,16,28,27],
"-62m"  => [17,18,8,32,31,13,23,24,16,28,27],
"-6m2"  => [17,18,8,32,31,15,25,26,14,22,21],
"6/mmm" => [17,18,6,30,29,13,23,24,14,22,21,2,19,20,8,32,31,15,25,26,16,28,27],
"23"    => [6,3,5,33,38,36,40,34,35,39,37],
"m-3"   => [6,3,5,33,38,36,40,34,35,39,37,2,8,4,7,41,46,44,48,42,43,47,45],
"432"   => [6,3,5,33,38,36,40,34,35,39,37,13,14,10,9,50,53,54,49,51,55,52,56],
"-43m"  => [6,3,5,33,38,36,40,34,35,39,37,16,15,11,12,62,57,58,61,64,60,63,59],
"m-3m"  => [6,3,5,33,38,36,40,34,35,39,37,13,14,10,9,50,53,54,49,51,55,52,56,2,8,4,7,41,46,44,48,42,43,47,45,15,16,12,11,58,61,62,57,59,63,60,64],
)
)