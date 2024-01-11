# generated from original `xyzt`-based data via:
#   begin
#   io = stdout
#   println(io, "PG_GENS_CODES_Ds = (")
#   for D in 1:3
#       pglabs = Crystalline.PG_IUCs[D]
#       gensv = generators_from_xyzt.(pglabs, PointGroup{D})
#       ROTATIONS = D == 3 ? Crystalline.ROTATIONS_3D : 
#                   D == 2 ? Crystalline.ROTATIONS_2D :
#                   D == 1 ? Crystalline.ROTATIONS_1D : error()
#       println(io, "# PointGroup{", D, "}")
#       println(io, "Dict{String,Vector{Int8}}(")
#       for (pglab, gens) in zip(pglabs, gensv)
#           Nops = length(gens)
#           pg_codes = Vector{Int8}(undef, Nops)
#           idx = 0
#           for (n,op) in enumerate(gens)
#               r = rotation(op)
#               r_idx = Int8(something(findfirst(==(r), ROTATIONS)))
#               pg_codes[idx+=1] = r_idx
#           end       
#           print(io, "\"", pglab, "\" => [")
#           join(io, pg_codes, ",")
#           println(io, "],")
#       end
#       println(io, "),")
#   end
#   println(io, ")")
#   end

const PG_GENS_CODES_Ds = (
# PointGroup{1}
Dict{String,Vector{Int8}}(
"1" => [1],
"m" => [2],
),
# PointGroup{2}
Dict{String,Vector{Int8}}(
"1" => [1],
"2" => [2],
"m" => [3],
"mm2" => [2,3],
"4" => [2,5],
"4mm" => [2,5,3],
"3" => [9],
"3m1" => [9,8],
"31m" => [9,7],
"6" => [9,2],
"6mm" => [9,2,8],
),
# PointGroup{3}
Dict{String,Vector{Int8}}(
"1" => [1],
"-1" => [2],
"2" => [3],
"m" => [4],
"2/m" => [3,2],
"222" => [6,3],
"mm2" => [6,4],
"mmm" => [6,3,2],
"4" => [6,9],
"-4" => [6,11],
"4/m" => [6,9,2],
"422" => [6,9,3],
"4mm" => [6,9,4],
"-42m" => [6,11,3],
"-4m2" => [6,11,4],
"4/mmm" => [6,9,3,2],
"3" => [17],
"-3" => [17,2],
"312" => [17,14],
"321" => [17,13],
"3m1" => [17,15],
"31m" => [17,16],
"-31m" => [17,14,2],
"-3m1" => [17,13,2],
"6" => [17,6],
"-6" => [17,8],
"6/m" => [17,6,2],
"622" => [17,6,13],
"6mm" => [17,6,15],
"-62m" => [17,8,13],
"-6m2" => [17,8,15],
"6/mmm" => [17,6,13,2],
"23" => [6,3,33],
"m-3" => [6,3,33,2],
"432" => [6,3,33,13],
"-43m" => [6,3,33,16],
"m-3m" => [6,3,33,13,2],
)
)