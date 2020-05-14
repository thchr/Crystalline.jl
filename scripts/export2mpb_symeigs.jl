cd((@__DIR__))

sgnum = 16
D = 3 
id = "test_symeigs"
res = 8
write_dir = (@__DIR__)*"/../../mpb-ctl/input/"
filename = Crystalline.mpb_calcname(D, sgnum, id, res, "all")
open(write_dir*filename*".sh", "w") do io
    str = Crystalline.gen_symeig_mpbcalc(sgnum, D; res=res, id=id)
    write(io, str)
end