using SGOps
write_to_file = false
Nk = 20
kvecs = get_kvpath([[0,0], [0.5,0], [0.5,0.5], [0,0]], Nk)

sgnum = 1
D = 2
flat = levelsetlattice(sgnum, D, ntuple(_->2, D))
mflat = modulate(flat)
smflat = normscale(mflat,0);
C = gen_crystal(sgnum, D)
filling = .5
epsin = 10.0
epsout = 1.0
res = 32
runtype = "tm"
id = 1

if write_to_file
    write_dir = (@__DIR__)*"/../../../mpb-ctl/input/"
    filename = SGOps.mpb_calcname(D, sgnum, id, res, runtype)
    open(write_dir*filename*".sh", "w") do io
        calcname = prepare_mpbcalc!(io, sgnum, mflat, C, filling, epsin, epsout, kvecs, id, res, runtype)
    end
end

plot(smflat, C; repeat=1)