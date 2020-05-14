using Crystalline
write_to_file = false
Nk = 20
kvecs = get_kvpath([[0,0], [0.5,0], [0.5,0.5], [0,0]], Nk)

sgnum = 1
D = 2
flat = levelsetlattice(sgnum, D, ntuple(_->2, D))
mflat = modulate(flat)
smflat = normscale(mflat,0);
Rs = directbasis(sgnum, D)
filling = .5
epsin = 10.0
epsout = 1.0
res = 32
runtype = "tm"
id = 1

if write_to_file
    write_dir = (@__DIR__)*"/../../mpb-ctl/input/"
    filename = Crystalline.mpb_calcname(D, sgnum, id, res, runtype)
    open(write_dir*filename*".sh", "w") do io
        prepare_mpbcalc!(io, sgnum, mflat, Rs, filling, epsin, epsout, runtype;
                             res=res, kvecs=kvecs, id=id)
    end
end

plot(smflat, Rs; repeat=1)