using Crystalline, PyPlot
write_to_file = true
Nk = 20
#kvecs = get_kvpath([[0,0], [0.5,0], [0.5,0.5], [0.0,0.5], [0,0]], Nk)
kvecs = get_kvpath([[0,0], [0.5,0], [0.5,0.5], [0,0]], Nk)

PyPlot.close("all")
sgnum = 17
dim = 2
flat = levelsetlattice(sgnum, dim, ntuple(_->2, dim))
mflat = modulate(flat)
smflat = normscale(mflat,1);
Rs = directbasis(sgnum, dim)

epsin = 10.0
epsout = 1.0
res = 64
id = 1

write_dir = (@__DIR__)*"/../../mpb-ctl/input/"

# create continuously varied variations of the initial lattice, by complex rotation
smflat′ = deepcopy(smflat)
fig = plt.figure()
for (step, filling) in pairs(range(0.2, .8, length=2))
    # plot
    plot(smflat′, Rs, filling=filling, fig=fig, repeat=1)
    pause(.01)
    
    # write to file to disk
    if false
        id′ = string(id)*"-fill"*string(step)
        for runtype in ("tm", "te")
            filename = Crystalline.mpb_calcname(dim, sgnum, id′, res, runtype)
            open(write_dir*filename*".sh", "w") do io
                calcname = prepare_mpbcalc!(io, sgnum, smflat, Rs, filling, epsin, epsout, 
                                            runtype; res=res, kvecs=kvecs, id=id′)
            end
        end
    end
end