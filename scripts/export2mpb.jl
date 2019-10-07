using SGOps
write_to_file = false

sgnum = 4
dim = 2
flat = levelsetlattice(sgnum, dim, ntuple(_->2, dim))
mflat = modulate(flat)
smflat = normscale(mflat,0);
C = gen_crystal(sgnum, dim)
filling = .5
epsin = 10.0
epsout = 1.0

if write_to_file
    write_dir = (@__DIR__)*"/../../../mpb-ctl/"
    open(write_dir*"test_lattice.sh", "w") do io
        prepare_mpbcalc!(io, sgnum, mflat, C, filling, epsin, epsout) 
    end
end

plot(smflat, C, repeat=1)
# TODO: write to a format loadable by Julia (write sgnum, mflat, C, isoval, epsin, epsout) later on...