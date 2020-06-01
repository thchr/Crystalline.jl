using Crystalline, PyPlot, MAT, LinearAlgebra

# -----------------------------------------------------------------------------------------
# 2D lattices
plt.close("all")
sgnum = 4; dim = 2; idxmax = (2,2)
N = 100 # number of points per dimension
expon = .5 # exponent for coefficient fall-off
filling = 0.5
repeat = 1
savefigs = false

for sgnum in 1:MAX_SGNUM[2]
    flat = levelsetlattice(sgnum, dim, idxmax)
    Rs   = directbasis(sgnum, dim) 
    for n = 1:3
        mflat = normscale!(modulate(flat), expon)
        plot(mflat, Rs; N=N, filling=filling, repeat=repeat)
        display(plt.gcf().show());
        if savefigs
            PyPlot.savefig((@__DIR__)*"\\figures\\planegroup$(sgnum)_$(n).png")
        end
    end
end