using SGOps, PyPlot, AbstractPlotting, MAT

# 2D lattices
plt.close("all")
sgnum = 4; dim = 2; idxmax = (2,2)
N = 100 # number of points per dimension
expon = .5 # exponent for coefficient fall-off
filling = 0.5
repeat = 1

for sgnum = 1:17
    flat = levelsetlattice(sgnum, dim, idxmax)
    C = gen_crystal(sgnum, dim) 
    for n = 1:3
        mflat = normscale!(modulate(flat), expon)
        plotfourier(mflat, C, N, filling, repeat)
        display(plt.gcf().show());
        if false
            PyPlot.savefig((@__DIR__)*"\\figures\\planegroup$(sgnum)_$(n).png")
        end
    end
end
#display(fig.show()) # force focus when run from script

if false
    save2matlab = true
    # 3D lattices
    sgnum = 139; dim = 3; idxmax = 2 .* (1,1,1)
    N = 50 # number of points per dimension
    expon = 2 # exponent for coefficient fall-off
    filling = 0.5

    C = gen_crystal(sgnum, dim)
    flat = levelsetlattice(sgnum, dim, idxmax)
    mflat = normscale!(modulate(flat), expon)
    xyz,vals,isoval=plotfourier(mflat, C, N, filling)
    #display(AbstractPlotting.current_scene());

    if save2matlab
        file = matopen("test/isosurfvis/isosurfdata.mat", "w")
        write(file, "isoval", isoval)
        write(file, "xyz", collect(xyz))
        write(file, "vals", vals)
        close(file)
    end
end