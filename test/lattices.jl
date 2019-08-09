using SGOps, PyPlot, AbstractPlotting, MAT

# 2D lattices
plt.close("all")
sgnum = 1; dim = 2; idxmax = (2,2)
N = 100 # number of points per dimension
expon = .5 # exponent for fall off in coefficients
filling = 0.5

ijkorbits, orbcoefs, R = levelsetlattice(sgnum, dim, idxmax)
plotfourier(ijkorbits, orbcoefs, R, N, expon, filling)
display(plt.gcf().show());
#display(fig.show()) # force focus when run from script

if true
    save2matlab = true
    # 3D lattices
    sgnum = 139; dim = 3; idxmax = 2 .*(1,1,1)
    N = 50 # number of points per dimension
    expon = 2 # exponent for fall off in coefficients
    filling = 0.5

    ijkorbits, orbcoefs, R = levelsetlattice(sgnum, dim, idxmax)
    xyz,vals,isoval=plotfourier(ijkorbits, orbcoefs, R, N, expon, filling)
    #display(AbstractPlotting.current_scene());

    if save2matlab
        file = matopen("test/isosurfvis/isosurfdata.mat", "w")
        write(file, "isoval", isoval)
        write(file, "xyz", collect(xyz))
        write(file, "vals", vals)
        close(file)
    end
end