using Crystalline, PyPlot, AbstractPlotting, MAT

# 2D lattices
if false
    plt.close("all")
    sgnum = 4; dim = 2; idxmax = (2,2)
    N = 100 # number of points per dimension
    expon = .5 # exponent for coefficient fall-off
    filling = 0.5
    repeat = 1

    for sgnum = 1:17
        flat = levelsetlattice(sgnum, dim, idxmax)
        Rs = directbasis(sgnum, dim) 
        for n = 1:3
            mflat = normscale!(modulate(flat), expon)
            plot(mflat, Rs; N=N, filling=filling, repeat=repeat)
            display(plt.gcf().show());
            if false
                PyPlot.savefig((@__DIR__)*"\\figures\\planegroup$(sgnum)_$(n).png")
            end
        end
    end
end
#display(fig.show()) # force focus when run from script


# 3D lattices
if true
    save2matlab = false
    sgnum = 139; dim = 3; idxmax = 2 .* (1,1,1)
    N = 50 # number of points per dimension
    expon = .5 # exponent for coefficient fall-off
    filling = 0.5

    Rs = directbasis(sgnum, dim)
    flat = levelsetlattice(sgnum, dim, idxmax)
    mflat = normscale!(modulate(flat), expon)
    xyz,vals,isoval=plot(mflat, Rs; N=N, filling=filling)
    #display(AbstractPlotting.current_scene());

    if save2matlab
        file = matopen((@__DIR__)*"/isosurfvis/isosurfdata.mat", "w")
        write(file, "isoval", isoval)
        write(file, "xyz", collect(xyz))
        write(file, "vals", vals)
        close(file)
    end
end