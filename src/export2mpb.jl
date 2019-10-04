# probably ought to do this with show/repr and a custom mime type to be idiomatic
function lattice2mpb(flat::AbstractFourierLattice)
    orbits = getorbits(flat); coefs = getcoefs(flat)
    Nterms = sum(length, coefs)
    orbits_mpb_vec = Vector{String}(undef, Nterms)
    coefs_mpb_vec  = Vector{String}(undef, Nterms)
    idx = 1
    for (G, c) in zip(Iterators.flatten(orbits), Iterators.flatten(coefs))
        orbits_mpb_vec[idx] = "(vector3"*mapfoldl(x-> " "*string(x), *, G)*")"
        c_re = real(c); c_im = imag(c);
        coefs_mpb_vec[idx]  = string(real(c))*signaschar(c_im)*string(abs(c_im))*"i"
        idx += 1
    end
    orbits_mpb = "(list"*mapfoldl(x->" "*x, *, orbits_mpb_vec)*")"
    coefs_mpb  = "(list"*mapfoldl(x->" "*x, *, coefs_mpb_vec)*")"
    return orbits_mpb, coefs_mpb
end

function filling2isoval(flat::AbstractFourierLattice{D}, filling::Real=0.5, nsamples::Int64=51) where D
    step = 1.0/nsamples
    samples = range(-0.5, 0.5-step, length=nsamples)
    if D == 2
        itr = (real(calcfourier((x,y), flat)) for x in samples for y in samples)
    elseif D == 3
        itr = (real(calcfourier((x,y,z), flat)) for x in samples for y in samples for z in samples)
    end
    return quantile(itr, filling)
end

function crystal2mpb(C::Crystal)
    io = IOBuffer()
    write(io, "(list")
    for R in basis(C)
        write(io, " (vector3 ")
        join(io, R, ' ')
        write(io, ')')
    end
    write(io, ')')
    return String(take!(io))
end

function mpb_calcname!(io, dim, sgnum, id, res, runtype="all")
    write(io, "dim",  string(dim),
              "-sg",  string(sgnum), 
              "-",    string(id),
              "-res", string(res))
    if runtype != "all"
        write(io, "-", runtype)
    end
    return nothing
end


function prepare_mpbcalc!(io::IO, sgnum::Integer, flat::AbstractFourierLattice{D}, C::Crystal{D}, 
                                  filling::Real=0.5, εin::Real=10.0, εout::Real=1.0, 
                                  id=1, res::Integer=32, runtype::String="all") where D

    # --- prep-work to define some shell variables for logs etc ---
    write(io, "calcname=\"")
    mpb_calcname!(io, D, sgnum, id, res, runtype)
    write(io, "\"\n\n")

    # --- work to actually call mpb ---
    call_folder = "~/postdoc/mpb-transform-dev/1.8-dev/bin/mpb"
    rvecs = crystal2mpb(C)
    uc_gvecs, uc_coefs = lattice2mpb(flat)
    uc_level = filling2isoval(flat, filling)
    # kvecs = kvecs2mpb(...)     TODO
    # symops = symops2mpb(...)   TODO

    # prepare all mpb param-inputs in a single tuple
    input_tuple = (# where we call mpb from (should maybe be a function input...)
                   call_folder, 
                   # run-type ("all", "te", or "tm")
                   "run-type="*"'\""*runtype*"\"'",
                   # dimension, space group, and prefix name
                   "dim="*string(D), "sgnum="*string(sgnum), "prefix='\"'\${calcname}'\"'",
                   # crystal (basis vectors)
                   "rvecs=\""*rvecs*"\"",  
                   # unitcell/lattice shape
                   "uc-gvecs=\""*uc_gvecs*"\"", "uc-coefs=\""*uc_coefs*"\"", 
                   "uc-level="*string(uc_level),
                   # permittivities
                   "epsin="*string(εin), "epsout="*string(εout),
                   # TODO: k-points
                   # "kvecs=\""*kvecs*"\"",
                   # TODO: little group symmetry operations at each k-point (a list of lists)
                   # "symops=\""*symops*"\""
                  ) 
    # write inputs to io, separated by newlines and tabs
    join(io, input_tuple, " \\\n\t\t")
    # position of .ctl file to execute & where to write log file (TODO)
    write(io, " \\\n\t", "ctl/fourier-lattice.ctl 2>&1 | tee logs/\${calcname}.log")
    return nothing
end

function prepare_mpbcalc(sgnum::Integer, flat::AbstractFourierLattice{D}, C::Crystal{D}, 
                  filling::Real, epsin::Real=10.0, epsout::Real=1.0) where D
    io = IOBuffer()
    prepare_mpbcalc!(io, sgnum, flat, C, filling, epsin, epsout)
    return String(take!(io))
end
