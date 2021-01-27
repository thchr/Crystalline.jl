# Utility to write a ::BandRepSet to a .csv file
function bandrep2csv(filename::String, BRS::BandRepSet)
    open(filename, "w") do io
        bandrep2csv(io, BRS)
    end
end

function bandrep2csv(io::IO, BRS::BandRepSet)
    print(io, "Wyckoff pos.")
    for br in BRS
        print(io, "|", br.wyckpos, "(", br.sitesym, ")")
    end
    println(io)

    print(io, "Band-Rep.")
    for br in BRS
        print(io, "|", br.label, "(", br.dim, ")")
    end
    println(io)

    print(io, "Decomposable")
    for br in BRS
        print(io, "|", br.decomposable)
    end
    println(io)

    irlabs = BRS.irlabs
    for (klab,kv) in zip(BRS.klabs, BRS.kvs)
        print(io, klab, ":")
        print(io, '('*strip(replace(string(kv), " "=>""), ('[', ']'))*')')

        iridxs = findall(irlab->occursin(klab, irlab), irlabs)
        for br in BRS
            print(io, "|", )
            has_multiple = false
            for iridx in iridxs
                v = br.irvec[iridx]
                if !iszero(v)
                    has_multiple && print(io, '⊕')
                    isone(v) || print(io, v)
                    print(io, irlabs[iridx])
                    # TODO: write irrep dimension; possible, but tedious & requires
                    #       loading irreps separate
                    has_multiple = true
                end
            end
        end
        println(io)
    end
end

# Load ::BandRepSet for all 2D plane groups
include("tbd-filename.jl") # should define BRS_vec

# for allpaths in (false, true)
# for elementary in ...
for (sgnum, BRS) in enumerate(BRS_vec)
    filename = "../data/2d/elementary/allpaths/$sgnum.csv"
    bandrep2csv(filename, BRS)
end