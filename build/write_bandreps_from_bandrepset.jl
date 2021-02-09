# Utility to write a ::BandRepSet to a .csv file
using Crystalline

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

# ---------------------------------------------------------------------------------------- #
# define `make_bandrep_set` that creates `::BandRepSet`
include("setup_2d_band_representations.jl")

basepath = joinpath((@__DIR__), "..", "data/bandreps/2d")
for allpaths in (false, true)
    path_tag = allpaths ? "allpaths" : "maxpaths"
    for timereversal in (false, true)
        tr_tag = timereversal ? "elementaryTR" : "elementary"
        for sgnum in 1:MAX_SGNUM[2]
            BRS = make_bandrep_set(sgnum, Val(2); allpaths=allpaths,
                                                  timereversal=timereversal)
            filename = joinpath(basepath, tr_tag, path_tag, string(sgnum)*".csv")
            bandrep2csv(filename, BRS)
        end
    end
end