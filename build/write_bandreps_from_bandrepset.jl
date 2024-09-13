# Utility to write a ::BandRepSet to a .csv file
using Crystalline

function bandrep2csv(filename::String, BRS::BandRepSet)
    open(filename, "w") do io
        bandrep2csv(io, BRS)
    end
end

function extract_klabel(s::String)
    # cheat a bit here and exploit that we know that the only greek letters
    # occuring in the k-label are Γ, Λ, Δ, Σ, and Ω
    stopidx = findfirst(c->c∉(('A':'Z'..., 'Γ', 'Λ', 'Δ', 'Σ', 'Ω')), s)::Int
    
    return s[1:prevind(s, stopidx)] # `prevind` needed since labels are not just ascii
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

    irlabs = BRS.irlabs
    for (idx, (klab,kv)) in enumerate(zip(BRS.klabs, BRS.kvs))
        print(io, klab, ":")
        print(io, '('*strip(replace(string(kv), " "=>""), ('[', ']'))*')')

        iridxs = findall(irlab -> extract_klabel(irlab) == klab, irlabs)
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
                    #       loading irreps separately
                    has_multiple = true
                end
            end
        end
        idx ≠ length(BRS.klabs) && println(io)
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
            BRS = calc_bandreps(sgnum, Val(2); allpaths=allpaths, timereversal=timereversal)
            filename = joinpath(basepath, tr_tag, path_tag, string(sgnum)*".csv")
            bandrep2csv(filename, BRS)
        end
    end
end