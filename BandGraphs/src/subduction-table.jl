function SubductionTable(
            c::Connection{D},
            sg::SpaceGroup{D},
            lgirsd::AbstractDict{String, <:AbstractVector{LGIrrep{D}}}
            ) where D

    kᴳ, kᴴ = c.kᴳ, c.kᴴ
    bare_kᴳ_label = rstrip(string(kᴳ.label), '′') # w/o monodromy indicator
    lgirsᴳ = lgirsd[bare_kᴳ_label]

    # determine the free parameters for which kᴳ and kᴴ are compatible/"intersect"
    cntr = centering(sg)
    bool, αβγ, αβγ′, G = Crystalline.can_intersect(primitivize(kᴳ.kv, cntr), 
                                                   primitivize(kᴴ.kv, cntr))
    bool || error("failed to prove compatibility of $(kᴳ) and $(kᴴ)")
    iszero(αβγ) || error("implementation cannot handle nonzero free parameters for G-group")

    lgirsᴴ = lgirsd[string(kᴴ.label)]
    cosetsᴴ = cosets(sg, group(first(lgirsᴴ)))
    lgirsᴴ = Crystalline.remap_lgirreps_to_point_in_kstar(lgirsᴴ, kᴴ.kv, cosetsᴴ)
    
    # compute subduction table of G- and H-group irreps at their αβ′-intersection
    table = [subduction_count(lgirᴳ, lgirᴴ, αβγ′) for lgirᴳ in lgirsᴳ, lgirᴴ in lgirsᴴ]
    
    monodromy = endswith(string(kᴳ.label), '′') ? true : false
    irlabsᴳ, irlabsᴴ = label.(lgirsᴳ), label.(lgirsᴴ)
    if monodromy
        # append '′' to monodromy-affected irrep-labels (subduction tables are already
        # right, since αβγ′ captures the difference from the standard "in-cell" case)
        # our convention is to add `′` to every occurence of the k-label; e.g.,
        # mapping as: "X₁ → X′₁", "X₁X₂ → X′₁X′₂", and "X₁⁺X₂⁻ → X′₁⁺X′₂⁻"
        irlabsᴳ = map(irlabsᴳ) do irlab
            io = IOBuffer()
            start = firstindex(irlab)
            while (r = findnext(bare_kᴳ_label, irlab, start)) !== nothing
                write(io, irlab[start:r[end]], "′")
                start = nextind(irlab, r[end])
            end
            start == firstindex(irlab) && error("missing irrep label-change!")
            write(io, irlab[start:end])
            String(take!(io))
        end
    end

    return SubductionTable{D}(num(sg), c, irlabsᴳ, irlabsᴴ, table, monodromy)
end