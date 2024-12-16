using Crystalline, Test

if !isdefined(Main, :LGIRS)
    LGIRS = lgirreps.(1:MAX_SGNUM[3], Val(3)) # loaded from our saved .jld2 files
end

@testset "k-vectors required by BandRepSet analysis" begin
allpaths = false
spinful  = false
debug = false

@testset "Complex (no TR) irreps" begin
# --- test complex-form irreps (not assuming time-reversal symmetry) ---
for (sgnum, lgirsd) in enumerate(LGIRS)
    brs = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=false)
    irlabs_brs = irreplabels(brs)
    klabs_brs = klabels(brs)

    irlabs_ISO = [label(lgir) for lgirs in values(lgirsd) for lgir in lgirs]
    klabs_ISO = keys(lgirsd)

    for (iridx_brs, irlab_brs) in enumerate(irlabs_brs)
        klab = klabel(irlab_brs)
        @test irlab_brs ∈ irlabs_ISO
        if debug && irlab_brs ∉ irlabs_ISO
            @info "Cannot find complex irrep $(irlab_brs) in ISOTROPY dataset (sgnum = $sgnum)"
        end

        # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
        kidx_brs = findfirst(==(klab), klabs_brs)
        @test brs.kvs[kidx_brs] == position(first(lgirsd[klab]))
        if debug && brs.kvs[kidx_brs] ≠ position(first(lgirsd[klab]))
            println("Different definitions of k-point labels in space group ", sgnum)
            println("   brs, ", klab, ": ",string(brs.kvs[kidx_brs]))
            println("   ISO, ", klab, ": ",string(position(first(lgirsd[klab]))), "\n")
        end

    end
end
end

@testset "Physically irreducible irreps/co-reps (with TR)" begin
# --- test physically irreducible irreps/co-reps (assuming time-reversal symmetry) ---
for (sgnum, lgirsd) in enumerate(LGIRS)
    brs = bandreps(sgnum, allpaths=allpaths, spinful=spinful, timereversal=true)
    irlabs_brs = irreplabels(brs)
    klabs_brs = klabels(brs)

    irlabs_ISO = Vector{String}()
    realirlabs_ISO = Vector{String}()
    klabs_ISO = keys(lgirsd)
    klabs_ISO = Vector{String}(undef, length(lgirsd))
    for lgirs in values(lgirsd)
        append!(irlabs_ISO,     [label(lgir) for lgir in lgirs])
        append!(realirlabs_ISO, label.(realify(lgirs)))
    end
    irlabs_ISO = irlabs_ISO
    realirlabs_ISO = realirlabs_ISO

    for (iridx_brs, irlab_brs) in enumerate(irlabs_brs)
        klab = klabel(irlab_brs)
        @test irlab_brs ∈ realirlabs_ISO
        if debug && irlab_brs ∉ realirlabs_ISO
            @info "Cannot find real irrep $(irlab_brs) in ISOTROPY dataset (sgnum = $sgnum)"
        end

        # test that ISOTROPY's labelling & representation of k-vectors agree with BCD
        kidx_brs = findfirst(==(klab), klabs_brs)
        @test brs.kvs[kidx_brs] == position(first(lgirsd[klab]))
    end
end
end
end

@testset "BandRepSet and BandRep" begin
    brs = bandreps(230)
    # iterated concatenation of vectors of `brs` should give `matrix`
    @test stack(brs) == stack(brs) == hcat(brs...)
    # length of BandRep as vectors should be = number of irreps + 1 (i.e. includes filling)
    @test length(brs[1]) == length(brs[1].irvec)+1
    @test brs[1] == vcat(brs[1].irvec, dim(brs[1]))
end


# NOTE/TODO: This would be nice to have, but is presently broken because the irrep labels
#   used by BANDREP are not quite the proper ones: e.g., BANDREP will call the 2D
#   glued-together" irrep E the ¹E²E irrep; that's not wrong per se, but it's not quite
#   right either.
#   Similarly, when there's only one A irrep, BANDREP will still include a redundant 
#   ₁-subscript; again, not wrong, but not quite right. The right fix seems to be to go
#   through the stored data we retrieve from BANDREP and then fix it there - but that's 
#   too annoying for now - so, we just don't test it at the moment.
#=
@testset "BandRepSet site-symmetry irreps" begin
    siteir_name(br) = replace(br.label, "↑G"=>"")
    for timereversal in (true)
        for sgnum in 1:230
            brs = bandreps(sgnum, 3; timereversal)
            wps = wyckoffs(sgnum)
            sitegd = Dict(label(wp)=>sitegroup(brs.sgnum, wp) for wp in wps)
            siteirsd = Dict(wp_str=>Crystalline.siteirreps(siteg) for (wp_str, siteg) in sitegd)
            timereversal && (siteirsd = Dict(wp_str => realify(siteirs) for (wp_str, siteirs) in siteirsd))
            for br in brs
                siteirs = siteirsd[br.wyckpos] 
                siteirs_labs = mulliken.(siteirs)
                @test siteir_name(br) ∈ siteirs_labs
            end
        end
    end
end
=#