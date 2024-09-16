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