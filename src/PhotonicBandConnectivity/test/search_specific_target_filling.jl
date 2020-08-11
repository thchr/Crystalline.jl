using Crystalline
if !isdefined(Main, :(PhotonicBandConnectivity))
    includet("../src/PhotonicBandConnectivity.jl")
    using Main.PhotonicBandConnectivity
    const PBC = PhotonicBandConnectivity
end
if !isdefined(Main, :(SymmetryBases))
    includet("../../SymmetryBases/SymmetryBases.jl")
    using Main.SymmetryBases
end
using PrettyTables

#------------------------------------------------------------------------------------------
# UTILITY FUNCTIONS

function add_content_to_symvec_str_at_kidx(str::String, kidx::Integer, insert::String)

    parts = replace.(split.(strip.(str, Ref(('[', ']'))), ", "), Ref(" "=>""))
    Γslot = parts[kidx]
    parts[kidx] = insert*(length(Γslot) > 0 ? '+'*Γslot : "")

    return '['*join(parts, ", ")*']'
end

#------------------------------------------------------------------------------------------
# SETUP/PREP-WORK

#sgnum        = 230
timereversal = true
verbose      = false
latex        = false
file         = "/mnt/c/Dropbox (MIT)/Web/thchr.github.io/_assets/connectivity_tables.md"
io           = open(file, "w+")
νᵗstart      = 3 # start target filling
Nsolutions   = 2 # how many νᵀ solutions we want

for sgnum in 1:230
t=@elapsed begin
lgirs_Γ  = PBC.get_lgirreps_at_Γ(sgnum, Val(3)); # Γ-irreps
lgirsvec = get_lgirreps(sgnum, Val(3))
timereversal && (lgirs_Γ  = realify(lgirs_Γ))
timereversal && (lgirsvec = realify.(lgirsvec))

# mostly copied stuff from inimal_expansion_of_zero_freq_bands and
# find_minimum_bandreps_regular1L: should probably be factored into parts
sb, Γidxs = PBC.compatibility_bases_and_Γidxs(sgnum, lgirs_Γ, timereversal)
Nⁱʳʳ      = length(first(sb))
notΓidxs  = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]

nontopo_sb, _ = nontopological_bases(sgnum, timereversal=timereversal)
nontopo_M     = SymmetryBases.matrix(nontopo_sb)
BRS           = bandreps(sgnum, timereversal=timereversal)

ms¹ᴸ = PBC.find_representation¹ᴸ(lgirs_Γ)
ms²ᵀ = PBC.find_representation²ᵀ(lgirs_Γ)
ms   = PBC.find_representation²ᵀ⁺¹ᴸ(lgirs_Γ)

ntidxs¹ᴸ  = PBC.find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
νsᴴ       = PBC.fillings(sb)
_, pick¹ᴸ = findmin(νsᴴ[ntidxs¹ᴸ])
idx¹ᴸ     = ntidxs¹ᴸ[pick¹ᴸ]

#------------------------------------------------------------------------------------------
# FINDING VALID SOLUTIONS

νᵗ               = νᵗstart-1
Nsolutions_found = 0
Γkv_idx_in_sb    = findfirst(==("Γ"), sb.klabs)
topos =[]
println(io, "## SG $sgnum")
while true
    #global νᵗ, Nsolutions_found
    νᵗ += 1
    verbose && print("   … νᵗ = ", νᵗ, ": ");

    cⁱs, νᵀ = check_target_filling_regular1L(νᵗ, ms¹ᴸ, ms, νsᴴ, sb, idx¹ᴸ, Γidxs, notΓidxs; 
                                             verbose=verbose)
    # go to next νᵗ if no solutions found
    isempty(cⁱs) && continue 

    # Extract the "physical parts" of the symmetry vector (i.e. sans singular Γ-irreps)
    nᵀ⁺ᴸs = PBC.sum_symbases.(Ref(sb), cⁱs)
    nᴸ    = sb[idx¹ᴸ]
    nᵀs   = nᵀ⁺ᴸs .-  Ref(nᴸ)
    for nᵀ in nᵀs
        # tricky indexing: this achieves the ordering-adjusted subtraction into the `sb`
        # basis (note the broadcasted assignment, which is what makes this work)
        nᵀ[Γidxs] .-= ms²ᵀ
    end
    unique!(nᵀs) # remove equivalent solutions (due to expansion non-uniqueness)
    
    # construct human-readable symmetry vectors
    ios = [IOBuffer() for _ in eachindex(nᵀs)]
    Crystalline.prettyprint_symmetryvector.(ios, nᵀs, Ref(sb.irlabs))
    nᵀs_str = String.(take!.(ios))

    # insert (▪)²ᵀ in Γ-irrep spot
    nᵀs_str .= add_content_to_symvec_str_at_kidx.(nᵀs_str, Γkv_idx_in_sb, Ref("(▪)²ᵀ"))
    latex && (nᵀs_str .= Crystalline.convert_irreplabel2latex.(nᵀs_str))

    # find "Z₁" factor-type topology of each solution
    topos = topology_from_2T1L_xor_1L.(nᵀs, Ref(nᴸ), Ref(ms²ᵀ), Ref(Γidxs), Ref(nontopo_M))

    # write results to screen
    plural = length(nᵀs) > 1
    print(io, "\n", "νᵀ = ", νᵀ, " (", length(nᵀs), " solution", plural ? "s" : "")
    if Crystalline.classification(BRS) ≠ "Z₁"
        println(io, ")\n")
        pretty_table(io,
            [nᵀs_str topos],    # contents
            ["nᵀ", "ℤ₂ index"]; # header row
            tf = unicode,
            vlines = :none, hlines = [:begin, 1, :end],
            alignment = :l,
            crop = :none,
            #tf = markdown,
            )
    else
        println(io, "; trivial symmetry indicator group, ℤ₁)\n")
        pretty_table(io,
            nᵀs_str,            # contents
            ["nᵀ"];             # header row
            tf = unicode,
            vlines = :none, hlines = [:begin, 1, :end],
            alignment = :l,
            crop = :none,
            #tf = markdown,
            )
    end

    # check if we are done
    Nsolutions_found += 1
    Nsolutions_found ≥ Nsolutions && break
    if Nsolutions_found ≥ 1 && (sgnum ∈ (2,10,47))  # skip computational quagmires 
        println("Skipped higher-order solutions for ", sgnum)
        break
    end
end # while
println(io)
end # begin (@elapsed)
println("SG ", sgnum, ": ", round(t, digits=1), " s")
end # for sgnum
close(io)
nothing