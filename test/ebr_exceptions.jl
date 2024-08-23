# general exceptions (Table 3 of https://doi.org/10.1103/PhysRevB.97.035139)
const ebr_exceptions = Dict{Int, Tuple{Bool, String, Vector{String}}}(
    # sgnum => (applies_with_tr, site group (Schoenflies notation), irreps (Mulliken notation))
    ([163, 165, 167, 228, 230, 223, 211, 208, 210, 228, 188, 190, 192, 193] 
            .=> Ref((false, "D₃"  #=312/321=#,    ["E"] #=Γ₃=#))) ...,
      192    =>    ((false, "D₆"  #=622=#,        ["E₁", "E₂"] #=Γ₅, Γ₆=#)),
    ([207, 211, 222, 124, 140] 
            .=> Ref((false, "D₄"  #=422=#,        ["E"] #=Γ₅=#))) ...,
    ([229, 226, 215, 217, 224, 131, 132, 139, 140, 223] 
            .=> Ref((true, "D₂d" #=-42m/-4m2=#,  ["E"] #=Γ₅=#))) ...
)

# time-reversal exceptions (Table 2 of https://doi.org/10.1103/PhysRevB.97.035139)
# always associated with site group S₄ (-4) and site irrep E (Γ₃Γ₄)
const ebr_tr_exceptions = (84, 87, 135, 136, 112, 116, 120, 121, 126, 130, 133, 138, 142,
                           218, 230, 222, 217, 219, 228)

# ---------------------------------------------------------------------------------------- #

function is_exceptional_br(sgnum::Integer, br::BandRep; timereversal::Bool=true)
    # first check (Table 2): a composite physically real complex bandrep ("realified" rep)?
    (timereversal && is_exceptional_tr_br(sgnum, br)) && return true      

    # second check (Table 3)
    haskey(ebr_exceptions, sgnum) || return false
    except_applies_with_tr, except_sitesym, except_siteirlabs = ebr_exceptions[sgnum]

    # if we have TR, table 3 entries only apply if `except_applies_with_tr = true`
    timereversal && (except_applies_with_tr || return false)
    
    # check if site symmetry group matches
    br_sitesym_iuc = br.sitesym == "32" ? "321" : br.sitesym
    br_sitesym_schoenflies = Crystalline.PG_IUC2SCHOENFLIES[br_sitesym_iuc]
    except_sitesym == br_sitesym_schoenflies || return false # didn't match

    # check if site irrep matches
    br_siteirlab = replace(br.label, "↑G"=>"")
    br_siteirlab ∈ except_siteirlabs || return false # didn't match
    
    return true
end

function is_exceptional_tr_br(sgnum::Integer, br::BandRep)
    sgnum ∈ ebr_tr_exceptions || return false

    # check if site symmetry group is S₄
    br_sitesym_iuc = br.sitesym == "32" ? "321" : br.sitesym
    br_sitesym_schoenflies = Crystalline.PG_IUC2SCHOENFLIES[br_sitesym_iuc]
    br_sitesym_schoenflies == "S₄" || return false

    # check if site irrep is E
    br_siteirlab = replace(br.label, "↑G"=>"")
    br_siteirlab == "E" || return false

    return true
end