using Crystalline, Test, LinearAlgebra

if !isdefined(Main, :LGIRS) # load complex little groups, if not already loaded
    LGIRS = get_all_lgirreps(Val(3)) # ≡ parselittlegroupirreps()
end

# this file tests the consistency between little group irreps at the Γ point and the 
# associated point group irreps of the isogonal point group (they should agree, up
# basis change transformations of the irreps; i.e. their characters must agree)

for lgirsd in LGIRS
    # grab the lgirreps at Γ point; then find associated point group of sgnum
    # [pgirreps are sorted by label; ISOTROPY is not (so we sort here)]
    lgirsΓ = sort(lgirsd["Γ"], by=label)
    lg     = group(first(lgirsΓ))
    sgnum  = num(lg)
    pg     = Crystalline.find_parent_pointgroup(lg)

    pgirs = get_pgirreps(label(pg), Val(3))

    # the little group irreps at Γ should equal the irreps of the associated point group 
    # (isogonal point group or little cogroup) of that space group because Γ implies k=0,
    # which in turn turns all the translation phases {E|t} = e^{-ikt}I = I trivial. In other
    # words, at Γ, the translation parts of the little group does not matter, only the 
    # rotational parts - and the associated irreps are then just that of the point group

    # --- test that we get the same (formatted) labels ---
    @test label.(pgirs) == Crystalline.formatirreplabel.(label.(lgirsΓ))

    # --- test that we get the same (point group) operations ---
    pgops = operations(pg)
    lgops = operations(lg)
    lgops_pgpart = SymOperation{3}.(hcat.(rotation.(lgops), Ref(zeros(3))))
    
    pg_sortperm = sortperm(pgops, by=seitz) # make sure everything's sorted identically
    pgops = pgops[pg_sortperm]

    lg_sortperm = sortperm(lgops_pgpart, by=seitz) # (we sort by the _point group_ part)
    lgops = lgops[lg_sortperm]
    lgops_pgpart = lgops_pgpart[lg_sortperm]

    # --- test that pgops and lgops contain the same (point-group-part) operators ---
    @test seitz.(pgops) == seitz.(lgops_pgpart)

    # --- test that the character tables are the same ---

    # special casing for SGs 38:41 (see details under (*)):
    if sgnum ∈ 38:41
        # overall, as it relates to our tests below, this effectively swaps the irreps of 
        # the m₁₀₀ and m₀₁₀ operators
        m₁₀₀_idx = findfirst(op->seitz(op)=="m₁₀₀", lgops_pgpart)
        m₀₁₀_idx = findfirst(op->seitz(op)=="m₀₁₀", lgops_pgpart)
        a, b = lg_sortperm[m₁₀₀_idx], lg_sortperm[m₀₁₀_idx]
        lg_sortperm[m₁₀₀_idx], lg_sortperm[m₀₁₀_idx] = b, a
    end

    for (lgir, pgir) in zip(lgirsΓ, pgirs)
        @test characters(pgir)[pg_sortperm] ≈ characters(lgir)[lg_sortperm]
    end
end

# (*) Special casing for SGs 38:41 w/ centering 'A':
# These four space groups are the only one-face centered lattices with centering 'A',
# and all have the isogonal point group mm2. There are several other SGs with isogonal
# point group mm2 (18 others), but they don't need special casing.
# For these four SGs, it seems the CDML irrep-labelling convention is to swap Γ₃ and Γ₄
# relative to the point group convention, or, equivalently, to match up with a
# differently oriented point group that has m₁₀₀ and m₀₁₀ swapped. I imagine that this
# is motivated by the fact that 'A' centered space groups in principle could be recast
# as 'C' centered groups by a basis change. Long story short, we just swap m₁₀₀ and m₀₁₀
# for the little group irreps here, just to get by. I checked the LGIrreps extracted from
# ISOTROPY indeed match those of Bilbao's (REPRES and REPRESENTATIONS SG), which they do:
# they also have this "odd" difference between the Γ-point irreps for SGs 38:41, versus all
# the other irreps with isogonal point group mm2 (e.g. 37).