# Corrections to data from `bandpaths.pl` on BCS
# ---------------------------------------------------------------------------------------- #
function subduction_table_corrections!(SUBDUCTIONSD_TR_3D, SUBDUCTIONSD_3D)
    subduction_table_corrections_sg200!(SUBDUCTIONSD_TR_3D, SUBDUCTIONSD_3D)
    subduction_table_corrections_sg110_tr!(SUBDUCTIONSD_TR_3D)
    subduction_table_corrections_sg110_notr!(SUBDUCTIONSD_3D)
end

# ---------------------------------------------------------------------------------------- #
# space group 200: missing ZA connections from X to M & Γ
function subduction_table_corrections_sg200!(SUBDUCTIONSD_TR_3D, SUBDUCTIONSD_3D)
    sgnum = 200
    subt_X_ZA = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:X, KVec{3}("0,1/2,0")), # actual connection is to X at [1/2,0,0], but we use [0,1/2,0] as representative since it is used elsewhere also
            LabeledKVec{3}(:ZA, KVec{3}("1/2,u,0"))),
        ["X₁⁺", "X₂⁺", "X₃⁺", "X₄⁺", "X₁⁻", "X₂⁻", "X₃⁻", "X₄⁻"],
        ["ZA₁", "ZA₂", "ZA₃", "ZA₄"],
        [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],
        false #=monodromy=#
    )
    subt_M_ZA = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:M, KVec{3}("1/2,1/2,0")), 
            LabeledKVec{3}(:ZA, KVec{3}("1/2,u,0"))),
        ["M₁⁺", "M₂⁺", "M₃⁺", "M₄⁺", "M₁⁻", "M₂⁻", "M₃⁻", "M₄⁻"],
        ["ZA₁", "ZA₂", "ZA₃", "ZA₄"],
        [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0],
        false #=monodromy=#
    )
    # with TR
    _throw_if_already_in_table(SUBDUCTIONSD_TR_3D[sgnum], [subt_X_ZA, subt_M_ZA])
    push!(SUBDUCTIONSD_TR_3D[sgnum], subt_X_ZA, subt_M_ZA)
    # without TR; same as above
    _throw_if_already_in_table(SUBDUCTIONSD_3D[sgnum], [subt_X_ZA, subt_M_ZA])
    push!(SUBDUCTIONSD_3D[sgnum], deepcopy(subt_X_ZA), deepcopy(subt_M_ZA))
    return nothing
end

# ---------------------------------------------------------------------------------------- #
# space group 110: missing B-connections from N to M & Γ 
#                  there is also a missing monodromy in N via B, but only for spinful case
#                  (we omit for now, since we only treat spinless cases)
function subduction_table_corrections_sg110_tr!(SUBDUCTIONSD_TR_3D) # with TR
    sgnum = 110
    subt_N_B = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:N, KVec{3}("1/2,0,1/2")), LabeledKVec{3}(:B, KVec{3}("u,0,w"))),
        ["N₁N₂"], ["B₁", "B₂"], [1 1;],
        false #=monodromy=#
    )
    subt_M_B = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:M, KVec{3}("1,1,1")), LabeledKVec{3}(:B, KVec{3}("u,0,w"))),
        ["M₁M₂", "M₃M₄", "M₅"], ["B₁", "B₂"], [0 2; 2 0; 1 1],
        false #=monodromy=#
    )
    subt_Γ_B = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:Γ, KVec{3}("0,0,0")), LabeledKVec{3}(:B, KVec{3}("u,0,w"))),
        ["Γ₁", "Γ₂", "Γ₃", "Γ₄", "Γ₅"], ["B₁", "B₂"], [1 0; 1 0; 0 1; 0 1; 1 1],
        false #=monodromy=#
    )
    subt_N′_B = SubductionTable{3}( # not technically necessary; only for spinful case & 
        sgnum,                      # for TR-broken, but include for consistency w/ no-TR
        Connection{3}(
            LabeledKVec{3}(:N′, KVec{3}("3/2,0,3/2")), LabeledKVec{3}(:B, KVec{3}("u,0,w"))
        ),
        ["N′₁N′₂"], ["B₁", "B₂"], [1 1;],
        true #=monodromy=#
    )
    # with TR
    _throw_if_already_in_table(SUBDUCTIONSD_TR_3D[sgnum], [subt_N_B, subt_M_B, subt_Γ_B, subt_N′_B])
    push!(SUBDUCTIONSD_TR_3D[sgnum], subt_N_B, subt_M_B, subt_Γ_B, subt_N′_B)
end
function subduction_table_corrections_sg110_notr!(SUBDUCTIONSD_3D) # without TR
    sgnum = 110
    # N-B and Γ-B are already in BCS tables, as is N-monodromy; but M-B is missing
    subt_M_B = SubductionTable{3}(
        sgnum,
        Connection{3}(
            LabeledKVec{3}(:M, KVec{3}("1,1,1")), LabeledKVec{3}(:B, KVec{3}("u,0,w"))),
        ["M₁", "M₂", "M₃", "M₄", "M₅"], ["B₁", "B₂"], [0 1; 0 1; 1 0; 1 0; 1 1],
        false #=monodromy=#
    )
    _throw_if_already_in_table(SUBDUCTIONSD_3D[sgnum], subt_M_B)
    push!(SUBDUCTIONSD_3D[sgnum], subt_M_B)
end

# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# some checking: make sure we don't add something already included in tables (in case they
# get fixed in the future)
function _has_connection(subts::AbstractVector{SubductionTable{D}}, c::Connection{D}) where D
    kᴳlab, kᴴlab = c.kᴳ.label, c.kᴴ.label
    for subtᵢ in subts
        cᵢ = subtᵢ.c
        kᴳlabᵢ, kᴴlabᵢ = cᵢ.kᴳ.label, cᵢ.kᴴ.label
        if kᴳlabᵢ == kᴳlab && kᴴlabᵢ == kᴴlab
            return true
        end
    end
    return false
end
function _has_connection(
        subts::AbstractVector{SubductionTable{D}},
        subt::SubductionTable{D}) where D
    _has_connection(subts, subt.c) 
end
function _throw_if_already_in_table(
    subts::AbstractVector{SubductionTable{D}},
    subt_to_add::SubductionTable{D}) where D
    if _has_connection(subts, subt_to_add)
        error(lazy"tabulated data already includes connection $(subt_to_add.c)")
    end
end
function _throw_if_already_in_table(
    subts::AbstractVector{SubductionTable{D}},
    subts_to_add::AbstractVector{SubductionTable{D}}) where D
    foreach(subt->_throw_if_already_in_table(subts, subt), subts_to_add)
end