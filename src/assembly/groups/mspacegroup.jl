function mspacegroup(BNS₁::Integer, BNS₂::Integer, Dᵛ::Val{D}=Val(3)) where D
    D == 3 || throw(DomainError("only 3D magnetic space groups are supported"))
    msgnum = (BNS₁, BNS₂)
    codes = MSG_CODES_D[msgnum]

    label = MSG_BNS_LABELs_D[msgnum]
    cntr = first(label)
    Ncntr = centering_volume_fraction(cntr, Dᵛ)
    Nop = (length(codes)+1) # number of operations ÷ by equiv. centering translations
    operations = Vector{MSymOperation{D}}(undef, Nop * Ncntr)

    operations[1] = MSymOperation{D}(one(SymOperation{D}), false)
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(ROTATIONS_3D[code[1]], TRANSLATIONS_3D[code[2]])
        operations[n+1] = MSymOperation{D}(op, code[3])
    end

    # add the centering-translation related operations (not included in `codes`)
    if Ncntr > 1
        cntr_translations = all_centeringtranslations(cntr, Dᵛ)
        for (i, t) in enumerate(cntr_translations)
            for n in 1:Nop
                mop = operations[n]
                t′ = reduce_translation_to_unitrange(translation(mop) + t)
                op′ = SymOperation{D}(mop.op.rotation, t′)
                mop′ = MSymOperation{D}(op′, mop.tr)
                operations[n+i*Nop] = mop′
            end
        end
    end

    return MSpaceGroup{D}(msgnum, operations)
end

function mspacegroup(OG₃::Integer, ::Val{D}=Val(3)) where D
    D == 3 || throw(DomainError("only 3-dimensional magnetic space groups are supported"))
    OG₃ ∈ 1:MAX_MSGNUM[D] || _throw_invalid_msgnum(OG₃, D)

    BNS₁, BNS₂ = MSG_OG₃2BNS_NUMs_V[OG₃]
    return mspacegroup(BNS₁, BNS₂, Val(3))
end
@noinline function _throw_invalid_msgnum(OG₃, D)
    throw(DomainError(OG₃, "the sequential magnetic space group number must be between 1 and $(MAX_MSGNUM[D]) in dimension $D"))
end