

"""
    spacegroup(sgnum::Integer, ::Val{D}=Val(3))
    spacegroup(sgnum::Integer, D::Integer)          --> SpaceGroup{D}

Return the space group symmetry operations for a given space group number `sgnum` and 
dimensionality `D` as a `SpaceGroup{D}`.
The returned symmetry operations are specified relative to the conventional basis vectors,
i.e. are not necessarily primitive (see [`centering`](@ref)).
If desired, operations for the primitive unit cell can subsequently be generated using 
[`primitivize`](@ref) or [`Crystalline.reduce_ops`](@ref).

The default choices for the conventional basis vectors follow the conventions of the Bilbao
Crystallographic Server (or, equivalently, the International Tables of Crystallography), 
which are:

- Unique axis *b* (cell choice 1) for monoclinic space groups.
- Obverse triple hexagonal unit cell for rhombohedral space groups.
- Origin choice 2: inversion centers are placed at (0,0,0). (relevant for certain
  centrosymmetric space groups with two possible choices; e.g., in the orthorhombic,
  tetragonal or cubic crystal systems).

See also [`directbasis`](@ref).

## Data sources

The symmetry operations returned by this function were originally retrieved from the [Bilbao
Crystallographic Server, SPACEGROUP GENPOS](https://www.cryst.ehu.es/cryst/get_gen.html).
The associated citation is: ([Aroyo et al., Z. Kristallogr. Cryst. Mater. **221**, 15
(2006).](https://doi.org/10.1524/zkri.2006.221.1.15)).
"""
function spacegroup(sgnum, Dᵛ::Val{D}=Val(3)) where D
    @boundscheck _check_valid_sgnum_and_dim(sgnum, D)

    codes = D == 3 ? SG_CODES_3D_V[sgnum] : 
            D == 2 ? SG_CODES_2D_V[sgnum] :
            D == 1 ? SG_CODES_1D_V[sgnum] :
                     error("unreachable; D is invalid but boundscheck assumed valid")

    cntr = centering(sgnum, D)
    Ncntr = centering_volume_fraction(cntr, Dᵛ)
    Nop = length(codes)+1 # number of operations ÷ by equiv. centering translations
    operations = Vector{SymOperation{D}}(undef, Nop * Ncntr)

    # convert `codes` to `SymOperation`s and add to `operations`
    _include_symops_from_codes!(operations, codes)

    # add the centering-translation related operations (not included in `codes`)
    if Ncntr > 1
        cntr_translations = all_centeringtranslations(cntr, Dᵛ)
        _include_symops_centering_related!(operations, cntr_translations, Nop)
    end

    return SpaceGroup{D}(sgnum, operations)
end
spacegroup(sgnum::Integer, D::Integer) = spacegroup(sgnum, Val(D))

function _include_symops_from_codes!(operations::Vector{SymOperation{D}}, codes) where D
    operations[1] = one(SymOperation{D}) # add trivial identity op separately and manually
    for (n, code) in enumerate(codes)
        op = SymOperation{D}(get_indexed_rotation(code[1], Val{D}()), 
                             get_indexed_translation(code[2], Val{D}()))
        operations[n+1] = op
    end
    return operations
end

function _include_symops_centering_related!(
            operations::Vector{SymOperation{D}}, cntr_translations, Nop) where D
    for (i, t) in enumerate(cntr_translations)
        for n in 1:Nop
            op = operations[n]
            t′ = reduce_translation_to_unitrange(translation(op) + t)
            op′ = SymOperation{D}(op.rotation, t′)
            operations[n+i*Nop] = op′
        end
    end
    return operations
end

function _check_valid_sgnum_and_dim(sgnum::Integer, D::Integer)
    if D == 3 
        sgnum > 230 && _throw_invalid_sgnum(sgnum, D)
    elseif D == 2
        sgnum > 17  && _throw_invalid_sgnum(sgnum, D)
    elseif D == 1
        sgnum > 2   && _throw_invalid_sgnum(sgnum, D)
    else
        _throw_invalid_dim(D)
    end
    sgnum < 1 && throw(DomainError(sgnum, "group number must be a positive integer"))
    return nothing
end