# NOTE: There's some unfortunate conventional choices regarding the setting of rod groups:
#       the periodicity is assumed to be along the z-direction - this is contrary to how the
#       layer groups (xy periodicity) and frieze groups (x periodicity) are built... One
#       possible motivation for this choice would be to have rod groups more closely
#       resemble the corresponding space groups; but this way, the frieze groups don't
#       resemble the rod groups (instead, the layer groups)... Quite annoying in a
#       computational setting ...: we now need to keep track of _which_ direction is the
#       periodic one... Seems more appealing to just establish some convention (say,
#       P = 1 ⇒ x periodicity; P = 2 ⇒ xy periodicity). Unfortunately, if we do that,
#       we need to change the labels for the rod groups (they are setting dependent).
#       This problem is discussed in e.g. ITE1 Section 1.2.6; one option is to indicate the
#       direction of periodicity with a subscript (e.g., ₐ for x-direction; but that doesn't
#       work too well with unicode that doesn't have {b,c}-subscripts.
"""
    subperiodicgroup(num::Integer, ::Val{D}=Val(3), ::Val{P}=Val(2))
    subperiodicgroup(num::Integer, D::Integer, P::Integer)
                                                            --> ::SubperiodicGroup{D,P}

Return the operations of the subperiodic group `num` of embedding dimension `D` and
periodicity dimension `P` as a `SubperiodicGroup{D,P}`.

The setting choices are those of the International Tables for Crystallography, Volume E.

Allowed combinations of `D` and `P` and their associated group names are:

- `D = 3`, `P = 2`: Layer groups (`num` = 1, …, 80).
- `D = 3`, `P = 1`: Rod groups (`num` = 1, …, 75).
- `D = 2`, `P = 1`: Frieze groups (`num` = 1, …, 7).

## Example

```jldoctest
julia> subperiodicgroup(7, Val(2), Val(1))
SubperiodicGroup{2, 1} ⋕7 (𝓅2mg) with 4 operations:
 1
 2
 {m₁₀|½,0}
 {m₀₁|½,0}
```

## Data sources

The symmetry operations returned by this function were originally retrieved from the [Bilbao
Crystallographic Database, SUBPERIODIC GENPOS](https://www.cryst.ehu.es/subperiodic/get_sub_gen.html).
"""
function subperiodicgroup(num::Integer, Dᵛ::Val{D}=Val(3), Pᵛ::Val{P}=Val(2)) where {D, P}
    @boundscheck _check_valid_subperiodic_num_and_dim(num, D, P)
    codes = SUBG_CODES_Vs[(D,P)][num]

    cntr = centering(num, D, P)
    Ncntr = centering_volume_fraction(cntr, Dᵛ, Pᵛ)
    Nop = (length(codes)+1) # number of operations ÷ by equiv. centering translations
    operations = Vector{SymOperation{D}}(undef, Nop * Ncntr)
    
    # convert `codes` to `SymOperation`s and add to `operations`
    _include_symops_from_codes!(operations, codes)

    # add the centering-translation related operations (not included in `codes`)
    if Ncntr > 1
        cntr_translations = all_centeringtranslations(cntr, Dᵛ, Pᵛ)
        _include_symops_centering_related!(operations, cntr_translations, Nop)
    end
    
    return SubperiodicGroup{D, P}(num, operations)
end
subperiodicgroup(num::Integer, D::Integer, P::Integer) = subperiodicgroup(num, Val(D), Val(P))