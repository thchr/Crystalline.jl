"""
    generators(num::Integer, ::Type{SubperiodicGroup{D,P}})  -->  ::Vector{SymOperation{D}}

Return a canonical set of generators for the subperiodic group `num` of embedding dimension
`D` and periodicity dimension `P`. See also [`subperiodicgroup`](@ref).

See also [`generators(::Integer, ::Type{SpaceGroup})`](@ref) and information therein.

## Example

```jldoctest
julia> generators(7, SubperiodicGroup{2, 1})
2-element Vector{SymOperation{2}}:
 2
 {m₁₀|½,0}
```

## Data sources

The generators returned by this function were originally retrieved from the [Bilbao
Crystallographic Database, SUBPERIODIC GENPOS](https://www.cryst.ehu.es/subperiodic/get_sub_gen.html).
"""
function generators(num::Integer, ::Type{SubperiodicGroup{D,P}}) where {D,P}
    @boundscheck _check_valid_subperiodic_num_and_dim(num, D, P)
    codes = SUBG_GENS_CODES_Vs[(D,P)][num]

    # convert `codes` to `SymOperation`s and add to `operations`
    operations = Vector{SymOperation{D}}(undef, length(codes))
    _include_symops_from_codes!(operations, codes; add_identity=false)

    return operations
end