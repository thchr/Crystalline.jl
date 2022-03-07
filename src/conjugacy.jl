"""
    classes(ops::AbstractVector{SymOperation{D}}, [cntr::Union{Char, Nothing}])
                                                    -->  Vector{Vector{SymOperation{D}}}

Return the conjugacy classes of a group ``G`` defined by symmetry operations `ops`.

## Definitions
Two elements ``a`` and ``b`` in ``G`` are considered conjugate if there exists a ``g ∈ G``
such that ``gag^{-1} = b``. This defines an equivalence relation ``\\sim``, i.e., we say
that ``a \\sim b`` if ``a`` and ``b`` are conjugate.
The conjugacy classes of ``G`` are the distinct equivalence classes that can be identified
under this equivalence relation, i.e. the grouping of ``G`` into subsets that are equivalent
under conjugacy.

## Extended help
If `ops` describe operations in a crystal system that is not primitive (i.e., if its
[`centering`](@ref) type is not `p` or `P`) but is presented in a conventional setting,
the centering symbol `cntr` _must_ be given. If `ops` is not in a centered crystal
system, or if `ops` is already reduced to a primitive setting, `cntr` should be given as
`nothing` (default behavior) or, alternatively, as `P` or `p` (depending on dimensionality).

A single-argument calls to `classes` with `SpaceGroup` or `LittleGroup` types will
assume that `ops` is provided in a conventional setting, i.e., will forward the method call
to `classes(ops, centering(ops, dim(ops)))`. To avoid this behavior (if `ops` was already
reduced to a primitive setting prior to calling `classes`), `cntr` should be provided
explicitly as `nothing`.
"""
function classes(
            ops::AbstractVector{SymOperation{D}},
            cntr::Union{Char, Nothing}=nothing
            ) where D

    ops⁻¹ = inv.(ops)
    cntr_ops = if cntr === nothing
        nothing
    else
        if cntr == 'P' || cntr == 'p'
            SymOperation{D}[]
        else
            SymOperation{D}.(all_centeringtranslations(cntr, Val(D)))
        end
    end

    conj_classes = Vector{Vector{SymOperation{D}}}()
    classified = sizehint!(BitSet(), length(ops))
    i = 0
    while length(classified) < length(ops)
        i += 1
        i ∈ classified && continue
        a = ops[i]
        conj_class = push!(conj_classes, [a])[end]
        push!(classified, i)
        for (g, g⁻¹) in zip(ops, ops⁻¹)
            b = g*a*g⁻¹
            if all(Base.Fix1(!≈, b), conj_class) # check that `b` is not already in class
                add_to_class!(classified, conj_class, b, ops)
            end

            # special treatment for cases where we have a space or little group that
            # is not in a primitive setting (but a conventional setting) and might
            # additionally not include "copies" of trivial centering translations among
            # `ops`; to account for this case, we just create translated copies of `b`
            # over all possible centers and check those against `ops` as well until
            # convergence (TODO: there's probably a much cleaner/more performant way of
            # doing this)
            if cntr_ops !== nothing && length(cntr_ops) > 0
                all_additions = 0
                new_additions = 1 # initialize to enter loop
                while new_additions != 0 && all_additions < length(cntr_ops)
                    # this while loop is here to make sure we don't depend on ordering of
                    # translation additions to `conj_class`, i.e. to ensure convergence
                    new_additions = 0
                    for cntr_op in cntr_ops
                        b′ = cntr_op*b
                        if all(Base.Fix1(!≈, b′), conj_class)
                            new_additions += add_to_class!(classified, conj_class, b′, ops)
                        end
                    end
                    all_additions += new_additions
                end
            end
        end
    end
    return conj_classes
end
classes(g::AbstractGroup) = classes(g, centering(g))

# adds `b` to `class` and index of `b` in `ops` to `classified`
function add_to_class!(classified, class, b, ops)
    i′ = findfirst(op -> isapprox(op, b, nothing, false), ops)
    if i′ !== nothing
        # `b` might not exist in `ops` if `ops` refers to a group of reduced symmetry
        # operations (i.e. without centering translation copies) that nevertheless is
        # still in a non-primitive setting; in that case, we just don't include `b`
        push!(classified, i′)
        push!(class, ops[i′]) # better to take `ops[i′]` than `b` cf. possible float errors
        return true
    else
        return false # added nothing
    end
end


"""
    is_abelian(ops::AbstractVector{SymOperation}, [cntr::Union{Char, Nothing}])  -->  Bool

Return the whether the group composed of the elements `ops` is Abelian.

A group ``G`` is Abelian if all its elements commute mutually, i.e., if
``g = hgh^{-1}`` for all ``g,h ∈ G``.

See discussion of the setting argument `cntr` in [`classes`](@ref).
"""
function is_abelian(ops::AbstractVector{SymOperation{D}}, cntr::Union{Char, Nothing}=nothing) where D
    return length(classes(ops, cntr)) == length(ops)
end
is_abelian(g::AbstractGroup) = is_abelian(g, centering(g))