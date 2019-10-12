"""
    __write_all_unique_irrep_matrices(LGIRS) --> Nothing (writes to disk)

Write all the **unique** (i.e. different by more than a scalar)
irrep matrices in 3D to a .jl file, for subsequent loading as a constant
`Tuple` of `Tuple`s of `Matrix{ComplexF64}`s with the variable name
`UNIQUE_MATRICES_3D` (with the number of such unique matrices per
irrep dimension stored in the tuple `NUMBER_OF_UNIQUE_MATRICES_3D`).

(The highest irrep dimensionality in 3D is 6.)

This approach is probably pretty awful on the compiler, but it 
has a small storage footprint {~46 kbyte) and is fast at runtime.
"""
function __write_all_unique_irrep_matrices(LGIRS)
    D = dim(kvec(first(first(first(LGIRS)))))
    unique_ms = find_all_unique_irrep_matrices(LGIRS)
    N = length(unique_ms)
    n_unique = zeros(Int64,6)

    # sort by dimensionality and write to const variable in .jl file
    writename = (@__DIR__)*"/../data/lgirreps/"*string(D)*"d/unique_matrices"

    open(writename*".jl", "w") do file
        file_context = IOContext(file, :typeinfo => ComplexF64, :compact => false)
        # --- write matrices to Tuple of Tuples of Matrix{ComplexF64} ---
        write(file, "const UNIQUE_MATRICES_3D = (\n")
        # loop over irrep matrix dimensions (highest is 6 in 3D)
        for irdim = 1:6
            write(file, "\t# irdim = ", string(irdim))
            # special-case irdims 1 and 5
            if irdim == 1  
                write(file, "\n\t(ones(ComplexF64,1,1),),\n")
                n_unique[irdim] = 1
            elseif irdim == 5
                write(file, ": no matrices\n\t(Matrix{ComplexF64}(undef,0,0),),\n")
            else
                write(file, "\n\t(")
                for m in unique_ms
                    if size(m, 1) == irdim
                        n_unique[irdim] += 1
                        # when writing with :compact = false, complex numbers are written
                        # with a space between real and imag parts; jump through some hoops
                        # to remove those spaces
                        tempio = IOBuffer()
                        tempio′ = IOContext(tempio, :typeinfo => ComplexF64, :compact => false)
                        Base._show_nonempty(tempio′, m, "")
                        array_str = replace(String(take!(tempio)), r" (\+|\-) "=>s"\1") # replaces " + " by "+" and " - " by "-"

                        write(file, "\n\t\t", array_str, ",")
                    end
                end
            
            write(file, "\n\t),\n")
            end
        end
        write(file, ")")

        # write tuple indicating how many unique matrices per dimension
        write(file, "\n\n\nconst NUMBER_OF_UNIQUE_MATRICES_3D = ")
        write(file, repr(tuple(n_unique...)))
    end
    sum(n_unique) ≠ N && throw(error("Not all unique matrices were written to disk ($(sum(n_unique))/$(N))"))
    return nothing
end
__write_all_unique_irrep_matrices() = __write_all_unique_irrep_matrices(parselittlegroupirreps())

function _unique_matrices_up_to_scalar(ms)#::AbstractVector{Matrix{ComplexF64}})
    unique_ms = Vector{Matrix{ComplexF64}}()
    for m in ms
        mdim = size(m, 1)   # always square
        hypothesis = true   # hypothesis: m is a "new" matrix 
        for m′ in unique_ms # test hypothesis against current list of "unique" matrices
            if mdim == size(m′,1)
                S = m\m′
                # If m is proportional to one of the existing "unique" matrices m′, 
                # we will have S≡m⁻¹m′ equal to a equal-diagonal matrix, i.e. a 
                # UniformScaling. In that case our initial "new m" hypothesis is false.
                if isapprox(UniformScaling(S[1,1]), S, rtol=1e-11) # check if S is a UniformScaling
                    hypothesis = false
                    continue
                end
            end
        end
        if hypothesis # our hypothesis was true; add m to list of unique matrices
            push!(unique_ms, m)
        end
    end
    return unique_ms
end

@noinline function match_parallel_irrep_matrix(m)
    irdim = size(m, 1)
    if irdim == 1
        return UInt8(1), m[1]
    else
        for (idx, m′) in enumerate(UNIQUE_MATRICES_3D[irdim])
            S = m\m′
            prefac = mean(@view S[diagind(S)])
            # If m is proportional to one of the existing "unique" matrices m′, 
            # we will have S≡m⁻¹m′ equal to a equal-diagonal matrix, i.e. a UniformScaling.
            if isapprox(UniformScaling(prefac), S, rtol=1e-11) # check if S is a UniformScaling
                return convert(UInt8, idx), prefac
            end
        end
    end
end

function find_all_unique_irrep_matrices(LGIRS)
    unique_ms = Vector{Matrix{ComplexF64}}()
    for lgirsvec in LGIRS
        ms = (m for lgirs in lgirsvec for lgir in lgirs for m in lgir.matrices)
        append!(unique_ms, _unique_matrices_up_to_scalar(ms))
    end
    unique_ms = _unique_matrices_up_to_scalar(unique_ms)
    return unique_ms
end