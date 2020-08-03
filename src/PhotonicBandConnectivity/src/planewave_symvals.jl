# (sum of) symmetry eigenvalues of ω=0 branches

# two transverse and two longitudinal plane waves (2T+1L)
function get_symvals²ᵀ⁺¹ᴸ(ops::AbstractVector{SymOperation{3}})
    symvals = Vector{ComplexF64}(undef, length(ops))
    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        n = abs(rotval)
        # This covers every case, including rotations, mirrors, rotoinversions, & inversion
        θ = 2π/n
        symvals[i] = sign(rotval)* (cis(θ) + cis(-θ) + one(ComplexF64))
    end
    
    return symvals
end

# single longitudinal plane wave (1L)
function get_symvals¹ᴸ(ops::AbstractVector{SymOperation{3}})
    symvals = ones(ComplexF64, length(ops)) 
end

# two transverse plane waves (2T)
function get_symvals²ᵀ(ops::AbstractVector{SymOperation{3}})
    symvals = Vector{ComplexF64}(undef, length(ops))

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        
        n = abs(rotval) # rotation order 

        if !signbit(rotval)                                     # ← Proper rotation
            # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
            # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
            θ = 2π/n
            symvals[i] = 2cos(θ) # eⁱᶿ + e⁻ⁱᶿ = 2cos(θ)

        else                                                    # ← Improper rotation
            # It is not generally possible to infer the all the symmetry eigenvalues of 
            # roto-inversions with rotval = (-1, -3, -4, -6) for the two transverse 
            # plane waves (2T) in isolation. This is because there are there no lines of
            # symmetry from Γ along which 2T could be symmetry-allowed eigenfunctions 
            # for the rotoinversions.
            # Instead, we pick a _possible_ choice for 1L and infer _possible_ symmetry
            # values from [2T+1L] - 1L

            # It _is_ possible for a simple mirror though (i.e., rotation followed by
            # inversion, i.e. -2 === m): the right choice is to pick the symmetry
            # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
            # each high-symmetry k-vector)
            if rotval == -2                     # ← Mirror
                symvals[i] = zero(ComplexF64) # (+1) + (-1)

            elseif rotval ∈ (-1, -3, -4, -6)    # ← Roto-inversions & inversion
                θ = 2π/n
                # In general, we have: 
                #   [2T+1L] - 1L = [(-eⁱᶿ) + (-e⁻ⁱᶿ) + (-1)] - (+1) = -2cos(θ) - 2.0
                # For inversion specifically, this is:
                #   [2T+1L] - 1L = [(-1) + (-1) + (-1)] - (+1) = -4
                symvals[i] = -2cos(θ) - 2.0 
                # SGs w/ non-mirror rotoinversions are 81:88, 111:142, 147:148, 162:167, 
                # 174:176, 187:194, 200:206, and 215:230
            end
        end
    end

    return symvals
end

# convenience accessors via space/little groups, ensuring primitive basis
for f in (:get_symvals²ᵀ⁺¹ᴸ, :get_symvals¹ᴸ, :get_symvals²ᵀ)
    @eval $f(sg::Union{LittleGroup{3}, SpaceGroup{3}}) = $f(operations(primitivize(sg)))
end