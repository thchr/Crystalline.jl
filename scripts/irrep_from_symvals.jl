using Pkg

# activate environment of file's folder if not already active
dirname(Pkg.project().path) == (@__DIR__) || Pkg.activate(@__DIR__)

using Crystalline
import Crystalline: DEFAULT_ATOL, rotation, irdim

if !isdefined(Main, :CrystallineHilbertBases)
    include((@__DIR__)*"/hilbert_basis.jl")
end
using Main.CrystallineHilbertBases

# NOTES TO SELF: SG23 seems to present a challenge: we can find a solution with ν=3, but Lu
# & Watanabe found that there must be at least 3 tranverse modes. It may be that we need to
# impose some other restrictions, that guarantee the prohibition of crossings between the
# transverse and longitudinal modes at k-points away from Γ. I am not really sure how we
# would manage this?
if !isdefined(Main, :MAX_FILL)
    const MAX_FILL = 20
end
# ------------------------------------------------------------------------------------------
function find_bandrepresentation_lowest_bands(
            sgnum::Integer; 
            timereversal::Bool=true, 
            verbose::Bool=true
            )
    # Irreps at Γ, irrep-multiplicities of ω=0 2T bands, and symmetry operations
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    timereversal && (lgirs = realify(lgirs))
    lg = group(first(lgirs))
    rotvals = map(op->(W=rotation(op); Crystalline.rotation_order_3d(det(W), tr(W))), lg)

    # 2T irreps; check if "simple treatment"/fast-path is applicable
    ms²ᵀ = find_zero_freq_gamma_transverse_representation(lgirs)
    has_nonmirror_improper = any(∈((-1, -3, -4, -6)), rotvals)
    is_regular²ᵀ = all(≥(0), ms²ᵀ)
    if !has_nonmirror_improper && is_regular²ᵀ # All symvals known & regular 2T irrep

        # Simple case: if there are no non-mirror improper rotations, we can directly infer 
        # the irrep of the 2T branches. If that irrep is regular (i.e. has no negative
        # coefficients), we don't need to invoke 1L at all, and can solve for just 2T alone.
        ms = ms²ᵀ
#
#        return nothing
#    elseif !has_nonmirror_improper # All symvals known & irregular 2T and regular 1L irreps
#        ms¹ᴸ = find_zero_freq_gamma_longitudinal_representation(lgirs)
#        ms = ms²ᵀ + ms¹ᴸ
#        @assert all(≥(0), ms¹ᴸ)                              # check: 1L irrep regular (Γ₁)
#        @assert ms == find_zero_freq_gamma_representation(lgirs) # →: [2T+1L] = 2T+1L
#
#        nsᴴ, Γidxs = _compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
#        allowable = allowable_nsᴴ_idxs(ms¹ᴸ, nsᴴ)
#        return nothing
    else # Not all symvals known; multiple irrep options
        # TODO: implement non-fast path treatment
        # Main idea is to explore every choice for 1L and 2T that is _feasible_ given the 
        # existing symmetry constraints.
        return nothing
    end
    nsᴴ, Γidxs = _compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    Nᴴ = size(nsᴴ, 2)
    # Try to find a nonnegative expansion in minimal "band-numbers" of basis of `nsᴴ` 
    # subject to constraints from `ms`; define as a constrained feasibility problem
    m = Model(GLPK.Optimizer)
    @variable(m, c[1:Nᴴ] >= 0, Int)
    @constraint(m, nsᴴ[Γidxs,:]*c .>= ms)

    # Define `filling_constraint` as a variable, so its value can be changed on the fly
    @variable(m, filling_constraint)
    @constraint(m, nsᴴ[end,:]'*c == filling_constraint)
    for filling = 2:MAX_FILL
        # Impose constraints on the number of included bands (band filling)
        fix(filling_constraint, filling) # Change/set `filling_constraint` value
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found
            verbose && println("   Found feasible solution with ν = ", filling, " bands")
            return m, filling
        end

        # I think this might be doable instead (or optimizable) by just looking at the 
        # individual bases, and excluding those that do not have any elements at _some_ of 
        # the required Γ irreps. Might still need an optimization step though. Wonder if we
        # can try to find all the solutions by excluding a solution once we've found it?
    end
    return m, -1
end

function _compatibility_bases_and_Γidxs(sgnum, lgirs, timereversal)
    # Find the Hilbert basis that respects the compatibility relations
    nsᴴ, _, BRS = compatibility_bases(sgnum, spinful=false, timereversal=timereversal)

    # Find the indices of the Γ irreps in `BRS::BandRepSet` (and hence in `nsᴴ`), and how  
    # they map to the corresponding irrep indices in `lgirs`
    irlabs_brs = BRS.irreplabs
    irlabs_lgirs = Crystalline.formatirreplabel.(label.(lgirs))
    Γidxs = map(irlab->findfirst(==(irlab), irlabs_brs), irlabs_lgirs)

    return nsᴴ, Γidxs
end


# ------------------------------------------------------------------------------------------
function get_lgirreps_at_Γ(sgnum::Integer, Dᵛ::Val=Val(3)) # small irreps at Γ
   lgirs = first(get_lgirreps(sgnum,  Dᵛ))
   kv = kvec(first(lgirs))
   @assert all(iszero, kv.k₀) && isspecial(kv) # Make sure that lgirs indeed is sampled at Γ

   return lgirs
end

function find_zero_freq_gamma_transverse_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_transverse_symvals(lg)

    return find_representation(symvals, lgirs, nothing)
end
function find_zero_freq_gamma_longitudinal_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_longitudinal_symvals(lg)

    return find_representation(symvals, lgirs, nothing)
end
function find_zero_freq_gamma_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))
    symvals = zero_freq_gamma_symvals(lg)

    return find_representation(symvals, lgirs)
end
# convenience accessors
for f in (:find_zero_freq_gamma_transverse_representation, 
          :find_zero_freq_gamma_longitudinal_representation,
          :find_zero_freq_gamma_representation)
    @eval function $f(sgnum::Integer; timereversal::Bool=false)
        lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
        timereversal && (lgirs = realify(lgirs))

        return $f(lgirs)
    end
end

# ------------------------------------------------------------------------------------------
# (sum of) symmetry eigenvalues of ω=0 branches
function zero_freq_gamma_transverse_symvals(ops::AbstractVector{SymOperation{3}})
    symvalsᵀ = Vector{ComplexF64}(undef, length(ops))

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)
        
        if rotval != -1     # not inversion
            n = abs(rotval) # rotation order 
            isimproper = signbit(rotval)

            if !isimproper  # ← "Proper" rotation
                # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
                # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
                θ = 2π/n
                symvalsᵀ[i] = 2cos(θ) # eⁱᶿ + e⁻ⁱᶿ = 2cos(θ)
            else            # ← Roto-inversion or mirror
                # This is not generally possible to infer for transverse plane waves alone. 
                # E.g., by direct example: there are
                # no lines of symmetry from Γ that contain -4₀₀₁ in sg 81; nor any with 
                # -3₀₀₁ in sg 147; nor any with -6₀₀₁ in sg 174. As such, there is no line 
                # of symmetry from Γ along which transverse plane waves could be 
                # symmetry-allowed eigenfunctions for the rotoinversions

                # It _is_ possible for a simple mirror though (i.e., rotation followed by
                # inversion, i.e. -2 === m): the right choice is to pick the symmetry
                # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
                # each high-symmetry k-vector)
                if rotval == -2
                    symvalsᵀ[i] = zero(ComplexF64) # (+1) + (-1)
                elseif rotval ∈ (-3, -4, -6)
                    # SGs with non-mirror rotoinversions 81:88, 111:142, 147:148, 162:167, 
                    # 174:176, 187:194, 200:206, and 215:230
                    θ = 2π/n
                    symvalsᵀ[i] = -2cos(θ) - 2.0 # [2T+1L] - 1L = [(-eⁱᶿ) + (-e⁻ⁱᶿ) + (-1)] - (+1)
                else
                    throw("Unexpected rotation value $rotval")
                end
            end
        else                # ← Inversion
            # This is a bit of a special case because, like rotoinversions, there are no
            # lines of symmetry from Γ that are invariant under inversion. However, if we 
            # recall that, overall, we are considering two transverse and a single
            # longitudinal plane wave branch (≡T²+L), which effectively transform like the
            # three Cartesian vectors at Γ, there should be no difficulty since, regardless
            # of partitioning into T² and L, each symmetry eigenvalue must be -1
            symvalsᵀ[i] = -4.0 # [2T+1L] - 1L = [(-1) + (-1) + (-1)] - (+1)
        end
    end

    return symvalsᵀ
end

function zero_freq_gamma_longitudinal_symvals(ops::AbstractVector{SymOperation{3}})
    symvalsᴸ = ones(ComplexF64, length(ops)) 
    #=
    Vector{ComplexF64}(undef, length(ops))
    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = Crystalline.rotation_order_3d(W)  
        if rotval != -1     # not inversion
            n = abs(rotval) # rotation order 
            isimproper = signbit(rotval)
            if !isimproper  # ← Ordinary rotation
                symvalsᴸ[i] = one(ComplexF64)
            else            # ← Roto-inversion or mirror
                if rotval == -2
                    symvalsᴸ[i] = one(ComplexF64)
                elseif rotval ∈ (-3, -4, -6)
                    symvalsᴸ[i] = ComplexF64(NaN) # Indeterminate
                else
                    throw("Unexpected rotation value $rotval")
                end
            end
        else                # ← Inversion
            symvalsᴸ[i] = -one(ComplexF64)
        end
    end
    return symvalsᴸ
    =#
end

function zero_freq_gamma_symvals(ops::AbstractVector{SymOperation{3}})
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
# convenience accessors via space/little groups, ensuring primitive basis
for f in (:zero_freq_gamma_transverse_symvals, :zero_freq_gamma_symvals)
    @eval $f(sg::Union{LittleGroup{3}, SpaceGroup{3}}) = $f(operations(primitivize(sg)))
end
# ------------------------------------------------------------------------------------------
function zero_freq_resolve_representation(sgnum::Integer)
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    lg = group(first(lgirs))

    symvals = sum.(zero_freq_gamma_symvals(lg))
    symvalsᵀ = sum.(zero_freq_gamma_transverse_symvals(lg)) # Transverse (T²) branches
    symvalsᴸ = zero_freq_gamma_longitudinal_symvals(lg)     # Longitudinal (L) branch


    ct = CharacterTable(lgirs, nothing)
    χs = ct.chartable # character table as matrix (irreps-as-columns & operations-as-rows)
    ms = find_representation(symvals, lgirs)
    if !any(isnan, symvalsᵀ) # "Trivial" case
        # We know how all the symvals split into T² and L branches
        
        # find multiplicities of T² and L branches, allowing non-integer multiplicities
        msᵀ = find_representation(symvalsᵀ, lgirs, nothing, Float64)
        msᴸ = find_representation(symvalsᴸ, lgirs, nothing, Float64)

        hasinv = any(==(SymOperation("-x,-y,-z")), operations(lg))
        println(hasinv ? "Has inversion" : "")
        println.(msᵀ, Ref(" | "), msᴸ)
        println()
        # check consistency of splitting
        @assert symvals ≈ symvalsᵀ + symvalsᴸ
        @assert msᴸ + msᵀ ≈ ms
        return true
        
        # TODO: Continue from here - try to understand if you can find M in these simpler cases
    else                     # Nontrivial case
        
        # We do not know how all the symvals split into T² and L branches
        # Bail out for now
        return false
        #=
        known_idxs = findall(!isnan, symvalsᵀ)
        symvalsᵀ′ = @views symvalsᵀ[known_idxs]
        symvalsᴸ′ = @views symvalsᴸ[known_idxs]
        χs′ = @views χs[known_idxs, :]
        mpᵀ = χs′\symvalsᵀ′ 
        mpᴸ = χs′\symvalsᴸ′
        N = nullspace(χs′)
        =#
    end
        
end

# ------------------------------------------------------------------------------------------
if false
νs = Vector{Int}(undef, MAX_SGNUM[3])
for sgnum in 1:MAX_SGNUM[3]
    print("SG", sgnum, ": ")
    mν = find_bandrepresentation_lowest_bands(sgnum, timereversal=true, verbose=false)
    if mν != nothing
        m, ν = mν
    else
        ν = -1
    end
    νs[sgnum] = ν

    println(" "^(4-ndigits(sgnum)), ν, " bands")
end

# Compare with Watanabe & Lu
include("scripts/watanabelu_results.jl") # loads Watanabe & Lu data (in `pairings`)
Q = [[sg, M, Mbound] for (sg, M, Mbound) ∈ zip(1:230, νs, getindex.(pairings, 2))]
Q′ = filter(x-> x[2]≠-1, Q) # filter out those sgs that are not currently implemented (i.e. allow only regular 2T)
issues = map(x->x[2]==-1 ? "─" : (x[2]≥(x[3]) ? " " : "!"), Q)
differences = map(x->x[2]==-1 ? "─" : (x[2]==(x[3]) ? " " : "≠"), Q)
map(vcat.(Q, issues, differences)) do x
    println("|", " "^(4-ndigits(x[1])), x[1], " |", " "^(3-ndigits(x[2])),  # SG no.
            x[2] ==-1 ? "─" : x[2], " | ",      # our M predictions
            x[3] == 2 ? "=" : "≥", x[3], " | ", # M-bounds from Watanabe & Lu
            x[4], " | ",                        # issues
            x[5], " |"                          # differences from bound?
    )
end;


end