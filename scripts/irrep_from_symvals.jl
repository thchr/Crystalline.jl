using SGOps, LinearAlgebra, Test, JuMP, GLPK
import SGOps: DEFAULT_ATOL, rotation, rotation_order_3d, irdim

if !isdefined(Main, :SGOpsHilbertBases)
    include((@__DIR__)*"/hilbert_basis.jl")
end
using Main.SGOpsHilbertBases

# NOTES TO SELF: SG23 seems to present a challenge: we can find a solution with ν=3, but Lu
# & Watanabe found that there must be at least 3 tranverse modes. It may be that we need to
# impose some other restrictions, that guarantee the prohibition of crossings between the
# transverse and longitudinal modes at k-points away from Γ. I am not really sure how we
# would manage this?
if !isdefined(Main, :MAX_FILL)
    const MAX_FILL = 20
end
function find_bandrepresentation_lowest_bands(sgnum::Integer; timereversal::Bool=true)
    # Find the Hilbert basis that respects the compatibility relations
    nsᴴ, _, BRS = compatibility_bases(sgnum, spinful=false, timereversal=timereversal)
    Nᴴ = size(nsᴴ, 2)

    # Irreps at Γ and irrep-multiplicities of ω=0 band-triad
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    timereversal && (lgirs = realify(lgirs))
    ms = find_zero_freq_gamma_representation(lgirs)

    # Find the position of the Γ irreps in the BandRepSet
    irlabs_brs = BRS.irreplabs
    irlabs_lgirs = SGOps.formatirreplabel.(label.(lgirs))
    Γiridxs_in_brs = map(irlab->findfirst(==(irlab), irlabs_brs), irlabs_lgirs)

    # Try to find a nonnegative expansion in minimal "band-numbers" of basis of `nsᴴ` 
    # subject to constraints from `ms`; define as a constrained feasibility problem
    m = Model(GLPK.Optimizer)
    @variable(m, c[1:Nᴴ] >= 0, Int)
    @constraint(m, nsᴴ[Γiridxs_in_brs,:]*c .>= ms)

    @variable(m, filling_constraint)                    # define `filling_constraint` as a variable,
    @constraint(m, nsᴴ[end,:]'*c == filling_constraint) # so its value can be changed on the fly
    for filling = 3:MAX_FILL
        # Impose constraints on the number of included bands (band filling)
        fix(filling_constraint, filling) # change/set `filling_constraint` value
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found
            println("Found feasible solution with ν = ", filling, " bands")
            return m, filling
        end

        # I think this might be doable instead (or optimizable) by just looking at the 
        # individual bases, and excluding those that do not have any elements at _some_ of 
        # the required Γ irreps. Might still need an optimization step though. Wonder if we
        # can try to find all the solutions by excluding a solution once we've found it?
    end
    return m, -1
end



function get_lgirreps_at_Γ(sgnum::Integer, Dᵛ::Val=Val(3)) # small irreps at Γ
   lgirs = first(get_lgirreps(sgnum,  Dᵛ))
   kv = kvec(first(lgirs))
   @assert all(iszero, kv.k₀) && isspecial(kv) # Make sure that lgirs indeed is sampled at Γ

   return lgirs
end

function find_zero_freq_gamma_representation(lgirs::AbstractVector{LGIrrep{3}})
    lg = group(first(lgirs))

    symvals = zero_freq_gamma_symvals(lg)

    return find_representation(symvals, lgirs)
end
function find_zero_freq_gamma_representation(sgnum::Integer; timereversal::Bool=false)   
    lgirs = get_lgirreps_at_Γ(sgnum, Val(3))
    timereversal && (lgirs = realify(lgirs))

    return find_zero_freq_gamma_representation(lgirs)
end

# symmetry eigenvalues of the two transverse branches
function zero_freq_gamma_transverse_symvals(ops::AbstractVector{SymOperation{3}})
    symvalsᵀ = [Vector{ComplexF64}(undef, 2) for _ in eachindex(ops)]

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = rotation_order_3d(det(W), tr(W))
        
        if rotval != -1                                                   # not inversion
            n = abs(rotval) # rotation order 
            hasmirror = signbit(rotval)

            if !hasmirror   # ← Ordinary rotation
                # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
                # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
                θ = 2π/n
                symvalsᵀ[i] .= (cis(θ), cis(-θ))
            
            else            # ← Roto-inversion or mirror
                # I don't believe this is generally possible to do for
                # the choice of transverse plane waves. E.g., by direct example: there are
                # no lines of symmetry from Γ that contain -4₀₀₁ in sg 81; nor any with 
                # -3₀₀₁ in sg 147; nor any with -6₀₀₁ in sg 174. As such, there is no line 
                # of symmetry from Γ along which transverse plane waves could be 
                # symmetry-allowed eigenfunctions for the rotoinversions

                # It _is_ possible for a simple mirror though (i.e., rotation followed by
                # inversion, i.e. -2 === m): the right choice is to pick the symmetry
                # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
                # each high-symmetry k-vector)
                if rotval == -2
                    symvalsᵀ[i] .= (+one(ComplexF64), -one(ComplexF64))

                elseif rotval ∈ (-3, -4, -6)
                    # Use NaN-sentinel value instead of "missing" for rotoinversions (arises
                    # for 82 SGs, namely for 81:88, 111:142, 147:148, 162:167, 174:176,
                    # 187:194, 200:206, and 215:230)
                    symvalsᵀ[i] .= ComplexF64(NaN, NaN)

                else
                    throw("Unexpected rotation value $rotval")
                end
            end

        else             # ← Inversion
            # This is a bit of a special case because, like rotoinversions, there are no
            # lines of symmetry from Γ that are invariant under inversion. However, if we 
            # recall that, overall, we are considering two transverse and a single
            # longitudinal plane wave branch (≡T²+L), which effectively transform like the
            # three Cartesian vectors at Γ, there should be no difficulty since, regardless
            # of partitioning into T² and L, each symmetry eigenvalue must be -1
            symvalsᵀ[i] .= (-one(ComplexF64), -one(ComplexF64))
            # TODO: I'm honestly not too sure about this. This is rather hand-wavy. I think
            #       it really may be analogous to the rotoinversions...
        end
    end

    return symvalsᵀ
end
function zero_freq_gamma_transverse_symvals(sg::Union{LittleGroup{3}, SpaceGroup{3}})
    ops = operations(primitivize(sg))
    return zero_freq_gamma_transverse_symvals(ops)
end


function zero_freq_gamma_symvals(ops::AbstractVector{SymOperation{3}})
    symvals = [Vector{ComplexF64}(undef, 3) for _ in eachindex(ops)]

    for (i, op) in enumerate(ops)
        W = rotation(op)
        rotval = rotation_order_3d(det(W), tr(W))
        n = abs(rotval)

        # This covers every case, including rotations, mirrors, rotoinversions, & inversion
        θ = 2π/n
        symvals[i] .= sign(rotval).*(cis(θ), cis(-θ), one(ComplexF64))
    end
    
    return symvals
end
function zero_freq_gamma_symvals(sg::Union{LittleGroup{3}, SpaceGroup{3}})
    ops = operations(primitivize(sg))
    return zero_freq_gamma_symvals(ops)
end


# ------------------------------------------------------------------------------------------
find_bandrepresentation_lowest_bands.(1:10);