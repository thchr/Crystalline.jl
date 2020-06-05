module SymmetryBases

using Crystalline, PyCall, SmithNormalForm, Test, JuMP, GLPK, PrettyTables
import Base: OneTo, show, size, getindex, firstindex, lastindex, IndexStyle, length
import Crystalline: matrix

export SymBasis, fillings, matrix
export compatibility_bases, nontopological_bases, split_fragiletrivial_bases

const PyNormaliz = pyimport("PyNormaliz") # import the PyNormaliz library

# -----------------------------------------------------------------------------------------

struct SymBasis <: AbstractVector{Vector{Int}}
    symvecs::Vector{Vector{Int}}
    irlabs::Vector{String}
    klabs::Vector{String}
    kvs::Vector{KVec}
    kv2ir_idxs::Vector{UnitRange{Int}} # pick k-point; find assoc. ir indices
    sgnum::Int
    spinful::Bool
    timeinvar::Bool
    allpaths::Bool
    compatbasis::Bool
end
function SymBasis(nsᴴ::AbstractMatrix{Int}, BRS::BandRepSet, compatbasis::Bool=true)
    kv2ir_idxs = [(f = irlab -> klabel(irlab)==klab; 
                   findfirst(f, BRS.irlabs):findlast(f, BRS.irlabs)) for klab in BRS.klabs]
    return SymBasis(collect(eachcol(nsᴴ)),
                    BRS.irlabs, BRS.klabs, BRS.kvs, kv2ir_idxs, 
                    BRS.sgnum, BRS.spinful, BRS.timeinvar, BRS.allpaths, compatbasis)
end

# accessors
matrix(sb::SymBasis) = hcat(sb.symvecs...)
vecs(sb::SymBasis)   = sb.symvecs
num(sb::SymBasis)    = sb.sgnum
irreplabels(sb::SymBasis) = sb.irlabs
klabels(sb::SymBasis)     = sb.klabs
isspinful(sb::SymBasis)   = sb.spinful
istimeinvar(sb::SymBasis) = sb.timeinvar
hasnonmax(sb::SymBasis)   = sb.allpaths
iscompatbasis(sb::SymBasis) = sb.compatbasis
fillings(sb::SymBasis)    = [nᴴ[end] for nᴴ in sb.symvecs]
length(sb::SymBasis)      = length(vecs(sb))

# define the AbstractArray interface for SymBasis
size(sb::SymBasis) = (length(sb),)
getindex(sb::SymBasis, keys...) = vecs(sb)[keys...]
firstindex(::SymBasis) = 1
lastindex(sb::SymBasis) = length(vecs(sb))
IndexStyle(::SymBasis) = IndexLinear()

# show method
function show(io::IO, ::MIME"text/plain", sb::SymBasis)
    Nⁱʳʳ = length(sb[1]) - 1

    # print a "title" line and the irrep labels
    println(io, iscompatbasis(sb) ? "Compatibility" : "Nontopological",
                " SymBasis (#", num(sb), "): ",
                length(sb), " Hilbert bases, sampling ",
                Nⁱʳʳ, " LGIrreps ",
                "(spin-", isspinful(sb) ? "½" : "1", " ",
                istimeinvar(sb) ? "w/" : "w/o", " TR)")

    k_idx = (i) -> findfirst(==(klabel(irreplabels(sb)[i])), klabels(sb)) # highlighters
    h_odd = Highlighter((data,i,j) -> i≤Nⁱʳʳ && isodd(k_idx(i)), crayon"light_blue")
    h_ν   = Highlighter((data,i,j) -> i==Nⁱʳʳ+1,                 crayon"light_yellow")

    pretty_table(io, 
        # table contents
        matrix(sb),
        # header
        eachindex(sb),
        # row names
        row_names = vcat(sb.irlabs, "ν"),
        # options/formatting/styling
        formatters = (v,i,j) -> iszero(v) ? "·" : string(v),
        vlines = [1,], hlines = [:begin, 1, Nⁱʳʳ+1, :end],
        row_name_alignment = :l,
        alignment = :c, 
        highlighters = (h_odd, h_ν), 
        header_crayon = crayon"bold"
        )

    # print k-vec labels
    print(io, "  KVecs (", hasnonmax(sb) ? "incl. non-maximal" : "maximal only", "): ")
    join(io, klabels(sb), ", ")
end

# -----------------------------------------------------------------------------------------

# All band structures can be written as 𝐧 = B𝐩 with pᵢ∈ℚ and nᵢ∈𝐍, and B a matrix whose 
# columns are EBRs. We can decompose B to a Smith normal form B = SΛT, such that all 
# allowable band structures can be written as 𝐧 = S𝐳. Here, S is an integer matrix with 
# elements Sᵢⱼ∈ℤ. To make nᵢ integer, we thus require zᵢ∈ℤ.

"""
    compatibility_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; algorithm)
    compatibility_bases(sgnum::Integer; kwargs...)

Computes the Hilbert bases associated with a Smith normal form `F` of the EBR matrix or from
a space group number `sgnum`, which respects all compatibility relations. The resulting 
bases form a non-negative integer coefficient basis for all possible band structures.

If the method is called with `sgnum::Integer`, the associated `BandRepSet` is also returned.

Several keyword arguments `kwargs` are possible:

    - `algorithm::String`: controls the algorithm used by Normaliz to compute the Hilbert
    basis. Choices are `"DualMode"` (default) and `"PrimalMode"`
    - `spinful::Bool`: Use single- (`false`, default) or double-valued (`true`) irreps.
    - `timereversal::Bool`: Assume presence (`true`, default) or absence (`false`) of
    time-reversal symmetry.
"""
function compatibility_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; 
                             algorithm::String="DualMode")
    # To restrict nᵢ to only positive integers, i.e. ℕ, the values of zᵢ must be such that 
    # ∑ⱼ Sᵢⱼzⱼ ≥ 0. This defines a set of inequalities, which in turn defines a polyhedral
    # integer cone. This is where (Py)Normaliz comes in.
    dᵇˢ = count(!iszero, F.SNF)           # "Dimensionality" of band structure
    S = @view F.S[:,OneTo(dᵇˢ)]           # All the nontrivial conditions on zⱼ

    C = PyNormaliz.Cone(inequalities = S) # Construct cone consistent with Sᵢⱼzⱼ ≥ 0
    C.Compute("HilbertBasis", algorithm)  # Compute¹ the Hilbert basis
    zsᴴ  = transpose(C.HilbertBasis())    # Columns are Hilbert basis vectors in 𝐳-space

    nsᴴ  = S*zsᴴ                          # Columns are Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ, BRS, true), zsᴴ  # Bases of all valid symmetry vectors in 𝐧- and 𝐲-space
end

"""
    nontopological_bases(F::SmithNormalForm.Smith, BRS::BandRepSet; algorithm)
    nontopological_bases(sgnum::Integer; kwargs...)

Computes the "non-topological" Hilbert bases associated with a Smith normal form `F` of the
EBR matrix or from a space group number `sgnum`, which forms a non-negative, integer span of
all non-topological band structures (i.e. both trivial and fragile-topological).

If the method is called with `sgnum::Integer`, the associated `BandRepSet` is also returned.

For possible keyword arguments `kwargs`, see `compatibility_bases(..)`.
"""
function nontopological_bases(F::SmithNormalForm.Smith, BRS::BandRepSet;
                              algorithm::String="DualMode")
    # To find _all_ nontopological bases we build a cone subject to the inequalities 
    # (SΛy)ᵢ ≥ 0 with yᵢ ∈ ℤ, which automatically excludes topological cases (since they
    # correspond to rational yᵢ)
    dᵇˢ = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    S = @view F.S[:,OneTo(dᵇˢ)]  # All the nontrivial conditions on zⱼ
    Λ = @view F.SNF[OneTo(dᵇˢ)]  # Nonzero invariant factors of Smith normal decomposition
    SΛ = S .* Λ' # Equivalent to S*diagm(dᵇˢ, dᵇˢ, Λ); all the nontrivial conditions on yᵢ
    
    C_nontopo = PyNormaliz.Cone(inequalities = SΛ)    # Cone consistent with SᵢⱼΛⱼyⱼ ≥ 0
    C_nontopo.Compute("HilbertBasis", algorithm)      # Compute¹ the Hilbert basis
    ysᴴ_nontopo = transpose(C_nontopo.HilbertBasis()) # Hilbert basis vectors in 𝐲-space

    nsᴴ_nontopo = SΛ*ysᴴ_nontopo                      # Hilbert basis vectors in 𝐧-space

    return SymBasis(nsᴴ_nontopo, BRS, false), ysᴴ_nontopo # Bases of nontopological states
end

# Convenience accessors from a space group number alone
for f in (:compatibility_bases, :nontopological_bases)
    @eval begin
        function $f(sgnum::Integer; algorithm::String="DualMode", spinful::Bool=false,
                                    timereversal::Bool=true)
            BRS = bandreps(sgnum, spinful=spinful, timereversal=timereversal)
            B   = matrix(BRS, true)    # Matrix with columns of EBRs.
            F   = Crystalline.smith(B) # Smith normal decomposition of B

            return $f(F, BRS, algorithm=algorithm)..., BRS
        end
    end
end

function split_fragiletrivial_bases(sb_nontopo::SymBasis, B::AbstractMatrix)
    if sb_nontopo.compatbasis
        throw(DomainError(sb_nontopo, "Specified SymBasis must have compatbasis=false"))
    end
    # Every vector of nsᴴ_nontopo that has a non-negative integer coefficient expansion in
    # EBRs, i.e. in B, represents a trivial basis element. All other elements represent 
    # fragile basis elements. We can just test whether such a solution exists through
    # constrained optimization, and then partition into trivial and fragile categories
    Nᴱᴮᴿ = size(B, 2)
    trivial_idxs = Int[]; fragile_idxs = Int[]
    for (j, nᴴ_nontopoʲ) in enumerate(sb_nontopo)
        m = Model(GLPK.Optimizer)
        @variable(m, c[1:Nᴱᴮᴿ] >= 0, Int)
        @constraint(m, B*c .== nᴴ_nontopoʲ)
        optimize!(m)

        # Check to see what the termination status (i.e. result) of the optimization was 
        status = termination_status(m)
        if status == MOI.OPTIMAL         # A feasible solution was found ⇒ trivial!
            push!(trivial_idxs, j)
        elseif status == MOI.INFEASIBLE  # No feasible solution exists    ⇒ fragile!
            push!(fragile_idxs, j)
        else
            throw("Unexpected termination status $status")
        end
    end
    nsᴴ_trivial = @view sb_nontopo[trivial_idxs]
    nsᴴ_fragile = @view sb_nontopo[fragile_idxs]

    return nsᴴ_trivial, nsᴴ_fragile
end

"""
    _test_hilbert_bases_consistency(BRS::BandRepSet, F::SmithNormalForm.Smith,
                    nsᴴ::AbstractMatrix, nsᴴ_nontopo::AbstractMatrix, zsᴴ::AbstractMatrix)

Test that the obtained "non-topological" bases indeed obey some of the conditions that we
know they must. Prints a checkmark (✓) if succesful; throws `Test.FallbackTestSetException`
otherwise. Returns `nothing`.
"""
function _test_hilbert_bases_consistency(BRS::BandRepSet, F::SmithNormalForm.Smith,
                nsᴴ::AbstractMatrix, nsᴴ_nontopo::AbstractMatrix, zsᴴ::AbstractMatrix)

    print("   ... checking consistency of non-topological bases: ")

    dᵇˢ        = count(!iszero, F.SNF)  # "Dimensionality" of band structure
    Nᴴ         = size(nsᴴ, 2)           # Number of compatibility Hilbert bases
    Nᴴ_nontopo = size(nsᴴ_nontopo, 2)   # Number of nontopological Hilbert bases

    # All zsᴴ that are not element divisible by Λ ≡ F.SNF correspond to "proper" topological
    # cases (i.e. not fragile). This is because those z ≡ Λy with yᵢ∈ℚ\ℤ are topological, 
    # whereas all those with yᵢ∈ℤ are either trivial or fragilely topological.
    # Note: The approach below is not valid in general: while it does find all the 
    #       non-topological elements among the Hilbert bases nsᴴ, it does not find a full
    #       Hilbert basis for all non-topological states. The easiest way to see this is to
    #       think in terms of a "unit cell" for the Hilbert bases, and then realize that 
    #       this unit cell may be far larger, when we take a subcone of the original cone.
    Λ = @view F.SNF[OneTo(dᵇˢ)] # Nonzero invariant factors of Smith normal decomposition
    nontopo_idxs_subset = findall(zᴴ -> all(zᴴΛᵢ -> mod(zᴴΛᵢ[1], zᴴΛᵢ[2]) == 0, zip(zᴴ, Λ)),
                           collect(eachcol(zsᴴ)))
    topo_idxs_subset   = findall(i -> i ∉ nontopo_idxs_subset, OneTo(Nᴴ))
    nsᴴ_nontopo_subset = @view nsᴴ[:, nontopo_idxs_subset]
    nsᴴ_topo_subset    = @view nsᴴ[:, topo_idxs_subset]

    # If classification is Z₁, nsᴴ and nsᴴ_nontopo must be equal
    if classification(BRS) == "Z₁"
        @test Set(eachcol(nsᴴ)) == Set(eachcol(nsᴴ_nontopo))
        println("✓ (trivially)")
    else 
        # Must be a superset of the naive extraction approach
        @test Set(eachcol(nsᴴ_nontopo_subset)) ⊆ Set(unique(eachcol(nsᴴ_nontopo)))

        # Check that every basis vector in nsᴴ_nontopo can be expanded in the compatibility 
        # basis nsᴴ using only non-negative integer coefficients
        if Nᴴ < Nᴴ_nontopo # no need to test unless there are more elements in nsᴴ_nontopo
            for (j,nsᴴ_nontopoʲ) in enumerate(eachcol(nsᴴ_nontopo))
                j ≠ 1 && print(stdout, "\b"^ndigits(j-1))
                print(stdout, j)
                flush(stdout)

                m = Model(GLPK.Optimizer)
                @variable(m, c[1:Nᴴ] >= 0, Int)
                @constraint(m, nsᴴ*c .== nsᴴ_nontopoʲ)
                optimize!(m)

                cvals = Int.(value.(c)) # extract optimized expansion coefficients
                @test nsᴴ*cvals == nsᴴ_nontopoʲ
                @test all(cvals .≥ 0)
            end
        end
        print("\b"^ndigits(Nᴴ_nontopo), "✓", " "^ndigits(Nᴴ_nontopo-1))
    end

    nothing
end


# Footnotes:
# ¹⁾ We want to control the algorithm used to calculate the Hilbert basis in Normaliz: I've
# found that Normaliz usually picks the "PrimalMode" algorithm (not always though), despite
# the documentation stating that "DualMode" usually is better for cones defined by
# inequalities. Through tests, I've found that there is usually a speedup of order 2-10,
# if the "DualMode" algorithm is picked instead; occassionally, up to ×100s. The speedup 
# seems to be especially large for the systems that (very) hard to solve (e.g. 131). 
# To force "DualMode", we use the method C.Compute(<quantity-to-compute>, <algorithm>), 
# which then calculates <quantity-to-compute> and stores it in C (here, "HilbertBasis");
# it can then subsequently be retrieved from C with C.HilbertBasis().

end # module