# Derive the basis vectors from the compatibility relations in SG 2, following Song's 2018 PRX

# The only compatibility relation is n⁺(kⱼ) + n⁻(kⱼ) - n⁺(kᵢ) - n⁻(kᵢ) = 0 for i, j varying
# over the distinct high-symmetry k-points kⱼ. We include the trivial combinations i = j, 
# since it doesn't change the outcome and makes the implementation simpler. 
# These relations reflect the required "invariance" of the occupation number within the BZ,
# here sampled at points with definite symmetry eigenvalues (i.e. the high-symmetry points).

# For space group 2, there are 8 high symmetry points kⱼ = (α, β, γ) corresponding to the 
# combinations of (α, β, γ) = ({0,½}, {0,½}, {0,½}).
# We find the solution space to these equations by rewriting them as a matrix equation
using Crystalline, LinearAlgebra, Nemo

Nk = 8
Nirr = 2*Nk
C = zeros(Int, Nk^2, Nirr) # 8*8 compatibility relation combinations & 16 irreps (kᵢ⁺ and kᵢ⁻)

# C is a matrix such that when acting on a valid symmetry vector 𝐧, we have C𝐧 = 0.
for r in 1:Nk
    for p in 1:Nk
        global C
        C[(r-1)*Nk+p, (2*r-1):2r] .+= 1
        C[(r-1)*Nk+p, (2*p-1):2p] .-= 1
    end
end


# Get nullspace of C with integer coefficients; i.e. find the nullspace in a field of 
# integers ℤ. 
# We can use the Smith Normal Form: For A = SΛT, we can obtain the nullspace of A from the
# last n columns of T⁻¹ with n denoting the number of zeros in Λ [i.e. this is n=nullity(A); 
# contrast this with r=rank(A). Note r+n=dim(V), cf. the rank-nullity theorem, where V is
# the mapping domain in A: V→W; since A is a matrix this is simply the number of columns of 
# A]. See e.g. https://core.ac.uk/download/pdf/82343294.pdf regarding the Smith normal form
# and its application to null spaces.
F = Crystalline._smith′(C) # Smith Normal Form (small wrapper around `SmithNormalForm.smith(..)`)
S, S⁻¹, T, T⁻¹, Λ = F.S, F.Sinv, F.T, F.Tinv, F.SNF
r = count(!iszero, Λ) # number of nonzeros in Smith normal diagonal matrix = rank(C)
zidxs  = r+1:length(Λ)
basis = T⁻¹[:, zidxs]
display(basis)

# This can also be done directly with Nemo.jl, using Nemo.nullspace (see documentation at
# http://nemocas.github.io/Nemo.jl/latest/)

Cᴺᴱᴹᴼ = Nemo.matrix(ZZ, size(C)..., C) # Convert C from ::Matrix{Int} to Nemo's ::fmpz_mat 
            # type, which e.g. contains info regarding the kind of element field (here, `ZZ
            # = Integer Ring =` ℤ)
n′ᴺᴱᴹᴼ, basis′ᴺᴱᴹᴼ = Nemo.nullspace(Cᴺᴱᴹᴼ) # get null space with element type `ZZ` = ℤ
basis′ = Matrix{Int}(basis′ᴺᴱᴹᴼ) # Convert back to a standard Julia matrix
display(basis′)

# We can verify that `basis` and `basis′` span the same space by testing that we can expand
# every column of `basis′` in `basis` and vice versa. We can use ordinary solve (\) or 
# Nemo's 
basisᴺᴱᴹᴼ = Nemo.matrix(ZZ, size(basis)..., basis)
_, x = Nemo.cansolve(basisᴺᴱᴹᴼ, basis′ᴺᴱᴹᴼ)  # Try to write every element of basis′ᴺᴱᴹᴼ as a linear combination of basis′ᴺᴱᴹᴼ
_, x′ = Nemo.cansolve(basis′ᴺᴱᴹᴼ, basisᴺᴱᴹᴼ) # ... vice versa
println("Equivalent span check 1: ", basisᴺᴱᴹᴼ*x == basis′ᴺᴱᴹᴼ)
println("Equivalent span check 2: ", basis′ᴺᴱᴹᴼ*x′ == basisᴺᴱᴹᴼ)

# Now we ought to compare this against the basis obtained from Crystalline.wyckbasis(...)
basis′′ = Crystalline.wyckbasis(bandreps(2))[1]
basis′′ᴺᴱᴹᴼ = Nemo.matrix(ZZ, size(basis′′)..., basis′′)
_, x′′ =  Nemo.cansolve(basisᴺᴱᴹᴼ, basis′′ᴺᴱᴹᴼ)
println("Crystalline.jl span check: ", basisᴺᴱᴹᴼ*x′′ == basis′′ᴺᴱᴹᴼ)


# Any valid symmetry indicator vector 𝐧 = [n(Γ⁺), n(Γ⁻), n(X⁺), n(X⁻), ..., n(Z⁺), n(Z⁻)] 
# must be spannable by `basis`, i.e. there must exist integer coefficients cᵢ such that 
#       𝐧 = c₁𝐫₁ + c₂𝐫₂ + ... + c₉𝐫₉
# where 𝐫ᵢ denotes the rows of `basis`. In practice, this means that if we have a vector 𝐧,
# we can test whether it fulfils the compatibility relations by seeing if there is a 
# solution [𝐫₁ ... 𝐫₉]𝐜 = 𝐧, i.e. by calling `𝐜 = basis\𝐧` and then subsequently checking if
# that indeed is a solution, e.g. via Nemo.cansolve(basis′, 𝐧). If we want to look for a 
# solution with fractional coefficients, we might promote basis and 𝐧 to QQ ≡ ℚ and then
# attempt a `cansolve` or maybe even simply `solve_rational` call?