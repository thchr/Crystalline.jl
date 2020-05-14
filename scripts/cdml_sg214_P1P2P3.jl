using Crystalline, LinearAlgebra, Test

sgnum = 214
kv = KVec([0.5,0.5,0.5])
k = kv()
cntr = centering(214)

# Constants from CDML
C = 1/√2                            # scalars
I = 1.0im
P = cis(1*π/12)
Q = cis(5*π/12)
U = cis(1*π/6)
V = cis(1*π/4)
W = cis(2*π/3)

M1  = ComplexF64[1.0 0.0; 1.0 0.0]  # matrices
M2  = ComplexF64[0.0 1.0; 1.0 0.0]
M4  = ComplexF64[1.0 0.0; 0.0 -1.0]
M16 = C*ComplexF64[P conj(Q); P -conj(Q)]

# Note that CDML has the same setting for 214 as ITA/Bilbao/ISOTROPY does; no need 
# to transform any of the operators
R1  = SymOperation{3}("x,y,z")         # rotations
R2  = SymOperation{3}("x,-y,-z")
R3  = SymOperation{3}("-x,y,-z")
R4  = SymOperation{3}("-x,-y,z")
R5  = SymOperation{3}("y,z,x")
R6  = SymOperation{3}("y,-z,-x")
R7  = SymOperation{3}("-y,z,-x")
R8  = SymOperation{3}("-y,-z,x")
R9  = SymOperation{3}("z,x,y")
R10 = SymOperation{3}("z,-x,-y")
R11 = SymOperation{3}("-z,x,-y")
R12 = SymOperation{3}("-z,-x,y")

T1  = SymOperation{3}("x,y+1/2,z")     # translations
T2  = SymOperation{3}("x,y,z+1/2")
T3  = SymOperation{3}("x+1/2,y,z")

# Generators (elements and irreps) at k-point P in sgnum = 214
opsgen = [R4∘T1, R2∘T2, R9]; 
irsgen = ([M4, M2, M16],         # P1
          [M4, M2, M16*W],       # P2
          [M4, M2, M16*conj(W)]) # P3

iridx = 2
# from the _generators_ compute all the elements and irreps of
# the irrep characterized by irsgen and opsgen
ops = deepcopy(opsgen)
irs = deepcopy(irsgen[iridx])
while true
    added = false
    Nₒₚ = length(ops)
    for row = 1:Nₒₚ
        for col =  1:Nₒₚ
            newop = ops[row]∘ops[col]
            # only want to include elements that don't include trivial primitive translations
            if Crystalline.conventionalize(primitivize(newop, cntr), cntr) != newop; continue; end 
            # create new irrep from multiplication properties of group elements
            if newop ∉ ops # only add if not already created
                push!(ops, newop)
                push!(irs, irs[row]*irs[col])

                t₀ = translation(ops[row]) + rotation(ops[row])*translation(ops[col]) - translation(newop)
                ϕ =  2π*dot(k,t₀) # phase is nonzero for ray irreps
                display(t₀)
                irs[end] .*= cis(-ϕ)
                for (idx, val) in enumerate(irs[end]) # small near-zero values are annoying, fix manually
                    if abs(real(val)) < 1e-12
                        irs[end][idx] = complex(0.0, imag(val))
                    elseif abs(imag(val)) < 1e-12
                        irs[end][idx] = complex(real(val), 0.0)
                    end
                end
                added = true
            end
        end
    end
    if !added; break; end
end
primitive_ops = primitivize.(ops, cntr)
Nₒₚ = length(ops)
dimirr = size(first(irs),1)

# populate our LGIrrep struct with computed values
lgir = LGIrrep(sgnum, "P"*string(iridx), kv, ops, irs, [zeros(3) for _=1:length(ops)], 1)

# check multiplication table (operator and irreps)
mt_ops = multtable(primitive_ops)
mt_irs = Crystalline.checkmulttable(mt_ops, lgir)
println("Irrep mult-table check: ", sum(mt_irs), "/", prod(size(mt_irs)), '\n')

# check great orthogonality theorem
kroncheck = sum(kron(conj.(irs[i]), irs[i]) for i in Base.OneTo(length(irs)))
display(kroncheck)
println("   |G|/l = $(Nₒₚ/dimirr)")

##
# MATCH OPERATOR SORTING TO THAT IN ISOTROPY
if !isdefined(Main, :LGIRS)
    LGIRS = parselittlegroupirreps()
end
lgir_iso = LGIRS[sgnum][7][iridx] # 7 is idx of P point

ops_iso_str = xyzt.(operations(lgir_iso))
perm = Vector{Int64}(undef, Nₒₚ)
for i = 1:Nₒₚ
    perm[i] = findfirst(xyzt(ops[i]) .== ops_iso_str)
end
ops′ = ops[invperm(perm)]
irs′ = irs[invperm(perm)]

if iridx == 2 # try to guess at appropriate transform between cdml and iso (only for P2)
    σx, σy = [0 1; 1 0], [0 -im; im 0]    
    T, S = eigen(σx).vectors, eigen(-σy).vectors
    P = S*inv(T)
    irs′ = [inv(P)*ir*P for ir in irs′]
end
@test xyzt.(ops′) == ops_iso_str
lgir_cdml = LGIrrep(sgnum, "P"*string(iridx), kv, ops′, irs′, [zeros(3) for _=1:length(ops)], 1)

display(irreps(lgir_cdml) .≈ irreps(lgir_iso))

display(characters(lgir_cdml) .≈ characters(lgir_iso))