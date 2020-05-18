using JuMP, GLPK, Crystalline

# see documentation for JuMP at http://www.juliaopt.org/JuMP.jl/v0.21.1
sgnum = 230
BRS = bandreps(sgnum)       # 
N = length(BRS)             # number of BandReps
B = matrix(BRS)             # "atomic" basis
B′ = Crystalline.wyckbasis(BRS)[1] # compat. rel. basis

# pick a symmetry vector out of a hat
c₀ = zeros(Int64, N)
c₀[1] = 1; c₀[3] = 2; c₀[5] = 4; c₀[end] = 1
#c₀ = [0,-1,0,-2,0,-4,7,7]; # for sg #3 this is an equivalent symmetry indicator to c₀=[1,0,3,0,4,0,0,0]
n = B*c₀

println("c_orig  = ", c₀)

# feasibility problem with positive-integer variables, subject to the condition Aᵀc = n
m = Model(GLPK.Optimizer)
@variable(m, c[1:N] >= 0, Int)
@constraint(m, B*c .== n)

# solve the model and print solution
optimize!(m)
println("c_optim = [", join(Int.(value.(c)), ", "), "]")

# try with Nemo.jl instead
using Nemo
B′ = Nemo.matrix(ZZ, size(B)..., B)
n′ = MatrixSpace(ZZ,size(n,1),1)(n)
_, x′ = cansolve(B′, n′)

println("c_nemo  = ", Int.(vec(Matrix(x′)))) # in general will not get a nice factorization from this