using SmithNormalForm
using BenchmarkTools
using SparseArrays
using Random

Random.seed!(2903872398473);
rows, cols = 64, 93
TMP = rand(1:10, rows*cols)
D = reshape(TMP, rows, cols)
i, j = 2, 3
a, b, c, d = 1, 2, 3, 4

suite = BenchmarkGroup()
suite["internal"] = BenchmarkGroup(["swap", "elimination"])

suite["internal"]["col_swap"] = @benchmarkable SmithNormalForm.cswap!(D, i, j) gctrial=true evals=1000 samples=100000
suite["internal"]["row_swap"] = @benchmarkable SmithNormalForm.rswap!(D, i, j) gctrial=true evals=1000 samples=100000

suite["internal"]["col_elim"] = @benchmarkable SmithNormalForm.colelimination(D, a, b, c, d, i, j) evals=500 gctrial=true
suite["internal"]["row_elim"] = @benchmarkable SmithNormalForm.rowelimination(D, a, b, c, d, i, j) evals=500 gctrial=true

suite["snf"] = BenchmarkGroup(["type"])
Z = rand(1:100, 100, 120)
Z[Z .> 66] .= 1
Z[Z .> 33] .= -1
Z[Z .> 0] .= 0
suite["snf"]["dense"] = @benchmarkable snf(Z)

Y = sprand(Int, 100, 120, 0.3)
Y[Y .< 0] .= -1
Y[Y .> 0] .= 1
suite["snf"]["sparse"] = @benchmarkable snf(Y)

# tune!(suite);
# BenchmarkTools.save("bmsetup.json", params(suite));
loadparams!(suite, BenchmarkTools.load("bmsetup.json")[1], :evals, :samples);
results = run(suite, verbose = true, seconds = 1)
BenchmarkTools.save("results.json", results)
