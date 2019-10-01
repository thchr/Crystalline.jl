using SGOps, Test

@testset "All tests" begin
    include("symops.jl")

    include("irreps.jl")

    include("multtable.jl")
    include("reality.jl")
end
