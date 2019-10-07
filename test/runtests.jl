using SGOps, Test

@testset "All tests" begin
    # basic symmetry operations
    include("symops.jl")

    # group checks
    include("littlegroup_orders.jl")

    # multiplication tables and irreps
    include("irreps_orthogonality.jl")
    include("multtable.jl")
    include("reality.jl")

    # band representations
    include("classification.jl")  # does topo classification agree w/ Adrian?
    include("bandrep.jl")         # do k-vectors match (BandRep vs ISOTROPY)?
end
