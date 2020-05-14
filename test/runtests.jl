using Crystalline, Test

@testset "Crystalline" begin
    # basic symmetry operations
    include("basisvecs.jl")
    include("symops.jl")

    # group checks
    include("littlegroup_orders.jl")
    include("pointgroup.jl")

    # loading irreps from .jld files vs parsing of ISOTROPY
    include("parsed_vs_loaded_littlegroup_irreps.jl")

    # multiplication tables and irreps
    include("irreps_orthogonality.jl")
    include("chartable.jl")
    include("multtable.jl")
    include("irreps_reality.jl")

    # additional k-vectors in Φ-Ω ("special" representation domain vectors)
    include("holosymmetric.jl")

    # band representations
    include("classification.jl")  # does topo classification agree w/ Adrian?
    include("bandrep.jl")         # do k-vectors match (Bilbao's bandreps vs ISOTROPY)?
end