using Crystalline, Test

@testset "Crystalline" begin
    # basic symmetry operations
    include("symops.jl")
    include("SquareStaticMatrices.jl")
    include("groups_xyzt_vs_coded.jl")
    include("generators_xyzt_vs_coded.jl")

    # Bravais.jl
    include("niggli.jl")
    include("basisvecs.jl")

    # abstractvecs
    include("kvecs.jl")
    include("wyckoff.jl")

    # symmetry vectors
    include("symmetryvectors.jl")

    # show, notation, and cached info
    include("show.jl")
    include("notation.jl")
    include("orders.jl")

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
    include("irreps_physical_reality.jl")
    include("lgirreps_vs_pgirreps_at_Gamma.jl")

    # compatibility
    include("compatibility.jl")

    # conjugacy classes and irreps
    include("conjugacy.jl")

    # lattices
    include("fourierlattices.jl")

    # additional k-vectors in Φ-Ω ("special" representation domain vectors)
    include("holosymmetric.jl")

    # band representations & site symmetry groups
    include("classification.jl")  # => does topo classification agree w/ Adrian?
    include("bandrep.jl")         # => do k-vectors match (Bilbao's bandreps vs ISOTROPY)?
    include("calc_bandreps.jl")   # => tests of /src/calc_bandreps.jl
    include("isomorphic_parent_pointgroup.jl") # => do we find the same "parent" point
                                               #    groups as Bilbao?

    # analysis of symmetry eigenvalue data (to symmetry vectors)
    include("symeigs_analysis.jl")
    
    # magnetic space groups
    include("mspacegroup.jl")

    # group relations (sub- & supergroups)
    include("grouprelations.jl")
end