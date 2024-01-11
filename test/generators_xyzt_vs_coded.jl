# ---------------------------------------------------------------------------------------- #
module CrystallineTestGeneratorsXYZT
# Module with functionality to load generators from stored xyzt-data .csv files
# (this was how the data was loaded in Crystalline ≤v0.4.22)

using Crystalline
using Crystalline: 
    _unmangle_pgiuclab, pointgroup_num2iuc, # point groups
    _check_valid_sgnum_and_dim, # space group
    _subperiodic_kind, _check_valid_subperiodic_num_and_dim, # subperiodic groups
    _throw_invalid_dim # general

export generators_from_xyzt

const XYZT_GEN_DATA_DIR = joinpath(@__DIR__, "data", "xyzt", "generators")

# ---------------------------------------------------------------------------------------- #
# Point groups
# ---------------------------------------------------------------------------------------- #
function read_pggens_xyzt(iuclab::String, D::Integer)
    @boundscheck D ∉ (1,2,3) && _throw_invalid_dim(D)
    @boundscheck iuclab ∉ PG_IUCs[D] && throw(DomainError(iuclab, "iuc label not found in database (see possible labels in PG_IUCs[D])"))

    filepath = joinpath(XYZT_GEN_DATA_DIR, "pgs/"*string(D)*"d/"*_unmangle_pgiuclab(iuclab)*".csv")

    return readlines(filepath)
end

function generators_from_xyzt(iuclab::String, ::Type{PointGroup{D}}=PointGroup{3}) where D
    ops_str = read_pggens_xyzt(iuclab, D)
    return SymOperation{D}.(ops_str)
end
function generators_from_xyzt(pgnum::Integer, ::Type{PointGroup{D}}, setting::Integer=1) where D
    iuclab = pointgroup_num2iuc(pgnum, Val(D), setting)
    return generators_from_xyzt(iuclab, PointGroup{D})
end

# ---------------------------------------------------------------------------------------- #
# Space groups
# ---------------------------------------------------------------------------------------- #

function generators_from_xyzt(sgnum::Integer, ::Type{SpaceGroup{D}}=SpaceGroup{3}) where D
    ops_str = read_sggens_xyzt(sgnum, D)

    return SymOperation{D}.(ops_str)
end

function read_sggens_xyzt(sgnum::Integer, D::Integer)
    @boundscheck _check_valid_sgnum_and_dim(sgnum, D)

    filepath = joinpath(XYZT_GEN_DATA_DIR, "sgs/"*string(D)*"d/"*string(sgnum)*".csv")

    return readlines(filepath)
end

# ---------------------------------------------------------------------------------------- #
# Subperiodic groups
# ---------------------------------------------------------------------------------------- #

function generators_from_xyzt(num::Integer, ::Type{SubperiodicGroup{D,P}}) where {D,P}
    ops_str = read_subperiodic_gens_xyzt(num, D, P)

    return SymOperation{D}.(ops_str)
end

function read_subperiodic_gens_xyzt(num::Integer, D::Integer, P::Integer)
    @boundscheck _check_valid_subperiodic_num_and_dim(num, D, P)

    kind = _subperiodic_kind(D, P)
    filepath = joinpath(XYZT_GEN_DATA_DIR, "subperiodic/"*kind*"/"*string(num)*".csv")

    return readlines(filepath)
end
end # module CrystallineTestGeneratorsXYZT

# ---------------------------------------------------------------------------------------- #

using Crystalline, Test
using .CrystallineTestGeneratorsXYZT

@testset "Agreement between \"coded\" group generators and xyzt-data" begin
    @testset "Point groups" begin
        for D in 1:3
            Dᵛ = Val(D)
            for (pgnum, pgiucs) in enumerate(Crystalline.PG_NUM2IUC[D])
                for (setting, pgiuc) in enumerate(pgiucs)
                    pg_gens = generators(pgiuc, PointGroup{D})
                    pg_gens_from_xyzt = generators_from_xyzt(pgiuc, PointGroup{D})
                    pg_gens′ = generators(pgnum, PointGroup{D}, setting)
                    @test pg_gens ≈ pg_gens_from_xyzt
                    @test pg_gens == pg_gens′
                end
            end
        end
    end

    @testset "Space groups" begin
        for D in 1:3
            Dᵛ = Val(D)
            for sgnum in 1:MAX_SGNUM[D]
                sg_gens = generators(sgnum, SpaceGroup{D})
                sg_gens_from_xyzt = generators_from_xyzt(sgnum, SpaceGroup{D})
                @test sg_gens ≈ sg_gens_from_xyzt
            end
        end
    end

    @testset "Subperiodic groups" begin
        for (D,P) in ((3, 2), (3, 1), (2, 1))
            Dᵛ, Pᵛ = Val(D), Val(P)
            for num in 1:MAX_SUBGNUM[(D,P)]
                subg_gens = generators(num, SubperiodicGroup{D,P})
                subg_gens_from_xyzt = generators_from_xyzt(num, SubperiodicGroup{D,P})
                @test subg_gens ≈ subg_gens_from_xyzt
            end
        end
    end
end
