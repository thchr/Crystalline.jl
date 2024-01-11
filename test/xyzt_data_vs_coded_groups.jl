# ---------------------------------------------------------------------------------------- #
module CrystallineTestXYZT
# Module with functionality to load `SymOperation`s from stored xyzt-data .csv files
# (this was how the data was loaded in Crystalline ≤v0.4.22)

using Crystalline
using Crystalline: 
    unmangle_pgiuclab, pointgroup_num2iuc, pointgroup_iuc2num, # point groups
    _check_valid_sgnum_and_dim, # space group
    subperiodic_kind, _check_valid_subperiodic_num_and_dim, # subperiodic groups
    _throw_invalid_dim # general

export pointgroup_from_xyzt, spacegroup_from_xyzt, subperiodicgroup_from_xyzt

const OPERATIONS_DATA_DIR = joinpath(@__DIR__, "data", "xyzt-operations")
# ---------------------------------------------------------------------------------------- #
# Point groups
# ---------------------------------------------------------------------------------------- #

@inline function pointgroup_from_xyzt(iuclab::String, Dᵛ::Val{D}=Val(3)) where D
    pgnum = pointgroup_iuc2num(iuclab, D) # this is not generally a particularly well-established numbering
    ops_str = read_pgops_xyzt(iuclab, D)
    
    return PointGroup{D}(pgnum, iuclab, SymOperation{D}.(ops_str))
end

# Here, and below for space groups and subperiodic groups, we obtain the symmetry operations
# via a set of stored .csv files listing the operations in xyzt format for each group: the
# data is stored in `test/data/xyzt-operations/`
function read_pgops_xyzt(iuclab::String, D::Integer)
    @boundscheck D ∉ (1,2,3) && _throw_invalid_dim(D)
    @boundscheck iuclab ∉ PG_IUCs[D] && throw(DomainError(iuclab, "iuc label not found in database (see possible labels in PG_IUCs[D])"))
    filepath = joinpath(XYZT_DATA_DIR, "pgs/"*string(D)*"d/"*unmangle_pgiuclab(iuclab)*".csv")

    return readlines(filepath)
end

# ---------------------------------------------------------------------------------------- #
# Space groups
# ---------------------------------------------------------------------------------------- #

function spacegroup_from_xyzt(sgnum::Integer, ::Val{D}=Val(3)) where D
    ops_str = read_sgops_xyzt(sgnum, D)
    ops = SymOperation{D}.(ops_str)

    return SpaceGroup{D}(sgnum, ops)
end

function read_sgops_xyzt(sgnum::Integer, D::Integer)
    @boundscheck _check_valid_sgnum_and_dim(sgnum, D)
    filepath = joinpath(XYZT_DATA_DIR, "sgs/"*string(D)*"d/"*string(sgnum)*".csv")

    return readlines(filepath)
end

# ---------------------------------------------------------------------------------------- #
# Subperiodic groups
# ---------------------------------------------------------------------------------------- #
function subperiodicgroup_from_xyzt(num::Integer, 
                                  ::Val{D}=Val(3), ::Val{P}=Val(2)) where {D,P}
    ops_str = read_subperiodic_ops_xyzt(num, D, P)
    ops = SymOperation{D}.(ops_str)

    return SubperiodicGroup{D,P}(num, ops)
end

function read_subperiodic_ops_xyzt(num::Integer, D::Integer, P::Integer)
    @boundscheck _check_valid_subperiodic_num_and_dim(num, D, P)
    kind = subperiodic_kind(D, P)
    filepath = joinpath(XYZT_DATA_DIR, "subperiodic/"*kind*"/"*string(num)*".csv")

    return readlines(filepath)
end

# ---------------------------------------------------------------------------------------- #
# Magnetic space groups
# ---------------------------------------------------------------------------------------- #
# Maybe TODO

end # module CrystallineTestXYZT

# ---------------------------------------------------------------------------------------- #

using Crystalline, Test
using .CrystallineTestXYZT

@testset "Agreement between \"coded\" group operations and xyzt-data" begin
    @testset "Point groups" begin
        for D in 1:3
            Dᵛ = Val(D)
            for pgiuc in PG_IUCs[D]
                pg = pointgroup(pgiuc, Dᵛ)
                pg_from_xyzt = pointgroup_from_xyzt(pgiuc, Dᵛ)
                @test all(pg .≈ pg_from_xyzt)
            end
        end
    end

    @testset "Space groups" begin
        for D in 1:3
            Dᵛ = Val(D)
            for sgnum in 1:MAX_SGNUM[D]
                sg = spacegroup(sgnum, Dᵛ)
                sg_from_xyzt = spacegroup_from_xyzt(sgnum, Dᵛ)
                @test all(sg .≈ sg_from_xyzt)
            end
        end
    end

    @testset "Subperiodic groups" begin
        for (D,P) in ((3, 2), (3, 1), (2, 1))
            Dᵛ, Pᵛ = Val(D), Val(P)
            for num in 1:MAX_SUBGNUM[(D,P)]
                subg = subperiodicgroup(num, Dᵛ, Pᵛ)
                subg_from_xyzt = subperiodicgroup_from_xyzt(num, Dᵛ, Pᵛ)
                @test all(subg .≈ subg_from_xyzt)
            end
        end
    end
end