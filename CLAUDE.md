# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Crystalline.jl is a Julia package providing tools for crystalline symmetry analysis, including:
- Symmetry operations for crystalline point groups, space groups, and subperiodic groups
- Wyckoff positions and their irreducible representations
- Band representations for topological classification
- Bravais lattice utilities through the embedded Bravais.jl subpackage
- Fourier lattices and magnetic space groups

## Commands

### Testing
```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

### Building Documentation
```bash
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
julia --project=docs docs/make.jl
```

### Running Individual Tests
```bash
julia --project=. -e "using Pkg; Pkg.test(); include(\"test/[testfile].jl\")"
```

### Package Development Setup
The repository includes a local Bravais.jl subpackage that must be developed locally:
```bash
julia --project=. -e "using Pkg; Pkg.develop(PackageSpec(path=\"Bravais\"))"
```

## Architecture

### Core Components

- **src/types.jl**: Core types (`SymOperation`, `SpaceGroup`, `PointGroup`, `LittleGroup`, `LGIrrep`, `BandRep`, etc.)
- **src/types_symmetryvectors.jl**: `SymmetryVector` and `SymmetryVectors` types for band representation analysis
- **SquareStaticMatrices.jl**: Submodule providing optimized static matrices for symmetry operations
- **JaggedVector**: Submodule for memory-efficient vectors of vectors

### Symmetry Operations
- **src/symops.jl**: Core symmetry operation functionality
- **src/tables/**: Precomputed symmetry operation tables for various groups
- **src/assembly/**: Group assembly from generators

### Irreducible Representations
- **src/littlegroup_irreps.jl**: Little group irreps loaded from JLD2 data files
- **src/pointgroup.jl**: Point group irreps
- **src/irreps_reality.jl**: Reality conditions for irreps
- **src/irreps_physical_reality.jl**: Physical reality considerations

### Band Representations & Topology  
- **src/bandrep.jl**: Elementary band representations
- **src/calc_bandreps.jl**: Band representation calculations
- **src/tqc_analysis.jl**: Topological quantum chemistry analysis
- **src/compatibility.jl**: Compatibility relations between k-points

### Data Files
- **data/**: Extensive precomputed crystallographic data (excluded from analysis per user request)
  - Irreducible representations stored in JLD2 format
  - Band representations and Wyckoff positions in CSV format
  - Data loaded lazily via `__init__()` function

### Subpackages
- **Bravais/**: Lightweight standalone package for Bravais lattice utilities
- **Vendored SmithNormalForm.jl**: Located in `.vendor/` directory

### Extensions
- **CrystallineMakieExt.jl**: Makie.jl plotting extensions
- **CrystallineGraphMakieExt.jl**: GraphMakie.jl extensions for visualizing symmetry relationships

## Development Notes

### Dependencies
- Core: LinearAlgebra, StaticArrays, JLD2, PrettyTables
- Optional plotting: Makie.jl, GraphMakie.jl (via extensions)
- Crystallographic data loaded from JLD2 files at package initialization

### Testing Strategy
- Comprehensive test suite covering all major functionality
- Tests compare against ISOTROPY and Bilbao databases
- Cross-validation between different computational methods
- Coverage includes edge cases and mathematical consistency checks

### CI/CD
- GitHub Actions for testing across Julia versions 1.10+ and multiple OS
- Documentation automatically built and deployed via Documenter.jl
- Coverage tracking via Codecov