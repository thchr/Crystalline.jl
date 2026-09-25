# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Crystalline.jl (v0.6.28, Julia ≥ 1.12) is a Julia package for programmatic access to and manipulation of crystalline symmetry. The author is Thomas Christensen. Core capabilities:

- Symmetry operations and groups for space groups (1D/2D/3D), point groups, site symmetry groups, subperiodic groups, and magnetic space groups
- Irreducible representations (irreps) of little groups (`LGIrrep`), point groups (`PGIrrep`), and site symmetry groups (`SiteIrrep`), loaded from precomputed JLD2 data files
- Wyckoff positions and site symmetry groups
- Band representations and topological quantum chemistry (TQC) analysis
- Compatibility relations between k-points
- Group–subgroup relations (maximal subgroups, minimal supergroups)
- Symmetry vectors and their TQC analysis
- Bravais lattice utilities (factored into the Bravais.jl subpackage)
- Fourier lattices for level-set models (photonic crystals)

Primary external data sources: **ISOTROPY** dataset (ISO-IR) for little group irreps; **Bilbao Crystallographic Server** for band representations and supplementary data.

## Commands

### Testing
```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

### Running individual test files
```bash
julia --project=. test/[testfile].jl
```

### Building documentation
```bash
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
julia --project=docs docs/make.jl
```

### Package development setup
The repository includes a local Bravais.jl subpackage that must be developed locally first:
```bash
julia --project=. -e "using Pkg; Pkg.develop(PackageSpec(path=\"Bravais\"))"
```

## Architecture

### Source file map

| File | Contents |
|------|----------|
| `src/Crystalline.jl` | Module entry point, imports, `__init__` (opens JLD2 files), exports |
| `src/types.jl` | All core types (see below) |
| `src/types_symmetryvectors.jl` | `SymmetryVector`, `SymmetryVectors`, `BandRep`, `CompositeBandRep`, `@composite` |
| `src/symops.jl` | String↔matrix (`@S_str`, `xyzt2components`), `compose`, `littlegroup`, `orbit`, `reduce_ops`, `issubgroup`, `cosets` |
| `src/littlegroup_irreps.jl` | `lgirreps`, `littlegroups` (load from JLD2); `(lgir::LGIrrep)(αβγ)` evaluator with phase factors; `israyrep` |
| `src/pointgroup.jl` | `pgirreps`, `find_isomorphic_parent_pointgroup` |
| `src/wyckoff.jl` | `wyckoffs`, `WyckoffPosition`, `SiteGroup`, `sitegroups`, `siteirreps`, `findmaximal` |
| `src/bandreps.jl` | `bandreps` (compute EBRs from site irreps; Cano et al. PRB 97, 035139 (2018)) |
| `src/tqc_analysis.jl` | `calc_topology`, `iscompatible`, `symmetry_indicators`, `indicator_group`, `basisdim`, `TopologyKind` (TRIVIAL/NONTRIVIAL/FRAGILE) |
| `src/compatibility.jl` | `subduction_count`, `remap_to_kstar` |
| `src/irreps_reality.jl` | `calc_reality`, `realify`, `realify!` (Herring criterion) |
| `src/irreps_physical_reality.jl` | `physical_realify`, `timereversal_unitary` (canonical form of co-representations under time reversal) |
| `src/grouprelations/` | `maximal_subgroups`, `minimal_supergroups`, `conjugacy_relations` |
| `src/fourierlattices.jl` | `ModulatedFourierLattice`, `levelsetlattice`, `modulate`, `filling2isoval`, etc. |
| `src/symeigs2irrep.jl` | `find_representation` (identify irrep from symmetry eigenvalues) |
| `src/symeigs_analysis.jl` | `collect_compatible`, `collect_irrep_annotations` |
| `src/conjugacy.jl` | `classes`, `is_abelian` |
| `src/notation.jl` | `schoenflies`, `iuc`, `seitz`, `mulliken` |
| `src/subperiodic.jl` | `SubperiodicGroup` (layer, rod, frieze groups) |
| `src/magnetic/` | `MSymOperation`, `MSpaceGroup` (type-IV magnetic space groups) |
| `src/double/` | `SU2`, `DSymOperation`, the `D`-prefixed groups and irreps, `su2`, `isbarred`, `doublegroup`, `isspinful` (double groups and double-valued/spinful irreps) |
| `src/assembly/` | `spacegroup`, `pointgroup`, `subperiodicgroup`, `mspacegroup`, `generate`, `generators` |
| `src/tables/` | Tabulated generator strings and rotation/translation tables for all group types |
| `src/collection_extensions.jl` | Additional methods on `Collection{T}` for specific `T` |
| `src/show.jl` | Pretty-printing for all custom types |
| `src/constants.jl` | `MAX_SGNUM` (230), `MAX_SUBGNUM`, `MAX_MSGNUM`, `MAX_MSUBGNUM`, `ENANTIOMORPHIC_PAIRS` |
| `src/utils.jl` | Internal utilities |
| `src/surface_unitcell.jl` | `surface_basis` (SmithNormalForm-based; TODO: move to Bravais) |
| `src/deprecations.jl` | Deprecated API aliases (e.g., `get_lgirreps` → `lgirreps`) |

> `src/special_representation_domain_kpoints.jl` exists but is currently **commented out** in `Crystalline.jl` due to a severe precompilation cost (~15 s for a `const`).

### Type hierarchy

```
AbstractOperation{D} <: AbstractMatrix{Float64}
  SymOperation{D}          — D×D rotation (SqSMatrix) + D-vector translation (SVector)
  MSymOperation{D}         — wraps SymOperation + time-reversal flag
  DSymOperation{D}         — wraps SymOperation + SU2 (a double group operation)

AbstractVec{D}
  KVec{D}                  — k₀ + kabc·(α,β,γ) in reciprocal coords (cnst + free·αβγ)
  RVec{D}                  — same structure, direct-space coords
  WyckoffPosition{D}       — wraps RVec + multiplicity + letter

AbstractGroup{D,O} <: AbstractVector{O}
  SpaceGroup{D}
  PointGroup{D}
  LittleGroup{D}           — also carries KVec and k-label
  SiteGroup{D}             — also carries WyckoffPosition and coset reps
  SubperiodicGroup{D,P}
  MSpaceGroup{D}
  GenericGroup{D}          — group from arbitrary operations; num = 0
  DSpaceGroup{D}, DPointGroup{D}, DLittleGroup{D}, DSiteGroup{D}
                           — double group counterparts, over DSymOperation{D}

AbstractIrrep{D}           — the `D`-prefixed types below are the double-valued (spinful)
                             irreps; they carry a `D`-prefixed group
  AbstractPGIrrep{D}       — PGIrrep{D}, DPGIrrep{D}: point group irreps
  AbstractLGIrrep{D}       — LGIrrep{D}, DLGIrrep{D}: little group irreps; have a
                             `translations` field for phase factors
  AbstractSiteIrrep{D}     — SiteIrrep{D}, DSiteIrrep{D}: site symmetry irreps; carry a
                             `pglabel` field

Collection{T} <: AbstractVector{T}    — thin wrapper around Vector{T}; same group for all T
CharacterTable{O} / ClassCharacterTable{O}  — characters vs operations/classes, over
                                              operations of type O
SymmetryVector{D,IR} <: AbstractSymmetryVector{D,IR} <: AbstractVector{Int}
BandRep{D,IR,SIR} <: AbstractSymmetryVector{D,IR}
CompositeBandRep{D,IR,SIR} <: AbstractSymmetryVector{D,IR}
                           — IR: little group irrep type; SIR: site symmetry irrep type
```

### Internal submodules
- **`SquareStaticMatrices`** (`src/SquareStaticMatrices.jl`): `SqSMatrix{D,T}` — stack-allocated D×D matrix; the rotation part of `SymOperation` is stored as one of these
- **`Jagged`** (`src/jaggedvector.jl`): `JaggedVector{T}` — memory-efficient vector-of-vectors; used in `SymmetryVector.multsv`
- **`SmithNormalForm`** (vendored in `.vendor/SmithNormalForm/`): Smith normal form; used in TQC analysis and `surface_basis`

### Subpackage: Bravais.jl
Located in `Bravais/`, a standalone registered package that Crystalline `@reexport`s in full. Provides: `DirectBasis`/`ReciprocalBasis`, `crystal`, `crystalsystem`, `bravaistype`, `centering`, `primitivebasismatrix`, `primitivize`/`conventionalize`/`transform`/`cartesianize`, Niggli reduction, and dual-basis utilities.

### Extensions (weak dependencies)
- **`CrystallineMakieExt`** (`ext/CrystallineMakieExt.jl`): Makie.jl plotting for Fourier lattice isosurfaces
- **`CrystallineGraphMakieExt`** (`ext/CrystallineGraphMakieExt.jl`): GraphMakie.jl visualization of group–subgroup graphs

### Data files (`data/`)

Do not attempt to read binary files (`.jld2`). Directory structure:

```
data/
  irreps/
    lgs/{1,2,3}d/irreps_data.jld2         — LGIrrep matrices, translations, realities
    lgs/{1,2,3}d/littlegroups_data.jld2   — LittleGroup operations and k-vectors
    pgs/3d/irreps_data.jld2               — PGIrrep data
  bandreps/                               — Bilbao's tabulated EBRs, as CSV; a lazy
                                            artifact, used only to validate `bandreps`
  wyckpos/                                — Wyckoff position data
  operations/                             — symmetry operation tables
  spacegroup_subgroups_data.jld2          — maximal subgroup / minimal supergroup graphs
  misc/
    ISOTROPY/                             — README.md has format docs; the *.txt data
                                              itself is a lazy artifact, not in the repo
    transformation_matrices_CDML2ITA.jl  — k-point transformation matrices between CDML and ITA settings
    CDML_RepresentationDomainSpecialKPoints_*.csv
```

All JLD2 files are opened once at package load (`__init__`) and kept open for the session via module-level `Ref`s (`LGIRREPS_JLDFILES`, `LGS_JLDFILES`, `PGIRREPS_JLDFILE`). This eliminates the per-call open/close overhead from `lgirreps`/`littlegroups`.

## Key Conventions and Domain Knowledge

### Coordinate conventions
All coordinates (rotation matrices, translation vectors, k-vectors, Wyckoff positions) are in **reduced/fractional** coordinates relative to the conventional unit cell, unless otherwise noted. `SymOperation{D}` stores (W, w): W is the D×D rotation and w the D-vector translation; the action on position **r** is W**r** + **w**. `matrix(op)` returns the D×(D+1) augmented matrix [W | w].

### LGIrrep phase convention (open issue #12)
`LGIrrep` stores matrices as loaded directly from ISOTROPY, without translation phase factors. When evaluated via `lgir(αβγ)`, the phase `exp(2πi k·τ)` is applied for each operation with translation factor τ (stored in `lgir.translations`). This follows the **Inui convention** and matches ISOTROPY (Acta Cryst. A69, 388, 2013):

```
D^k({R|τ}) = exp(+2πi k·τ) · D^k({R|0})
```

This is **opposite in sign** to Herring (1937), Kovalev, and arguably more physically natural (where one usually writes Bloch functions as exp(+ik·r) and symmetry operations act inversely). The sign discrepancy is documented in detail at `src/littlegroup_irreps.jl:139` and is tracked as **issue #12**.

### `littlegroups` vs. `spacegroup`
`littlegroups` returns little groups **without** centering-copy operations (e.g., for a body-centered group, `{1|½,½,½}` is excluded). `spacegroup` includes all centering copies. This difference is intentional and important when working with both.

### k-vectors and free parameters
`KVec{D}` represents **k** = k₀ + kabc·(α,β,γ). For special (high-symmetry) k-points, `free(kv)` is zero. For lines/planes, it is nonzero; evaluate at a specific point with `kv([α,β,γ])` or `kv(α, β, γ)`. The same functor interface applies to `LGIrrep`.

### Labels
Irrep labels follow **CDML** conventions. The k-point label is the leading alphabetic part of the irrep label; `klabel("X₃")` → `"X"`. Group labels use **IUC** (Hermann-Mauguin) notation via `iuc(...)`, or Schoenflies via `schoenflies(...)`.

### Reality
`Reality` is `@enum Reality::Int8` with values `REAL=1`, `PSEUDOREAL=-1`, `COMPLEX=0`, `UNDEF=2`. Characterizes irreps under time-reversal (Herring criterion). Use `realify` to form physically real irreps (co-representations).

### Dimensionality
Most APIs accept dimension `D` either as `Val{D}()` (preferred internally) or as a plain `Integer` (convenience). Supported: 1D, 2D, 3D for space groups, little groups, and irreps; 3D only for point groups; subperiodic groups cover layer (D=3, P=2), rod (D=3, P=1), and frieze (D=2, P=1) groups.

## Dependencies

**Core:** `LinearAlgebra`, `StaticArrays`, `JLD2`, `PrettyTables`, `DelimitedFiles`, `Reexport`, `DocStringExtensions`, `Graphs`, `Combinatorics`, `PrecompileTools`, `Statistics`, `SparseArrays`, `Contour`, `Meshing`, `LayeredLayouts`

**Optional (weak deps):** `Makie`, `GraphMakie`

**Subpackages/vendored:** `Bravais` (local, registered separately), `SmithNormalForm` (vendored at `.vendor/`)

## Testing

47 test files in `test/`, all run via `test/runtests.jl`. Selected highlights:

| Test file(s) | Coverage |
|---|---|
| `symops.jl`, `groups_xyzt_vs_coded.jl`, `generators_xyzt_vs_coded.jl` | Symmetry operations, group assembly from generators |
| `parsed_vs_loaded_littlegroup_irreps.jl` | JLD2-loaded irreps vs. freshly parsed ISOTROPY data |
| `irreps_orthogonality.jl`, `multtable.jl`, `chartable.jl` | Great orthogonality theorem, multiplication tables, character tables |
| `irreps_reality.jl`, `irreps_physical_reality.jl` | Herring criterion and co-rep construction; `physical_realify` and its time-reversal convention |
| `lgirreps_vs_pgirreps_at_Gamma.jl` | LGIrreps at Γ must match PGIrreps |
| `compatibility.jl` | Compatibility relations between k-points |
| `bandreps.jl`, `bilbao_bandreps.jl`, `classification.jl` | EBRs (computed, and vs. Bilbao's tables) and topological classification |
| `wyckoff.jl`, `isomorphic_parent_pointgroup.jl` | Wyckoff positions, site groups, isomorphic point groups |
| `primitivize_irreps.jl` | Irrep transformation to primitive basis |
| `grouprelations.jl` | Sub-/supergroup data integrity |
| `mspacegroup.jl` | Magnetic space groups |
| `double/` | Double groups and double-valued (spinful) irreps |
| `symeigs_analysis.jl` | Symmetry eigenvalue → irrep assignment |

## CI/CD

- GitHub Actions: tests on Julia 1.12 and latest, across multiple OS
- Documentation: Documenter.jl, auto-deployed
- Coverage: Codecov
