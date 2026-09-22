# Changelog

## v0.7.0 (unreleased)

Adds spinful (double-group) irreps. Breaking, hence the minor version bump.

### Breaking changes, and how to update

- **The k-points of `calc_bandreps` are sorted deterministically**, and hence so are those
  of the `SymmetryVector`s derived from them: by decreasing little group order, with ties
  broken alphabetically by k-label (Greek letters first; e.g., `[Γ, R, M, X]` for space
  group 221). Previously, the sorting followed the iteration order of the `Dict` returned by
  `lgirreps`, and so could change between Julia versions (as it did in Julia 1.13). Code
  that assumes a specific k-point sorting — e.g., by indexing into `irreps(brs)`, or via
  string comparisons against printed `SymmetryVector`s — may need updating.
- **`SymmetryVector{D}` → `SymmetryVector{D, IR}`**, `NewBandRep{D}` → `NewBandRep{D, IR, SIR}`,
  `CompositeBandRep{D}` → `CompositeBandRep{D, IR, SIR}` (`IR`: little group irrep type, `SIR`:
  site irrep type), and `AbstractSymmetryVector{D}` → `AbstractSymmetryVector{D, IR}`.
  - Signatures: `Vector{SymmetryVector{D}}` → `Vector{<:SymmetryVector{D}}` (and likewise for
    `Collection{NewBandRep{D}}` etc.). Scalar arguments `::SymmetryVector{D}` are unaffected.
  - Struct fields: `::SymmetryVector{D}` is now abstract; make the struct parametric, or use
    `SymmetryVector{D, LGIrrep{D}}` for spinless-only code.
  - Construction: `SymmetryVector{D}(lgirsv, multsv, μ)` → `SymmetryVector(lgirsv, multsv, μ)`
    (same for `NewBandRep`, `CompositeBandRep`).
  - Serialized (JLD2) instances of these types do not load directly; convert them, e.g. by
    loading into a stand-in struct of the old layout and rebuilding.
- **`NewBandRep` has no `spinful` field**: use `isspinful(br)`, and construct with
  `NewBandRep(siteir, n, timereversal)`.
- **`CharacterTable{D}` → `CharacterTable{O}`**, and `ClassCharacterTable{D}` →
  `ClassCharacterTable{O}`, with `O` the operation type (e.g. `SymOperation{3}`).
  Construct with `CharacterTable(ops, irlabs, table[, tag])`.
- **Positional argument order**: `littlegroups(sgnum, Val(D), jldfile)` and
  `lgirreps(sgnum, Val(D), lgs_jldfile, irs_jldfile)` now take a spinful argument
  (`Val(false)`) before the JLD2 files. Calls without explicit files are unaffected.
- **Printing**: band representations and symmetry vectors print with a spin tag,
  e.g. `SymmetryVector{3} (spinless)`, and the spin label in `BandRepSet` and
  `Collection{NewBandRep}` summaries is `spinless`/`spinful` rather than `spin-1`/`spin-½`.
  Update string comparisons against printed output.

### New

- Double groups: `SU2`, `DSymOperation`, `DSpaceGroup`, `DLittleGroup`, `DPointGroup`,
  `DSiteGroup`, `su2`, `isbarred`; `spacegroup`, `pointgroup` accept a spinful argument.
- Double-valued irreps: `DLGIrrep`, `DPGIrrep`, via `lgirreps(sgnum, Val(3), Val(true))` and
  `pgirreps(iuclab, Val(3), Val(true))` (or `lgirreps(sgnum, 3, true)`, etc.).
- Abstract supertypes `AbstractLGIrrep` and `AbstractPGIrrep`; `realify`, `calc_reality`,
  `characters`, `classes`, `subduction_count`, `primitivize`, `collect_compatible`,
  `collect_irrep_annotations` accept double-valued irreps.
- `isspinful` for irreps, symmetry vectors, and band representations.
