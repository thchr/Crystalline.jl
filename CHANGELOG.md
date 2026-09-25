# Changelog

## v0.7.0 (unreleased)

Adds spinful (double-group) irreps. Breaking, hence the minor version bump.

### Breaking changes, and how to update

- **`bandreps` now returns the band representations that Crystalline computes itself** —
  what `calc_bandreps` returned (`calc_bandreps` is deprecated to it) — rather than the
  Bilbao Crystallographic Server's tabulated EBRs.
  - This is the one change here that does **not** announce itself: `bandreps(sgnum, D; …)`
    still works, but returns a `Collection{<:BandRep}` in place of a `BandRepSet`.
  - The returned set may be larger: `bandreps` can include non-elementary "exceptional"
    band representations. This makes no difference for band connectivity or topology.
  - Code that reaches for `BandRepSet` fields, or that relies on the set being strictly
    elementary, needs review.
- **The `BandRep` and `BandRepSet` types that held Bilbao's tables were removed**, together
  with their parser; `NewBandRep` takes over the name `BandRep`. The `data/bandreps/` tables
  are now a lazy artifact, used only to validate `bandreps` in Crystalline's test suite.
- **`basisdim` no longer accepts a `BandRepSet`**; it takes a collection of band
  representations, an integer matrix, or a `Smith` factorization, and now lives alongside
  `indicator_group`.
- **`matching_littlegroups` and `matching_lgirreps` were removed** (unexported and unused),
  as were the deprecations `wyck(::BandRep)` and `matrix(::BandRepSet)`. `wyckbasis` is
  superseded by `Crystalline.smith_column_bases` (private API).
- **Long-standing deprecations (2021–2024) were removed**:
  - `kvec`, `wyck`, `kstar` → `position`, `parent`, `orbit`
  - `WyckPos` → `WyckoffPosition`
  - `get_littlegroups`, `get_lgirreps`, `get_pgirreps`, `get_wycks` → `littlegroups`,
    `lgirreps`, `pgirreps`, `wyckoffs`
  - `CharacterTable(::AbstractVector{<:AbstractIrrep})` → `characters`
  - `SiteGroup(::SpaceGroup, ::WyckoffPosition)` → `sitegroup`
  - `IrrepCollection` → `Collection{<:AbstractIrrep}`

  The more recent `classification`, `nontrivial_factors` and `symeigs_analysis` are kept.
- **The k-points of `bandreps` are sorted deterministically**, and hence so are those
  of the `SymmetryVector`s derived from them: by decreasing little group order, with ties
  broken alphabetically by k-label (Greek letters first; e.g., `[Γ, R, M, X]` for space
  group 221). Previously, the sorting followed the iteration order of the `Dict` returned by
  `lgirreps`, and so could change between Julia versions (as it did in Julia 1.13). Code
  that assumes a specific k-point sorting — e.g., by indexing into `irreps(brs)`, or via
  string comparisons against printed `SymmetryVector`s — may need updating.
- **`SymmetryVector{D}` → `SymmetryVector{D, IR}`**, `BandRep{D}` →
  `BandRep{D, IR, SIR}`, `CompositeBandRep{D}` → `CompositeBandRep{D, IR, SIR}` (`IR`:
  little group irrep type, `SIR`: site irrep type), and `AbstractSymmetryVector{D}` →
  `AbstractSymmetryVector{D, IR}`.
  - Signatures: `Vector{SymmetryVector{D}}` → `Vector{<:SymmetryVector{D}}` (and likewise
    for `Collection{BandRep{D}}` etc.). Scalar `::SymmetryVector{D}` args are unaffected.
  - Struct fields: `::SymmetryVector{D}` is now abstract; make the struct parametric, or use
    `SymmetryVector{D, LGIrrep{D}}` for spinless-only code.
  - Construction: `SymmetryVector{D}(lgirsv, multsv, μ)` →
    `SymmetryVector(lgirsv, multsv, μ)` (same for `BandRep`, `CompositeBandRep`).
  - Serialized (JLD2) instances of these types do not load directly; convert them, e.g. by
    loading into a stand-in struct of the old layout and rebuilding.
- **`BandRep` has no `spinful` field**: use `isspinful(br)`, and construct with
  `BandRep(siteir, n, timereversal)`.
- **`CharacterTable{D}` → `CharacterTable{O}`**, and `ClassCharacterTable{D}` →
  `ClassCharacterTable{O}`, with `O` the operation type (e.g. `SymOperation{3}`).
  Construct with `CharacterTable(ops, irlabs, table[, tag])`.
- **`physical_realify` returns a different (still explicitly real) basis for some spinless
  irreps.** It now solves for the same intertwiner direction as in the spinful case — the
  unitary part `Γ` of time reversal, mapping `conj(D)` onto `D`, rather than its inverse —
  which flips the sign of the transform for 18 of the 42 point group coreps that are not
  already real. Only the basis changes: the matrices remain real, equivalent to the input,
  and of unchanged characters and reality. Code that hardcodes specific physically real
  matrices, or a specific basis of tight-binding terms derived from them, may need updating.
- **Julia 1.12 is now the minimum supported version** (was 1.10). Bravais.jl, which is
  versioned separately, continues to support 1.10.
- **Printing**: band representations and symmetry vectors print with a spin tag,
  e.g. `SymmetryVector{3} (spinless)`, and the spin label in
  `Collection{BandRep}` summaries is `spinless`/`spinful` rather than `spin-1`/`spin-½`.
  Character tables print as `CharacterTable for ⋕9 (4) (spinless):`, dropping the type
  parameter. Update string comparisons against printed output.

### New

- Double groups: `SU2`, `DSymOperation`, `DSpaceGroup`, `DLittleGroup`, `DPointGroup`,
  `DSiteGroup`, `isbarred`.
- A `spinful` keyword argument on `spacegroup`, `pointgroup`, `littlegroups`, `lgirreps`,
  `pgirreps`, `sitegroups` and `bandreps` selects the double group or its
  double-valued irreps. As for the dimension, `spinful = Val(true)` keeps the return type
  inferrable, while `spinful = true` is a convenience that does not.
- Double-valued irreps: `DLGIrrep`, `DPGIrrep`, via `lgirreps(sgnum, Val(3);
  spinful=Val(true))` and `pgirreps(iuclab, Val(3); spinful=Val(true))`. Point group and
  site symmetry irreps also have Mulliken labels (`mulliken=true`), as in Bilbao.
- Double-valued site symmetry irreps: `DSiteIrrep`, via `siteirreps` of a double site
  symmetry group: `sitegroup(sg::DSpaceGroup, wp)`, `sitegroups(sgnum; spinful=Val(true))`.
- Spinful band representations: `bandreps(sgnum, Val(3); spinful=Val(true))`, induced
  from the double-valued site symmetry irreps.
- `physical_realify` also accepts double-valued irreps. A real form does not generally
  exist for them, since time reversal squares to `-1`; instead, the returned matrices obey
  `J*conj(D)*J' = D` with `J = iσʸ ⊗ 𝟙ₙ`, and are real where that is possible
  (i.e., unless the irrep is pseudoreal).
- `timereversal_unitary` returns the unitary part `Γ` of time reversal `T = ΓK` in the basis
  that `physical_realify` establishes: the identity for spinless irreps, and `J` for spinful
  ones. Models built on these irreps must adopt the same convention.
- `doublegroup`: the double group of a space, little, point, or site symmetry group.
- SU(2) elements: `SU2(op, sgnum)` returns the SU(2) element of a spatial operation, and
  `SU2(dop)` the one carried by a double group operation. The tabulated elements follow
  Bilbao and Altmann & Herzig, whose Cartesian frame differs from Crystalline's for
  hexagonal and trigonal lattices (see the `SU2` docstring).
- `realify`, `calc_reality`, `characters`, `classes`, `subduction_count`, `primitivize`,
  `collect_compatible` and `collect_irrep_annotations` accept double-valued irreps.
- Abstract supertypes to dispatch on, all exported: `AbstractLGIrrep`, `AbstractPGIrrep`
  and `AbstractSiteIrrep` span an irrep kind and its double-valued counterpart, while
  `AbstractSpaceGroup`, `AbstractPointGroup`, `AbstractLittleGroup` and `AbstractSiteGroup`
  span a kind of group and its ordinary, double, magnetic and subperiodic variants. With
  the latter, `iuc`, `centering`, `issymmorph`, `reduce_ops` and `primitivize` extend to the
  double groups. (The older supertypes — `Crystalline.AbstractIrrep`, `AbstractGroup`,
  `AbstractVec`, `AbstractOperation`, `AbstractSymmetryVector` — stay unexported for now;
  they can still be imported by name.)
- `reduce_ops` (and so `primitivize`) works for any operation type, including
  `DSymOperation`, so a `DSpaceGroup` can now be primitivized.
- `iuc(::SubperiodicGroup)`, which returns its subperiodic label.
- `isspinful` for irreps, character tables, symmetry vectors, and band representations.
- Documentation: a new manual page, *Double groups & spinful irreps*, covering the double
  group types, the SU(2) conventions, double-valued irreps, and their behavior under time
  reversal.
