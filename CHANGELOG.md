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
- **`SymmetryVector{D}` → `SymmetryVector{D, IR}`**, `NewBandRep{D}` →
  `NewBandRep{D, IR, SIR}`, `CompositeBandRep{D}` → `CompositeBandRep{D, IR, SIR}` (`IR`:
  little group irrep type, `SIR`: site irrep type), and `AbstractSymmetryVector{D}` →
  `AbstractSymmetryVector{D, IR}`.
  - Signatures: `Vector{SymmetryVector{D}}` → `Vector{<:SymmetryVector{D}}` (and likewise
    for `Collection{NewBandRep{D}}` etc.). Scalar `::SymmetryVector{D}` args are unaffected.
  - Struct fields: `::SymmetryVector{D}` is now abstract; make the struct parametric, or use
    `SymmetryVector{D, LGIrrep{D}}` for spinless-only code.
  - Construction: `SymmetryVector{D}(lgirsv, multsv, μ)` →
    `SymmetryVector(lgirsv, multsv, μ)` (same for `NewBandRep`, `CompositeBandRep`).
  - Serialized (JLD2) instances of these types do not load directly; convert them, e.g. by
    loading into a stand-in struct of the old layout and rebuilding.
- **`NewBandRep` has no `spinful` field**: use `isspinful(br)`, and construct with
  `NewBandRep(siteir, n, timereversal)`.
- **`CharacterTable{D}` → `CharacterTable{O}`**, and `ClassCharacterTable{D}` →
  `ClassCharacterTable{O}`, with `O` the operation type (e.g. `SymOperation{3}`).
  Construct with `CharacterTable(ops, irlabs, table[, tag])`.
- **Printing**: band representations and symmetry vectors print with a spin tag,
  e.g. `SymmetryVector{3} (spinless)`, and the spin label in `BandRepSet` and
  `Collection{NewBandRep}` summaries is `spinless`/`spinful` rather than `spin-1`/`spin-½`.
  Character tables print as `CharacterTable for ⋕9 (4) (spinless):`, dropping the type
  parameter. Update string comparisons against printed output.

### New

- Double groups: `SU2`, `DSymOperation`, `DSpaceGroup`, `DLittleGroup`, `DPointGroup`,
  `DSiteGroup`, `su2`, `isbarred`.
- A `spinful` keyword argument on `spacegroup`, `pointgroup`, `littlegroups`, `lgirreps`,
  `pgirreps`, `sitegroups` and `calc_bandreps` selects the double group or its
  double-valued irreps. As for the dimension, `spinful = Val(true)` keeps the return type
  inferrable, while `spinful = true` is a convenience that does not.
- Double-valued irreps: `DLGIrrep`, `DPGIrrep`, via `lgirreps(sgnum, Val(3);
  spinful=Val(true))` and `pgirreps(iuclab, Val(3); spinful=Val(true))`. Point group and
  site symmetry irreps also have Mulliken labels (`mulliken=true`), as in Bilbao.
- Double-valued site symmetry irreps: `DSiteIrrep`, via `siteirreps` of a double site
  symmetry group: `sitegroup(sg::DSpaceGroup, wp)`, `sitegroups(sgnum; spinful=Val(true))`.
- Spinful band representations: `calc_bandreps(sgnum, Val(3); spinful=Val(true))`, induced
  from the double-valued site symmetry irreps.
- `physical_realify` also accepts double-valued irreps. A real form does not generally
  exist for them, since time reversal squares to `-1`; instead, the returned matrices obey
  `J*conj(D)*J' = D` with `J = iσʸ ⊗ 𝟙ₙ`, and are real where that is possible
  (i.e., unless the irrep is pseudoreal).
- `timereversal_unitary` returns the unitary part `Γ` of time reversal `T = ΓK` in the basis
  that `physical_realify` establishes: the identity for spinless irreps, and `J` for spinful
  ones. Models built on these irreps must adopt the same convention.
- **`physical_realify` returns a different (still explicitly real) basis for some spinless
  irreps.** It now solves for the same intertwiner direction as in the spinful case — the
  unitary part `Γ` of time reversal, mapping `conj(D)` onto `D`, rather than its inverse —
  which flips the sign of the transform for 18 of the 42 point group coreps that are not
  already real. Only the basis changes: the matrices remain real, equivalent to the input,
  and of unchanged characters and reality. Code that hardcodes specific physically real
  matrices, or a specific basis of tight-binding terms derived from them, may need updating.
- `doublegroup`: the double group of a space, little, point, or site symmetry group.
- SU(2) elements (`su2`) follow Bilbao and Altmann & Herzig, whose Cartesian frame differs
  from Crystalline's for hexagonal and trigonal lattices (see the `su2` docstring).
- Abstract supertypes `AbstractLGIrrep`, `AbstractPGIrrep`, and `AbstractSiteIrrep`;
  `realify`, `calc_reality`, `characters`, `classes`, `subduction_count`, `primitivize`,
  `collect_compatible`, `collect_irrep_annotations` accept double-valued irreps.
- Internal abstract supertypes `Crystalline.AbstractSpaceGroup`, `AbstractPointGroup`,
  `AbstractLittleGroup` and `AbstractSiteGroup`, one per kind of group and spanning its
  ordinary, double, magnetic and subperiodic variants. With them, `iuc`, `centering`,
  `issymmorph`, `reduce_ops` and `primitivize` extend to the double groups.
- `reduce_ops` (and so `primitivize`) works for any operation type, including
  `DSymOperation`, so a `DSpaceGroup` can now be primitivized.
- `iuc(::SubperiodicGroup)`, which returns its subperiodic label.
- `isspinful` for irreps, character tables, symmetry vectors, and band representations.
