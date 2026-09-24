# Changelog

## v0.7.0 (unreleased)

### Breaking changes, and how to update

- **The k-points of `calc_bandreps` are sorted deterministically**, and hence so are those
  of the `SymmetryVector`s derived from them: by decreasing little group order, with ties
  broken alphabetically by k-label (Greek letters first; e.g., `[Γ, R, M, X]` for space
  group 221). Previously, the sorting followed the iteration order of the `Dict` returned by
  `lgirreps`, and so could change between Julia versions (as it did in Julia 1.13). Code
  that assumes a specific k-point sorting — e.g., by indexing into `irreps(brs)`, or via
  string comparisons against printed `SymmetryVector`s — may need updating.
