# Publishing data sets

Some of Crystalline's data is not kept in this repository but published as assets of a
`data-v*` [GitHub release](https://github.com/thchr/Crystalline.jl/releases) and consumed as
Julia artifacts, declared in `Artifacts.toml` at the repository root. This file is the
procedure for changing that data; `build/data_release.jl` does the mechanical parts.

## What is published, and why

A data set belongs in a release rather than in the repository when users do not need it —
when it exists to *regenerate* or *validate* what Crystalline ships, rather than to run it.
Such data would otherwise be downloaded by everyone, on every release, and read by no one.

| artifact | contents | used by |
|---|---|---|
| `isotropy` | ISOTROPY's `CIR_data.txt`, `PIR_data.txt` (55 MB → 1.8 MB) | `build/write_littlegroup_irreps.jl`, `test/parsed_vs_loaded_littlegroup_irreps.jl` |
| `bilbao_spinless_irreps` | Bilbao's single-valued little group irreps (16 MB → 1.8 MB) | `test/bilbao_vs_isotropy.jl` |

Every artifact is `lazy`: nothing is downloaded when Crystalline is installed, only on first
use, and the result is cached in the depot (`~/.julia/artifacts`) under its tree hash. An
unchanged data set is therefore never downloaded twice, no matter how many Crystalline
versions refer to it — which is the point of the exercise.

Each artifact unpacks to **its files directly**, under their own names: an artifact is
already its own namespace, so rebuilding the `data/` directories a file happens to sit in
would only add structure that nothing reads. A data set whose files do need a hierarchy of
their own — several dimensions or space groups, say — declares it by giving
`source => path-in-artifact` pairs in `DATASETS`. The eventual `data` artifact is the one
case where the hierarchy is the data set's own, since it *is* `data/`.

## Amending a data set

**1. Get the current data.** If it is still in `data/` — because you have just rebuilt it
with one of the `build/write_*.jl` scripts — there is nothing to do. If it is not (the usual
case, since these files are gitignored), fetch what is currently published:

```julia
julia --project=build
julia> include("build/data_release.jl")
julia> dir = stage_from_artifact("isotropy")   # downloads if needed; returns a temp dir
```

**2. Change it** — edit, replace or add files in that directory, or rebuild the data into
`data/`. If you are adding a file, list it in the `DATASETS` table at the top of
`build/data_release.jl` first.

**3. Package it and read off the hashes.**

```julia
julia> tarball, tree_hash, sha256sum = package("isotropy", dir)
```

or, for everything at once, from the shell:

```bash
julia --project=build build/data_release.jl data-v0.0.2
```

which stages, packages, prints the `Artifacts.toml` entries and the `gh` command to run.
Tarballs land in `build/data-release/` (gitignored). The two hashes are of different things
and both are needed: `git-tree-sha1` is the git hash of the *unpacked tree*, which Pkg
verifies after unpacking, and `sha256` is of the *tarball bytes*, verified on download.

**4. Upload, to a new tag.**

```bash
gh release create data-v0.0.2 --title "..." --notes "..."
gh release upload data-v0.0.2 build/data-release/isotropy.tar.gz
```

Only the data sets that actually changed need to be uploaded; the `Artifacts.toml` entries
of the others keep pointing at the tag they were published under. Assets may also be
replaced within a tag (`gh release upload --clobber`) or removed
(`gh release delete-asset`), but see the rule below before doing so.

**5. Update `Artifacts.toml`** with the printed entries, and commit it together with
whatever change prompted the new data.

**6. Check it from a clean slate**, since a stale copy in your depot will hide a wrong hash:

```bash
rm -rf ~/.julia/artifacts/<old tree hash>
julia --project=test -e 'using Pkg.Artifacts: ensure_artifact_installed
                         ensure_artifact_installed("isotropy", "Artifacts.toml")'
```

A wrong `sha256` fails at download, a wrong `git-tree-sha1` after unpacking; both are loud.

## The one rule that matters

**Never delete or overwrite an asset that a released `Artifacts.toml` points at.** Those
URLs are baked into versions of Crystalline that are already in the registry, and a user
installing such a version will fetch them. Always publish changed data under a *new* tag.
(Before the first release referring to a tag, it is of course still free to move.)

## Adding a new data set

1. Add it to `DATASETS` in `build/data_release.jl`, with its paths relative to `data/`.
2. Add those paths to `.gitignore`, so the data does not also end up committed.
3. Package and upload as above, and add the printed stanza to `Artifacts.toml`.
4. Consume it. From `src/`, `artifact"name"` works directly. From `test/` or `build/` it
   does **not**: `@artifact_str` searches upwards from the calling file for an
   `Artifacts.toml` and gives up at the first `Project.toml` it meets, which in both cases
   is their own. Name the file explicitly there instead:

   ```julia
   using Pkg.Artifacts: ensure_artifact_installed
   dir = ensure_artifact_installed("name", joinpath(pkgdir(Crystalline), "Artifacts.toml"))
   ```

5. In a test, wrap that call in a `try`/`catch` and skip on failure, so that a machine
   without a network connection does not fail the suite (see `test/bilbao_vs_isotropy.jl`).

## Notes

- **Packaging is reproducible.** `Tar.create` normalizes timestamps, ownership and
  permissions and `gzip -n` writes no timestamp of its own, so the same input gives a
  byte-identical tarball and the same two hashes. Re-packaging unchanged data is therefore
  safe and produces no spurious diff.
- **Stage on a POSIX filesystem.** A git tree hash records the executable bit, and a Windows
  drive mounted under WSL reports every file as executable while silently ignoring `chmod`,
  which would give a different tree hash for identical data. `stage` uses a temporary
  directory for this reason, and checks the result.
- **Versioning.** `data-vX.Y.Z` tags are independent of Crystalline's own versions; bump the
  patch for a repackaging, the minor for new or extended data, the major for data that
  changes meaning under an unchanged name.
- **Removing something from `data/` does not shrink `.git`**: the history still holds it.
  The saving is in what users download from the registry, which is what matters here.
