# ISOTROPY ISO-IR dataset
These file `CIR_data.txt` and `PIR_data.txt` are sourced from the
[ISOTROPY software's ISO-IR dataset](https://stokes.byu.edu/iso/irtables.php)
and contain data needed to generate *space group* irreps.
In ISOTROPY, the data files are extracted in Fortran using associated files
`CIR_data.f` and `PIR_data.f` (which we do not include here).
We extract the data files using Julia instead (see `build/parse_isotropy.jl`).

Two datasets are included in ISOTROPY, one for "ordinary" irreps (`CIR_data.txt`)
and one for "physically real" irreps/coreps in a real form (`PIR_data.txt`).
In practice, we only use the `CIR_data.txt` dataset in Crystalline.
(The `PIR_data.txt` dataset can be used via `build/parse_isotropy.jl` to obtain
physically real _space group_ (i.e. not little group) irreps, however.)

## References
Included below is the original header of the `CIR_data.txt` and `PIR_data.txt`
files (omitted in the included files for ease of parsing):

> ISO-IR: Complex Irreducible Representations of the 230 Crystallographic Space Groups
> 2011 Version
> Harold T. Stokes and Branton J. Campbell, 2013

The corresponding journal reference is:
- H. T. Stokes, B. J. Campbell, and R. Cordes, "Tabulation of Irreducible Representations of the Crystallographic Space Groups and Their Superspace Extensions", [Acta Cryst. A. **69**, 388-395 (2013)](https://doi.org/10.1107/S0108767313007538).