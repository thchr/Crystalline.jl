# Irreducible representations

Crystalline.jl provides easy access to crystallographic point group irreps, "small" irreps of space groups, and site symmetry group irreps.
Currently, we only provide access to spinless (or "single-valued") irreps.

## Data sources
Point group irreps are obtained from the Bilbao Crystallographic Server's [Representations PG program](https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=spg) and small irreps of space groups are obtained from [ISOTROPY's 2011 ISO-IR dataset](https://stokes.byu.edu/iso/irtables.php).

If these functionalities are used in published research, please cite the original publications (listed in associated function docstrings).

## Point group irreps
Irreps for the crystallographic point groups are accessible via [`pgirreps`](@ref), with the point group specified either by IUC label and dimensionality.
As an example, we may load the irreps of the 6mm (C₆ᵥ in Schoenflies notation; see also [`schoenflies(::PointGroup)`](@ref)) point group in 3D.
```@example pgirs
using Crystalline

pgirs = pgirreps("6mm", Val(3))
```
Frequently, the character table of the associated irreps is more informative than the irrep matrices themselves. We can construct this table using [`characters`](@ref), which returns a `CharacterTable`:
```@example pgirs
characters(pgirs)
```

The characters are functions of the conjugacy class (i.e., the characters of operations in the same conjugacy class are equal). Thus, a more compact representation of the character table can be achieved by a class-resolved table, achievable via [`classcharacters`](@ref):
```@example pgirs
classcharacters(pgirs)
```

### Notation 
The default point group irrep labeling follows the Bilbao Crystallographic Server's labeling, which in turn follows the 1963 labelling of Koster, Dimmock, Wheeler, & Statz [^2] (which is also followed e.g. by CDML [^1] labeling as well as Bradley and Cracknell's book).
Associated Muliken (or "spectroscopist's") notation can be obtained via `mulliken`.


## Small little group irreps
Little group irreps are accessible via [`lgirreps`](@ref) and provided with CDML [^1] labels (courtesy of ISOTROPY).
As an example, we can obtain the irreps of space group 183 (P6mm; the trivial 3D extension of plane group 17, which in turn is the space group extension of point group 6mm from above) by:
```@example lgirs
using Crystalline

lgirsd = lgirreps(183, Val(3))
```
which returns a dictionary of `LGIrrep`s, indexed by k-labels given as `String`s, corresponding to different little groups.
In general, we include all the little groups included in ISOTROPY; unfortunately, there is no strict guarantee that this includes a full listing of all inequivalent irreps (although it is typically true). The listing typically contains both special points, lines, and planes (and also always the general point).

We can inspect the small irreps of any particular **k**-point, accessing it via its canonical label.
As before, it is most instructive to inspects the associated character tables:
```@example lgirs
lgirs = lgirsd["Γ"] # small irreps at the Γ point
characters(lgirs)
```

## Space group irreps
We currently do not provide access to "full" space group irreps. They can, however, be readily built by induction from small irreps. Specifically, every small irrep $D_{\mathbf{k}}^\alpha$ associated with the little group $G_{\mathbf{k}}$, induces a space group irrep, sometimes denoted ${}^*D_{\mathbf{k}}^{\alpha}$ or $D^{\alpha}_{\mathbf{k}}\uparrow G$, in the full space group $G$:[^Inui]

```math
    [{}^*D_{\mathbf{k}}^{\alpha}(g)]_{ij}
    =
    \begin{cases}
    D_{\mathbf{k}}^{\alpha}(h_i^{-1}gh_j) & \text{if }h_i^{-1}gh_j \in G_{\mathbf{k}}\\
    \boldsymbol{0}_{d_{\mathbf{k}}^{\alpha}\times d_{\mathbf{k}}^{\alpha}} & \text{otherwise}
    \end{cases},
```

where $d_{\mathbf{k}}^{\alpha}$ is the dimension of the small irrep $D^{\alpha}_{\mathbf{k}}$, $\boldsymbol{0}_{d_{\mathbf{k}}^{\alpha}\times d_{\mathbf{k}}^{\alpha}}$ is a $d_{\mathbf{k}}^{\alpha}\times d_{\mathbf{k}}^{\alpha}$ zero matrix, and $h_i$ and $h_j$ iterate over the (left) coset representatives of $G_{\mathbf{k}}$ in $G$ (of which there are $|\mathrm{star}\{\mathbf{k}\}|$, i.e., the order of the star of $\mathbf{k}$). The induced irrep ${}^*D_{\mathbf{k}}^{\alpha}$ is consequently a $d_{\mathbf{k}}^{\alpha}|\mathrm{star}\{\mathbf{k}\}|\times d_{\mathbf{k}}^{\alpha}|\mathrm{star}\{\mathbf{k}\}|$ matrix.

[^Inui]: Inui, Tanabe, & Onodera, *Group Theory and its Applications in Physics*, Springer (1990). Section 11.9.

## Site symmetry irreps
To obtain irreps associated with a given site symmetry group (see [`SiteGroup`](@ref), use the [`find_isomorphic_parent_pointgroup`](@ref) to obtain the "parent" point group and any relevant permutation differences between the two groups. Associated irreps can then be inferred by a suitable permuation of results obtained from [`pgirreps`](@ref).

## Time-reversal symmetry and "physically real" irreps
Irreps returned in Crystalline.jl do not assume time-reversal symmetry by default. 
To incorporate time-reversal symmetry (or, equivalently, to obtain associated "physically real" irreps), which may cause irreps to "stick together", see [`realify`](@ref) (which takes a vector of `PGIrrep`s or `LGIrrep`s).

As an example, the Γ₃, Γ₄, Γ₅, and Γ₆ irreps of point group 6 (C₆) are intrinsically complex in the absence of time-reversal symmetry:
```@example realirs
using Crystalline

pgirs = pgirreps("6", Val(3))
characters(pgirs)
```
When time-reversal symmetry is incorporated, the irreps stick together pairwise and have real characters:
```@example realirs
pgirs′ = realify(pgirs)
characters(pgirs′)
```

To inspect the reality type of a given irrep, see [`reality`](@ref).
Possible types are `REAL`, `COMPLEX`, and `PSEUDOREAL` (the latter does not arise for point groups):
```@example realirs
label.(pgirs) .=> reality.(pgirs)
```
The reality type can be computed ab initio via `calc_reality`, using the Frobenius criterion for `PGIrrep`s and the Herring criterion for `LGIrrep`s.

## References
[^1]: Cracknell, A.P., Davies, B.L., Miller, S.C., & Love, W.F., *Kronecker Product Tables, Vol. 1. General Introduction and Tables of Irreducible Representations of Space Groups*, New York: IFI/Plenum (1979).

[^2]: Koster, G.F., Dimmock, J.O., Wheeler, R.G., & Statz, H., *Properties of the Thirty-two Point Groups*, Cambridge: MIT Press (1963).