# Double groups & spinful irreps

```@meta
CurrentModule = Crystalline
```

Spin-½ degrees of freedom do not transform under a crystallographic group itself, but under its _double group_.
The reason is that a rotation by $2\pi$ acts as $-1$ on a spinor rather than as the identity: each spatial operation $g$ is lifted to two distinct operations, conventionally written $g$ and $\bar{E}g$, which differ by that $2\pi$ rotation $\bar{E}$.
The double group of a group $G$ consequently has $2|G|$ elements.
Its irreps fall in two families, according to how they represent $\bar{E}$.
Those representing it by $+1$ cannot distinguish $g$ from $\bar{E}g$, and are just the irreps of $G$ itself: the _single-valued_ (spinless) irreps.
Those representing it by $-1$ are the _double-valued_ (spinful) irreps, and are the ones that describe spin-½ degrees of freedom.

Crystalline.jl provides an interface to [Bilbao Crystallographic Server](https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg)'s tables of the double groups of space, point, and little groups and their double-valued irreps, as well as constructors to build the associated site symmetry groups, irreps, and the associated band representations.
Throughout, the word "spinful" refers to the double-valued (odd-integer multiples of spin-½) case and "spinless" to the single-valued (integer spin) one; [`isspinful`](@ref) distinguishes them.

!!! note
    Double groups are currently implemented in 3D only.

## Double group operations
A double group operation is a [`DSymOperation{D}`](@ref): a spatial [`SymOperation{D}`](@ref) together with the [`SU2`](@ref) element by which it acts on spin-½ degrees of freedom.
The SU(2) element of a spatial operation is obtained from the [`SU2`](@ref) constructor, which additionally requires the crystal system (see below):
```@example doublegroups
using Crystalline

SU2(S"-y,x-y,z", 183) # the SU(2) element of 3₀₀₁⁺ in space group 183 (P6mm)
```

The two operations that share a spatial part carry the SU(2) elements $U$ and $-U$, differing by a further $2\pi$ rotation; the one whose rotation angle stays within half a turn is the unbarred one, and [`isbarred`](@ref) tells the two apart, with [`seitz`](@ref) marking the barred one by a superscripted `ᵈ`.
Because composition of SU(2) elements is what keeps track of the $2\pi$ rotation, a double group multiplication table looks quite different from its spinless counterpart — e.g., a two-fold rotation no longer squares to the identity:
```@example doublegroups
dpg = pointgroup("2", Val(3); spinful=Val(true))
MultTable(dpg)
```

!!! details "SU(2) conventions: the Cartesian frame, and changes of basis or setting"
    The SU(2) element is not fixed by the rotation part alone: a rotation matrix in fractional coordinates does not determine the Cartesian rotation axis.
    [`SU2`](@ref) therefore also takes the crystal system, and returns elements tabulated by Altmann & Herzig[^Altmann], as used by the Bilbao Crystallographic Server.

    Those tables fix the orientation of the conventional basis $(\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3)$ relative to a Cartesian frame.
    For all but hexagonal and trigonal lattices, this agrees with Crystalline's own convention (as used by e.g. [`crystal`](@ref)), i.e., $\mathbf{a}_1 \parallel x$ with $\mathbf{a}_2$ in the $xy$-plane; for hexagonal and trigonal lattices, the tables instead use $x \parallel \mathbf{a}_1 + 2\mathbf{a}_2$, $y \parallel \mathbf{a}_1$, and $z \parallel -\mathbf{a}_3$.
    We keep the tabulated elements unchanged, so that they agree with Bilbao's.
    The frame has no bearing on irreps, characters, or band representations; it matters only if SU(2) elements are combined with quantities referred to Cartesian axes, such as spin operators taken as physical directions, orbitals transforming under Cartesian rotations, or external fields.

    A change of lattice basis does not change the physical operation, so [`transform`](@ref), [`primitivize`](@ref), and [`conventionalize`](@ref) carry the SU(2) element over unchanged.
    The same holds when `transform` is used for a change of _setting_ (as in [`conjugacy_relations`](@ref)) — but there the new basis is generally _not_ in the orientation assumed above, so the transformed elements need not equal those tabulated for the new setting (e.g., $2_{001}$ rewritten as $2_{010}$ retains $-\mathrm{i}\sigma_z$, where the table gives $-\mathrm{i}\sigma_y$ for $2_{010}$).
    To obtain the double group operations _of the new setting_, transform the spatial operations and attach the elements tabulated for it, i.e., `SU2(op′, sgnum′)`.

## Double groups
The double group of a space, point, little, or site symmetry group can be constructed from a single-valued counterpart via [`doublegroup`](@ref):
```@example doublegroups
doublegroup(pointgroup("6mm", Val(3)))
```
The returned group contains all $2|G|$ operations, with barred operations listed after their unbarred partners.
Including the barred operations explicitly means that the double group is a group in its own right, so that multiplication tables, conjugacy classes, and the Herring criterion require no special casing.

Equivalently, and usually more conveniently, the group constructors take a `spinful` keyword argument:
```@example doublegroups
sg = spacegroup(183, Val(3); spinful=Val(true))  # a `DSpaceGroup{3}`
nothing # hide
```
and likewise for [`pointgroup`](@ref), [`littlegroups`](@ref), and [`sitegroups`](@ref).
As with the dimension, passing `spinful` as a `Val` keeps the return type inferrable, whereas a plain `Bool` does not; `spinful = Val(false)` is the default throughout.

## Double-valued irreps
Double-valued irreps are obtained from the same accessors as the single-valued ones — [`pgirreps`](@ref), [`lgirreps`](@ref), and [`siteirreps`](@ref) — and are returned as [`DPGIrrep`](@ref), [`DLGIrrep`](@ref), and [`DSiteIrrep`](@ref), respectively.
The point group and little group accessors take the same `spinful` keyword argument as their group counterparts:
```@example doublegroups
dpgirs = pgirreps("321", Val(3); spinful=Val(true))
classcharacters(dpgirs)
```
Site symmetry irreps are instead obtained by passing a `DSiteGroup` to `siteirreps`:
```@example doublegroups
dsiteg = sitegroups(183, Val(3); spinful=Val(true))[5] # site group of Wyckoff position 2b
classcharacters(siteirreps(dsiteg))
```

The irrep matrices are defined over all $2|G|$ operations, with $D(\bar{E}g) = -D(g)$, so that the barred half of the character table is the negative of the unbarred half — as evident also in the table above.
Labels follow CDML for little group irreps and the Bilbao/Koster convention for point group irreps, in both cases with an appended `ˢ` marking the irrep as double-valued (where Bilbao writes an overline); [`mulliken`](@ref) gives the associated Mulliken labels:
```@example doublegroups
label.(dpgirs) .=> mulliken.(dpgirs)
```

## Time-reversal symmetry
Time reversal is incorporated exactly as in the spinless case, via [`realify`](@ref) (see also the [Irreps](irreps.md) page).
It applies wherever time reversal maps the group back to itself: always for point groups and site symmetry groups, and, for little groups, at a time-reversal invariant momentum (TRIM), i.e., where $\mathbf{k} \equiv -\mathbf{k}$ up to a reciprocal lattice vector.
Its _effect_, however, is qualitatively different from the spinless case, because time reversal then squares to $-1$ rather than to $+1$.
In terms of the reality type (see [`reality`](@ref) and [`calc_reality`](@ref)), the roles of `REAL` and `PSEUDOREAL` are interchanged:

| [`Reality`](@ref) | spinless    | spinful     |
|:------------------|:------------|:------------|
| `REAL`            | unchanged   | **doubled** |
| `PSEUDOREAL`      | **doubled** | unchanged   |
| `COMPLEX`         | glued to its complex conjugate partner | glued to its complex conjugate partner |

The physical consequence is Kramers degeneracy, and it is common to all three cases: a spinful corep is always even-dimensional — by doubling if `REAL`, by gluing to a partner if `COMPLEX`, and, if `PSEUDOREAL`, because a pseudoreal irrep is already even-dimensional on its own[^1].
The `REAL` case is the starkest: the double group of point group 1 has a single double-valued irrep, one-dimensional, which time reversal nonetheless sticks to a copy of itself:
```@example doublegroups
pgirs¹ = pgirreps("1", Val(3); spinful=Val(true))
label.(pgirs¹) .=> reality.(pgirs¹)
```
```@example doublegroups
pgirs¹′ = realify(pgirs¹)
label.(pgirs¹′) .=> irdim.(pgirs¹′) # a two-fold degeneracy, i.e., a Kramers pair
```
Conversely, a `PSEUDOREAL` double-valued irrep is left untouched by `realify`: it is already its own time-reversal partner.
This is the common case for double-valued irreps — unlike for single-valued point group irreps, where the pseudoreal type does not occur at all.

Away from a TRIM there is no Kramers degeneracy, and odd-dimensional spinful coreps do occur.
E.g., in space group 183 (P6mm), the coreps at Γ are all two-dimensional, while those at K — which is not a TRIM, since $-\mathbf{k}_{\text{K}}$ differs from $\mathbf{k}_{\text{K}}$ by no reciprocal lattice vector — are not:
```@example doublegroups
lgirsd = lgirreps(183, Val(3); spinful=Val(true))
irdim.(realify(lgirsd["Γ"])), irdim.(realify(lgirsd["K"]))
```

### Physically real form
As in the spinless case, [`physical_realify`](@ref) brings the matrices of a time-reversal-invariant irrep to a canonical form.
The convention is stated in terms of the unitary part $\Gamma$ of time reversal, $T = \Gamma K$ with $K$ complex conjugation: the returned matrices obey $\Gamma D^*(g) \Gamma^\dagger = D(g)$, with $\Gamma$ given by [`timereversal_unitary`](@ref).
For spinless irreps $\Gamma$ is chosen as the identity, so this is just the statement that the matrices are real.
For spinful irreps of dimension $2n$ that choice is impossible, since $T^2 = -1$ requires $\Gamma\Gamma^* = -\mathbf{1}$; the canonical choice is instead $\Gamma = \mathrm{i}\sigma_y \otimes \mathbf{1}_n$[^2], requiring the matrices take the block form
```math
    D = \begin{pmatrix} A & B \\ -B^* & A^* \end{pmatrix},
```
the $n$-block analogue of the SU(2) form $[a\ b;\ -b^*\ a^*]$[^3].
The matrices are additionally real whenever that is possible, i.e., unless the irrep is pseudoreal:
```@example doublegroups
ir = realify(pgirreps("3", Val(3); spinful=Val(true)))[2] # a COMPLEX corep: a real form exists
physical_realify(ir)
```
```@example doublegroups
ir′ = realify(pgirreps("222", Val(3); spinful=Val(true)))[1] # PSEUDOREAL: no real form
pir′ = physical_realify(ir′)
```
The defining relation holds in either case — here, for the pseudoreal one:
```@example doublegroups
Γ′ = timereversal_unitary(ir′)
all(D -> Γ′ * conj(D) * Γ′' ≈ D, pir′(nothing))
```

Any subsequent use of the irrep matrices alongside time reversal — e.g., in constructing symmetry-constrained tight-binding models — must adopt the same convention for $\Gamma$ (as returned by [`timereversal_unitary`](@ref)).

## Band representations
[`calc_bandreps`](@ref) takes the same `spinful` keyword argument, and returns the band representations induced from the double-valued site symmetry irreps:
```@example doublegroups
calc_bandreps(183, Val(3); spinful=Val(true))
```
In the presence of time-reversal symmetry, every such band representation has even filling since Kramers' theorem pairs the states of a spinful system.

## Data sources
The provided double-valued little group and point group irreps are obtained from the Bilbao Crystallographic Server's [Representations DSG](https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg) and Representations DPG programs; ISOTROPY, the source of Crystalline's spinless little group irreps, provides no double-valued data.
If used in research, please cite the original reference[^Elcoro] for [Representations DSG](https://cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=dbg).
The SU(2) elements are those of Altmann & Herzig[^Altmann], as also used by Bilbao.

[^Elcoro]: Elcoro et al. *Double crystallographic groups and their representations on the Bilbao Crystallographic Server*, [J. Appl. Cryst. **50**, 1457 (2017)](https://doi.org/10.1107/S1600576717011712).

[^Altmann]: Altmann, S.L. & Herzig, P., *Point-Group Theory Tables*, Oxford: Clarendon Press (1994).

[^1]: If $D \cong D^*$, i.e., if $UD^*U^\dagger = D$ for some unitary $U$, then conjugating and resubstituting shows that $UU^*$ commutes with every $D(g)$; Schur's lemma then gives $UU^* = c\mathbf{1}$, and unitarity forces $c = \pm 1$, i.e., $U = \pm U^{\mathrm{T}}$. `PSEUDOREAL` is the case $U = -U^{\mathrm{T}}$. An invertible antisymmetric $n\times n$ matrix has $\det U = \det U^{\mathrm{T}} = (-1)^n\det U$, so $n$ must be even.

[^2]: Since $K\Gamma K = \Gamma^*$, we have $T^2 = \Gamma K\Gamma K = \Gamma\Gamma^*$, so $T^2 = -1$ demands $\Gamma\Gamma^* = -\mathbf{1}$ — which $\Gamma = \mathbf{1}$ cannot meet, and which reduces to $\Gamma^2 = -\mathbf{1}$ for real $\Gamma$. The block matrix $\mathrm{i}\sigma_y \otimes \mathbf{1}_n = \left(\begin{smallmatrix} \mathbf{0} & \mathbf{1}_n \\ -\mathbf{1}_n & \mathbf{0}\end{smallmatrix}\right)$ is the simplest real choice that does, and exists only in even dimension — consistent with [^1]. Its phase is a convention: $\mathrm{e}^{\mathrm{i}\alpha}\mathbf{1}$ leaves every $D$ unchanged but sends $\Gamma \to \mathrm{e}^{2\mathrm{i}\alpha}\Gamma$.

[^3]: Writing $D$ in $n\times n$ blocks $D_{ij}$ and using $\Gamma = \left(\begin{smallmatrix} \mathbf{0} & \mathbf{1}_n \\ -\mathbf{1}_n & \mathbf{0}\end{smallmatrix}\right)$, we have $\Gamma D^*\Gamma^\dagger = \left(\begin{smallmatrix} D_{22}^* & -D_{21}^* \\ -D_{12}^* & D_{11}^*\end{smallmatrix}\right)$; equating this to $D$ gives $D_{22} = D_{11}^*$ and $D_{21} = -D_{12}^*$, i.e., the quoted form with $A = D_{11}$ and $B = D_{12}$.
