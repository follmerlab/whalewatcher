# whalewatcher

Reads ORCA output files and shows you two things: animated vibrational normal modes, and
how much each part of your molecule contributes to the frontier molecular orbitals.

One file, one window, two tabs. No server, no notebook, no config.

- **Vibrational Modes** — lists every frequency in a `.out`, animates the selected normal
  mode as a looping 3D motion with bonds drawn and an orientation triad.
- **Orbital Analysis** — parses Loewdin per-MO populations, lets you bin basis functions
  into named groups (metal 3d, ligand π, whatever you care about), and plots stacked group
  character across the frontier MOs. Table view with CSV copy-out, plus a HOMO–LUMO gap
  readout.

---

## Install

```bash
pip install numpy matplotlib pandas
```

Tkinter ships with most CPython builds. On bare Linux you may need it separately:

```bash
sudo apt install python3-tk
```

Tested on Python 3.12.4 with numpy 1.26.4, matplotlib 3.8.4, pandas 2.2.2.

pandas is only needed for the Orbital Analysis tab. Without it the vibrational viewer still
works and the orbital tab shows an install prompt instead of crashing.

## Run

```bash
python orca_vib_viewer.py
```

Or hand it a frequency file directly:

```bash
python orca_vib_viewer.py mycomplex.out
```

The positional argument loads into the **Vibrational Modes** tab only. Population logs are
opened from inside the Orbital Analysis tab.

---

## Vibrational Modes tab

Point it at any ORCA output containing a `VIBRATIONAL FREQUENCIES` block — `! Freq`,
`! Opt Freq`, or `! NumFreq` all work.

What the parser pulls out:

| Block | Used for |
|---|---|
| `CARTESIAN COORDINATES (ANGSTROEM)` | geometry (the **last** block in the file, so opt+freq gives the optimized structure) |
| `VIBRATIONAL FREQUENCIES` | the frequency list, in cm⁻¹ |
| `NORMAL MODES` | mass-weighted Cartesian displacement vectors |

Reading the list:

- Grey rows are the near-zero translations and rotations (|ν| < 5 cm⁻¹).
- Red rows are imaginary frequencies, printed negative by ORCA. A saddle point, or a
  geometry that never converged.
- The filter box does a substring match on the formatted number. Typing `16` matches
  `1620.55` and also `216.30` — it is a text filter, not a range query.

Double-click or press Return to animate. Controls under the canvas:

- **Speed** — milliseconds per frame, 5 to 200. Applies live.
- **Amplitude** — peak displacement in Å, 0.05 to 1.5. Changing it restarts the animation.
  The mode vector is rescaled so the largest single-atom displacement equals this value,
  which makes modes visually comparable but means the amplitude is **not** physical.
- **Pause / Resume**.
- **Axis** — toggles the XYZ triad.

Bonds are geometric, not from ORCA. Two atoms are bonded if their separation is under the
sum of covalent radii plus 0.4 Å. Long metal–ligand bonds and dative interactions sometimes
fall outside that; the tolerance is `BOND_TOLERANCE` near the top of the file.

Atom colours and radii cover H, C, N, O, F, P, S, Cl, Br, I and the first-row transition
metals Mn through Zn. Anything else falls back to a tan sphere at 1.0 Å.

---

## Orbital Analysis tab

### Two sources, and why the default is not ORCA's own table

The **Source** selector at the top of the tab picks where populations come from.

| Source | Reads | Columns sum to |
|---|---|---|
| `exact (S,C)` | `OVERLAP MATRIX` + `MOLECULAR ORBITALS` | **100%** |
| `printed table` | `LOEWDIN ORBITAL POPULATIONS PER MO` | 85–92% |
| `auto` (default) | exact if the matrices are present, else the printed table | — |

**Prefer `exact`.** ORCA's own per-MO table is truncated, and not slightly. Its header says
`THRESHOLD FOR PRINTING IS 0.1%`, and that is not a rounding effect: a basis-function row is
omitted *entirely* unless it clears 0.1% for at least one MO in the printed six-column
block. Measured on a Cu dimer at CP(PPP)/def2-TZVPP, 2062 basis functions:

- only 486 of 2062 rows appear in the frontier block
- the HOMO column sums to **88.6%**, not 100
- at that MO, **1904 functions each contribute under 0.1% and together account for 15.6%**
- the deficit varies **8–16% between neighbouring MOs**, so it distorts comparisons *between*
  MOs as well as absolute values

Raising `Print[P_OrbPopMO_L]` to 2 does not help — verified, byte-identical output. The
threshold is hard-coded.

Exact mode sidesteps it by computing the populations rather than reading them:

$$P_{\mu i} = \left[(S^{1/2}C)_{\mu i}\right]^2 \times 100$$

with $S$ the AO overlap and $C$ the MO coefficients. No threshold anywhere, so every column
sums to 100% by construction.

### Input for exact mode

```
%output
  Print[P_Overlap]    1     # AO overlap matrix
  Print[P_MOs]        2     # MO coefficients
  Print[P_OrbPopMO_L] 1     # optional: ORCA's own table, for comparison
end
```

A single-point job on an existing `.gbw` is enough — no SCF needed:

```
! ... UKS moread
%moinp "yourjob.gbw"
%scf MaxIter 0 end
```

Cost: the overlap matrix dominates the file size. The Cu dimer log was 205 MB with it and
40 MB without. Parsing plus the $S^{1/2}$ diagonalisation took **13 s** for 2062 basis
functions; it scales as $O(N^3)$, so expect minutes rather than seconds past ~5000.

If the matrices are absent, `auto` falls back to the printed table and says so in a dialog
rather than silently handing you numbers that read 12% low.

### Accuracy of exact mode

Occupied and frontier MOs come out at 100.000% (mean deviation 1×10⁻⁴). High virtuals
deviate up to ±2%, because ORCA prints coefficients to six decimals and this basis is
near-linearly-dependent — the smallest overlap eigenvalue is 2×10⁻⁶, which amplifies that
truncation. It affects only virtuals hundreds of eV above the LUMO. The **Total** column in
the Table tab shows the sum for every MO, so this is visible rather than assumed.

### Reading either source

Population data is keyed on `SPIN UP` / `SPIN DOWN` in the printed table. Exact mode instead
detects the spin boundary as a restart of the MO column indices, so **restricted
calculations work in exact mode** even though the printed-table parser still requires
unrestricted output ([issue #3](https://github.com/follmerlab/whalewatcher/issues/3)).

One parsing hazard worth knowing if you touch that code: MO coefficients are **fixed-width
and can run together** — `-10.084917-10.367585` is two values, not one — so those rows are
sliced by column position, not `split()`.

The `.pop.log` suffix is a lab convention, not an ORCA default. The file picker accepts
`*.log` and `*.out`, and nothing depends on the name.

Energies are stored in Hartree and converted to eV for display at 27.2114 eV/Ha.

### Workflow

1. **Open .pop.log…** and wait for the progress bar. The window stays responsive; parsing
   runs off the main thread.
2. Column 1 fills with every basis function label found, naturally sorted so `MO_2` comes
   before `MO_10`. The filter box is a case-insensitive substring match — type `Cu` for
   every copper function, `3d` for all d functions.
3. In column 2, type a name and hit **+** to create a group. The name auto-increments, so
   clicking + repeatedly gives you Group 1, Group 2, Group 3. Each group gets a colour from
   a fixed 12-colour cycle; colours are assigned monotonically and are not reused after a
   delete. Double-click a group to rename it, **✕** to delete it.
4. Select a group, select orbitals in column 1, then **→ Add to Group** (or double-click a
   single orbital). Duplicates within a group are dropped silently. Nothing stops you from
   putting the same basis function in two different groups — if you do, its population is
   counted twice and the stack overshoots.
5. Set **Spin** and **n MOs each side**, then **Update Plot**.

### Controls

**Spin** — `up`, `down`, or `both`. `both` draws side-by-side panels on a shared y-axis and
reports both gaps. Note that in `both` mode the Table tab shows the up channel only
([issue #5](https://github.com/follmerlab/whalewatcher/issues/5)).

**n MOs each side** — how deep to reach on either side of the gap. It counts *inclusive* of
the frontier pair: `n = 10` gives HOMO−9 through HOMO and LUMO through LUMO+9, so 20 bars.
Clamped at the ends of the MO range.

HOMO is taken as the last MO with occupation > 0.5 and LUMO as the next one down the list,
which assumes ORCA printed the MOs in aufbau order. Fractional-occupation and
broken-symmetry cases can put the marker in the wrong place — check the Occ column in the
table if a result looks off.

### Reading the plot

Stacked bars, one per frontier MO, x-axis running HOMO−n → LUMO+n. Bar height is summed
Loewdin percentage for that group. The red dashed line sits in the HOMO/LUMO gap.

Anything you did not assign to a group is simply not drawn, so a short bar means either
genuinely low character or basis functions you left out. There is no "remainder" bar.

In **exact** mode a fully assigned stack reaches 100%, so a shortfall is unassigned basis
functions and nothing else — check the **Total** column, which reads 100.0.

In **printed table** mode a fully assigned stack still lands near 88%, and the shortfall
varies by MO, so absolute heights are not interpretable and cross-MO comparisons are skewed.
The Total column is the tell: 100.0 means nothing is missing, ~88 means ORCA truncated.

### Table tab

Same numbers as the plot, one row per MO: label, MO number, energy in eV, occupation, one
percentage column per group, and **Total** — the summed population over *every* basis
function in the file, not just the grouped ones. Total is the quickest check on whether a
short stack means low group character or missing data: 100.0 in exact mode, ~85–92 with
ORCA's printed table. HOMO and LUMO rows are tinted.

- **Show all MOs** switches from the frontier window to every MO in the file, recomputing
  group sums over the full range. Rows outside the frontier are labelled by MO number.
  On a large basis this is a lot of rows.
- **Copy as CSV** puts the whole table on the clipboard, header row included. Paste into
  Excel, Origin, or a plotting script.

**HOMO–LUMO gap** above the tabs reports the gap in eV with both edge energies. In `both`
mode it reports each spin channel separately. This is a bare eigenvalue difference — not a
TDDFT excitation energy, and not something to quote as an optical gap.

---

## Known limitations

Tracked as GitHub issues. The ones most likely to bite you:

| # | Problem |
|---|---|
| [3](https://github.com/follmerlab/whalewatcher/issues/3) | Restricted (closed-shell) population output is not recognised at all |
| [4](https://github.com/follmerlab/whalewatcher/issues/4) | Energies and occupations can silently desync from MO numbers across column blocks |
| [5](https://github.com/follmerlab/whalewatcher/issues/5) | `spin=both` plots both channels but tables only the up channel |
| [6](https://github.com/follmerlab/whalewatcher/issues/6) | Basis-function rows are dropped silently when the column count does not match |
| [7](https://github.com/follmerlab/whalewatcher/issues/7) | No test fixtures, so no ORCA version is verified |
| [8](https://github.com/follmerlab/whalewatcher/issues/8) | Header-role detection can bind the energy row to occupations |
| [15](https://github.com/follmerlab/whalewatcher/issues/15) | Section detection also matches the `ORBITAL ENERGIES` block, wasting a full scan of the file |

Of these, only #15 is confirmed to occur on real output. #4, #6, and #8 are latent paths that
did not trigger on the test file — see the validation comments on each.

Full list: <https://github.com/follmerlab/whalewatcher/issues>

## ORCA compatibility

The frequency and geometry blocks have been stable across ORCA 4, 5, and 6, so the
Vibrational Modes tab is expected to work broadly. It has not been run against a fixture
here.

The Orbital Analysis tab has been validated end to end on one file:

| | |
|---|---|
| Method | `! BP86 def2-TZVPP def2/J def2-TZVPP/C UKS RIJCOSX D3 CPCM` |
| System | Cu dimer, 2062 basis functions, 2062 MOs per spin |
| Size | 205 MB, 2.66M lines, 692 column blocks |
| Parse time | 38.8 s |
| Result | both spin channels recovered, all 2062 basis functions, MO numbers / energies / occupations in sync, HOMO–LUMO 2.374 eV |

What that does **not** cover:

- **Restricted output.** The test file is UKS. Closed-shell is known broken
  ([#3](https://github.com/follmerlab/whalewatcher/issues/3)), not merely unverified.
- **Spin-polarised output.** The test file is multiplicity 1 and its two channels are
  identical, so it confirms the channels parse independently but not that a genuinely
  polarised case is handled.
- **Other ORCA releases.** The four-header-line block layout and the two-token row labels
  (`0Cu  6s` → `0Cu_6s`) hold for this file. Both are assumptions, not documented guarantees.

If your file parses to zero groups or loses rows, that is the first thing to suspect.
Attaching a trimmed sample to
[issue #7](https://github.com/follmerlab/whalewatcher/issues/7) is the fastest way to get
your version supported.

One thing not to mistake for a bug: the top virtuals of a decontracted auxiliary basis can
carry eigenvalues in the millions of Hartree. In the test file MO 2061 sits at
2.36 × 10⁷ Eh, verbatim from the log. That is real, and any future range-checking on
energies has to tolerate it.

## Layout

```
orca_vib_viewer.py    everything — parsers, Tk UI, matplotlib canvases
LICENSE               MIT
```

Roughly: `parse_orca_output` handles geometry, frequencies, and normal modes.
`parse_orca_loewdin_populations_streaming` and its helpers handle populations.
`OrcaVibViewer` is the Tk application; `_build_modes_tab` and `_build_orbital_tab` build the
two halves.

## Credits

Built by the [Follmer Lab](https://github.com/follmerlab).

Contributors: Leland Gee, Alec Follmer.

MIT licensed. See [LICENSE](LICENSE).
