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

### Population mode: what the tab actually consumes

This tab does **not** read a plain ORCA `.out`. It reads a log containing ORCA's Loewdin
**reduced orbital populations per MO** — the table of per-basis-function percentage
contributions to each molecular orbital. That block is not printed by default. Request it
in your input:

```
%output
  Print[P_OrbPopMO_L] 1
end
```

Older and newer ORCA releases differ on keyword spelling; `! LargePrint` or
`! PrintBasis PrintMOs` will also get the block out. Confirm your output contains a section
headed like

```
LOEWDIN ORBITAL POPULATIONS PER MO
```

before assuming the file is usable.

These logs get large fast — hundreds of MB for a big basis on a metal complex — which is
why `*.pop.log` is gitignored and why the parser streams the file in column blocks on a
background thread with a progress bar rather than loading it whole.

The `.pop.log` suffix is a lab convention, not an ORCA default. The file picker accepts
`*.log` and `*.out`, and nothing in the parser depends on the name.

### What the parser expects to find

Population data is read from a `SPIN UP` / `SPIN DOWN` section header. Each column block
carries four header lines before the dashes — MO numbers, MO energies in Hartree,
occupation numbers, then the separator — followed by one row per basis function.

Two consequences worth knowing up front:

- **Unrestricted output only.** The parser keys on `SPIN UP` / `SPIN DOWN`. A closed-shell
  restricted calculation prints the population table with no spin header, and the tab
  reports "No Loewdin sections found". Run `UKS`/`UHF` if you need this tab, or see
  [issue #3](https://github.com/follmerlab/whalewatcher/issues/3).
- **Below-threshold contributions are missing.** ORCA truncates the table at a print
  threshold (commonly 0.1%). Group sums therefore undercount slightly, and a column will
  not add to exactly 100%. Treat the numbers as percentages of the printed population.

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
genuinely low character or basis functions you left out. There is no "remainder" bar. If
your stacks are consistently coming up around 60%, you are probably missing functions
rather than looking at a real result.

### Table tab

Same numbers as the plot, one row per MO: label, MO number, energy in eV, occupation, and
one percentage column per group. HOMO and LUMO rows are tinted.

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

Full list: <https://github.com/follmerlab/whalewatcher/issues>

## ORCA compatibility

Developed against ORCA 5.x output. The frequency and geometry blocks have been stable
across ORCA 4, 5, and 6, so the Vibrational Modes tab is expected to work broadly.

The population table's exact column layout and row-label spelling have changed between
major releases, and the repo carries no fixtures to pin this down — so the Orbital Analysis
tab is only known-good on the files it was written against. If your file parses to zero
groups or drops rows, that is the first thing to suspect. Attaching a trimmed sample to
[issue #7](https://github.com/follmerlab/whalewatcher/issues/7) is the fastest way to get
your version supported.

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
