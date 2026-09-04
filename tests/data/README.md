# Test fixtures

All files here were produced by **ORCA 6.1.1** (macOS, serial) from the `.inp` next to
each `.out`, and are committed unmodified. Regenerate with:

```bash
orca orca611_h2o_rks.inp > orca611_h2o_rks.out
```

| File | What it exercises |
|---|---|
| `orca611_h2o_rks.out` | Restricted closed-shell (RKS). `OVERLAP MATRIX`, `MOLECULAR ORBITALS (RHF, ROHF)`, and a `LOEWDIN ORBITAL POPULATIONS PER MO` table with **no** `SPIN UP` / `SPIN DOWN` header. 24 basis functions. |
| `orca611_oh_uks.out` | Unrestricted doublet (UKS, OH radical). Two MO matrices, `SPIN UP` / `SPIN DOWN` population sections, and an `ORBITAL ENERGIES` block whose `SPIN UP ORBITALS` line is the false positive of issue #15. 19 basis functions, genuinely spin-polarised. |
| `orca611_h2o_optfreq.out` | `Opt Freq`. Several `CARTESIAN COORDINATES (ANGSTROEM)` blocks (the parser must take the last), `VIBRATIONAL FREQUENCIES`, `NORMAL MODES`. 3 atoms, 9 modes, 3 real. |

`.pop.log` is only a lab naming convention. These are raw ORCA stdout, renamed to `.out`,
which is what the parsers are hardened against.
