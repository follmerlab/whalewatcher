"""ORCA frequency / geometry parser and geometric bond detection."""

import re
import numpy as np

# Covalent radii (Angstrom) for bond detection
COV_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Cu": 1.12,
    "Zn": 1.18, "DEFAULT": 1.0,
}
BOND_TOLERANCE = 0.4   # Angstrom tolerance added to sum of cov. radii


def parse_orca_output(path):
    """Return (atoms, coords, freqs, modes).
    atoms : list of element symbols  (len N)
    coords: np.ndarray shape (N, 3)  Angstrom, final geometry
    freqs : list of floats           (len n_modes)
    modes : np.ndarray shape (n_modes, N, 3)  mass-weighted Cartesian displacements
    """
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()

    # ---- 1. Grab ALL geometry blocks; keep the last one ----
    atoms = []
    coords = []
    i = 0
    last_geom_start = None
    while i < len(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[i]:
            last_geom_start = i
        i += 1

    if last_geom_start is None:
        raise ValueError("No CARTESIAN COORDINATES block found.")

    i = last_geom_start + 2          # skip header + dashes line
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("-"):
            break
        parts = line.split()
        if len(parts) == 4:
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        i += 1

    if not atoms:
        raise ValueError("Could not parse atom coordinates.")

    coords = np.array(coords)
    n_atoms = len(atoms)
    n_dof = 3 * n_atoms

    # ---- 2. Vibrational frequencies ----
    freqs = []
    freq_line_re = re.compile(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1")
    vib_section = False
    for line in lines:
        if "VIBRATIONAL FREQUENCIES" in line:
            vib_section = True
            continue
        if vib_section:
            m = freq_line_re.match(line)
            if m:
                freqs.append(float(m.group(1)))
            elif freqs and line.strip() == "":
                break           # blank line after last frequency → done
    if not freqs:
        raise ValueError("No vibrational frequencies found.")

    n_modes = len(freqs)

    # ---- 3. Normal modes matrix ----
    # Format: printed in blocks of 6 columns, rows 0..n_dof-1
    # Header line looks like:  "     6     7     8     9    10    11"
    modes_flat = np.zeros((n_modes, n_dof))   # modes_flat[mode_idx, dof_idx]

    in_nm_section = False
    current_col_indices = []
    # Header lines look like "                  6          7    ..." (only integers after strip)
    nm_header_re = re.compile(r"^(\d+)(\s+\d+)+$")

    for line in lines:
        if "NORMAL MODES" in line and "--------" not in line:
            in_nm_section = True
            continue
        if in_nm_section:
            if "IR SPECTRUM" in line or "RAMAN SPECTRUM" in line:
                break
            # Detect a column-header line (only integers, no decimals)
            stripped = line.strip()
            if stripped and nm_header_re.match(stripped):
                try:
                    current_col_indices = [int(x) for x in stripped.split()]
                    # Filter to valid mode indices
                    current_col_indices = [c for c in current_col_indices if c < n_modes]
                except ValueError:
                    pass
                continue
            # Detect a data row:  "   ROW   val val val ..."
            parts = stripped.split()
            if len(parts) >= 2 and current_col_indices:
                try:
                    row_idx = int(parts[0])
                    vals = [float(x) for x in parts[1:]]
                    for j, col in enumerate(current_col_indices):
                        if j < len(vals) and row_idx < n_dof and col < n_modes:
                            modes_flat[col, row_idx] = vals[j]
                except ValueError:
                    pass

    # Reshape: modes[mode, atom, xyz]
    modes = modes_flat.reshape(n_modes, n_atoms, 3)

    return atoms, coords, freqs, modes


def detect_bonds(atoms, coords, tol=BOND_TOLERANCE):
    bonds = []
    n = len(atoms)
    for i in range(n):
        ri = COV_RADII.get(atoms[i], COV_RADII["DEFAULT"])
        for j in range(i + 1, n):
            rj = COV_RADII.get(atoms[j], COV_RADII["DEFAULT"])
            d = np.linalg.norm(coords[i] - coords[j])
            if d < (ri + rj + tol):
                bonds.append((i, j))
    return bonds
