"""Exact Loewdin per-MO populations computed from the overlap and MO matrices."""

import re
import numpy as np

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


# ---------- Exact Loewdin populations from S and C ----------
#
# ORCA's LOEWDIN ORBITAL POPULATIONS PER MO table is truncated. The header says
# "THRESHOLD FOR PRINTING IS 0.1%", and that is not a rounding effect: a
# basis-function row is omitted entirely unless it clears 0.1% for at least one
# MO in the printed 6-column block. Measured on a Cu dimer at CP(PPP)/def2-TZVPP
# (2062 basis functions): only 486 of 2062 rows appear in the frontier block and
# the HOMO column sums to 88.6%, not 100%. At that MO, 1904 functions each
# contribute under 0.1% and together account for 15.6%. The deficit also varies
# 8-16% between neighbouring MOs, so it distorts comparisons between MOs and not
# just absolute values. Raising Print[P_OrbPopMO_L] to 2 does not change it; the
# threshold is hard-coded.
#
# So compute the populations directly instead:
#
#     P[u,i] = [ (S^1/2 C)[u,i] ]^2 * 100
#
# with S the AO overlap and C the MO coefficients. No threshold anywhere, so
# each column sums to exactly 100% by construction. Both matrices are printed by
#
#     %output Print[P_Overlap] 1  Print[P_MOs] 2 end
#
# Symmetry-equivalent atoms coming out with identical populations is a good check
# that a given file parsed correctly.

_INT_ROW = re.compile(r"^\s*\d+(\s+\d+)*\s*$")


def _is_rule(s):
    """True for a separator line. The MO block writes these as spaced groups
    ('--------  --------  ...'), the overlap block as one run ('------------'),
    so spaces have to be stripped before testing - not doing that made the
    separator fall through to the end-of-section branch."""
    t = s.replace(" ", "").replace("\t", "")
    return bool(t) and set(t) <= set("-=")


def _is_float(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _read_nbas(path, max_lines=200000):
    """Contracted basis-function count from the ORCA header."""
    pats = (re.compile(r"Number of basis functions\s*\.*\s*(\d+)"),
            re.compile(r"Basis Dimension\s+Dim\s*\.*\s*(\d+)"))
    with open(path, "r", errors="replace") as fh:
        for i, line in enumerate(fh):
            for p in pats:
                m = p.search(line)
                if m:
                    return int(m.group(1))
            if i > max_lines:
                break
    raise ValueError("Could not find the basis-function count in the header. "
                     "Is this a full ORCA output?")


def _read_overlap(path, nbas):
    """Parse the OVERLAP MATRIX block into an (nbas, nbas) array.

    Rows look like:  '      0       1.000000   0.858812  ...'
    """
    S = np.zeros((nbas, nbas))
    seen = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("OVERLAP MATRIX"):
                break
        else:
            raise ValueError("No OVERLAP MATRIX block found. Add "
                             "Print[P_Overlap] 1 to the %output block.")
        cols = None
        for line in fh:
            s = line.strip()
            if not s or _is_rule(s):
                continue
            if _INT_ROW.match(line):
                cols = [int(x) for x in s.split()]
                continue
            parts = s.split()
            try:
                r = int(parts[0])
                vals = [float(x) for x in parts[1:]]
            except (ValueError, IndexError):
                break            # left the matrix, next section reached
            if cols is None or r >= nbas:
                continue
            n = min(len(vals), len(cols))
            S[r, cols[:n]] = vals[:n]
            seen += n
    if seen < nbas:
        raise ValueError(f"OVERLAP MATRIX looks truncated ({seen} elements).")
    return S


def _dash_runs(s):
    """Spans of consecutive '-' in a line, as (start, end_exclusive)."""
    runs, j = [], 0
    while j < len(s):
        if s[j] == "-":
            k = j
            while k < len(s) and s[k] == "-":
                k += 1
            runs.append((j, k))
            j = k
        else:
            j += 1
    return runs


def _read_mo_matrices(path, nbas):
    """Parse MOLECULAR ORBITALS into one (C, labels, energies, occs) per spin.

    Block layout is four header lines then one row per basis function:

        <column indices>
        <orbital energies, Eh>
        <occupation numbers>
        --------  --------
        0Cu  1s    0.000000  0.000002 ...

    A spin boundary is a column-index line that restarts at or below the highest
    index already seen. ORCA prints no spin marker there, just a blank line.

    Columns are FIXED WIDTH and must be sliced, not split. ORCA writes each
    coefficient right-aligned in a 10-character field with no guaranteed
    separator, so adjacent values run together once one is wide enough:

        55C   5s       -10.084917-10.367585 -8.599721-30.294783

    str.split() turns that into '-10.084917-10.367585', which is not a float.
    Splitting silently skipped such rows, which shifted every later row in the
    same block and corrupted the coefficients. The field geometry is taken from
    the dashes line of each block rather than hard-coded, so it follows ORCA if
    the widths ever change. The same treatment is applied to the energy and
    occupation header rows, which can glue for deep core levels.
    """
    out = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("MOLECULAR ORBITALS"):
                break
        else:
            raise ValueError("No MOLECULAR ORBITALS block found. Add "
                             "Print[P_MOs] 2 to the %output block.")

        def _new():
            return (np.zeros((nbas, nbas)), [None] * nbas, [0.0] * nbas,
                    [0.0] * nbas)

        C, labels, energies, occs = _new()
        cols, row, seen_max, started = None, 0, -1, False
        slices, label_end, pending = None, 0, []

        def _slice_vals(raw):
            """Values for one row, or None if any field is not a number."""
            vals = []
            for a, b in slices:
                seg = raw[a:b].strip() if a < len(raw) else ""
                if not seg:
                    vals.append(0.0)
                    continue
                try:
                    vals.append(float(seg))
                except ValueError:
                    return None
            return vals

        for raw in fh:
            raw = raw.rstrip("\n")
            s = raw.strip()
            if not s:
                continue

            if _INT_ROW.match(raw):
                new_cols = [int(x) for x in s.split()]
                if started and new_cols and new_cols[0] <= seen_max:
                    out.append((C, labels, energies, occs))
                    C, labels, energies, occs = _new()
                    seen_max = -1
                cols, row = new_cols, 0
                slices, pending = None, []
                seen_max = max(seen_max, max(new_cols))
                started = True
                continue

            if cols is None:
                continue

            # The block's dashes line defines the field geometry.
            if slices is None and _is_rule(s):
                runs = _dash_runs(raw)
                if len(runs) != len(cols):
                    continue          # a section underline, not a column ruler
                pitch = (runs[1][0] - runs[0][0]) if len(runs) > 1 else \
                        (runs[0][1] - runs[0][0] + 2)
                slices = [(max(b - pitch, 0), b) for _, b in runs]
                label_end = max(runs[0][1] - pitch, 0)
                # Energies then occupations, buffered before the geometry was
                # known. These rows are space-separated and do NOT share the
                # coefficient rows' field alignment, so split() first and only
                # fall back to slicing if the token count disagrees (which
                # happens if a value is wide enough to glue, e.g. the 1e7 Eh
                # virtuals of a decontracted auxiliary basis).
                for k, hraw in enumerate(pending[-2:]):
                    toks = hraw.split()
                    if len(toks) == len(cols) and all(_is_float(t) for t in toks):
                        vals = [float(t) for t in toks]
                    else:
                        vals = _slice_vals(hraw)
                    if vals is None:
                        continue
                    tgt = energies if k == 0 else occs
                    for mo, v in zip(cols, vals):
                        if mo < nbas:
                            tgt[mo] = v
                pending = []
                continue

            if slices is None:
                pending.append(raw)   # energy / occupation rows
                continue

            vals = _slice_vals(raw)
            if vals is None:
                break                 # next section's banner: stop
            name = raw[:label_end].split()
            if len(name) < 2:
                continue
            if row < nbas:
                labels[row] = f"{name[0]}_{name[1]}"
                C[row, cols] = vals
            row += 1

        if started:
            out.append((C, labels, energies, occs))
    if not out:
        raise ValueError("MOLECULAR ORBITALS block present but no coefficients "
                         "were parsed.")
    return out


def parse_orca_exact_loewdin(path, progress=None):
    """Exact Loewdin per-MO populations, computed from S and C.

    Returns the same shape as parse_orca_loewdin_populations_streaming so the
    rest of the app is agnostic: {'spin_up': df, 'spin_down': df} with columns
    MO_n, rows = basis-function labels, and mo_numbers / mo_energies /
    mo_occupations in DataFrame.attrs.

    Unlike the printed-table parser every column sums to exactly 100%.
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for Loewdin analysis")

    def say(msg):
        if progress:
            progress(msg)

    say("reading header…")
    nbas = _read_nbas(path)

    say(f"parsing overlap matrix ({nbas}x{nbas})…")
    S = _read_overlap(path, nbas)

    say("building S^1/2…")
    w, V = np.linalg.eigh(S)
    del S
    # Large decontracted bases are near-linearly-dependent; the smallest
    # eigenvalues can be ~1e-6. Clip at zero so the square root stays real
    # rather than letting round-off produce NaNs.
    np.clip(w, 0.0, None, out=w)
    S_half = (V * np.sqrt(w)) @ V.T
    del V, w

    say(f"parsing MO coefficients ({nbas}x{nbas})…")
    mats = _read_mo_matrices(path, nbas)

    keys = ["spin_up", "spin_down"] if len(mats) > 1 else ["spin_up"]
    results = {}
    for key, (C, labels, energies, occs) in zip(keys, mats):
        say(f"computing populations ({key.replace('spin_', '')} spin)…")
        P = (S_half @ C) ** 2 * 100.0
        del C
        idx = [l if l else f"?_{i}" for i, l in enumerate(labels)]
        df = pd.DataFrame(P, index=idx,
                          columns=[f"MO_{n}" for n in range(nbas)])
        del P
        df.attrs["mo_numbers"] = list(range(nbas))
        df.attrs["mo_energies"] = energies
        df.attrs["mo_occupations"] = occs
        df.attrs["exact"] = True
        df.attrs["restricted"] = len(mats) == 1
        results[key] = df
    say("done")
    return results
