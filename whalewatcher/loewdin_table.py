"""Parser for ORCA's printed LOEWDIN ORBITAL POPULATIONS PER MO table.

The table is truncated by ORCA at a hard-coded 0.1% print threshold, so
columns do not sum to 100%. Prefer :func:`whalewatcher.parse_orca_exact_loewdin`
when the overlap and MO matrices are in the file. This parser exists for
output that only has the printed table.

Layout being parsed (restricted output has no SPIN UP / SPIN DOWN lines):

    ----------------------------------
    LOEWDIN ORBITAL POPULATIONS PER MO
    ----------------------------------
    THRESHOLD FOR PRINTING IS 0.1%%
    SPIN UP
                          0         1         2         3         4         5
                     -18.83109  -0.89963  -0.40556  -0.34343  -0.26669   0.03857
                       1.00000   1.00000   1.00000   1.00000   1.00000   0.00000
                      --------  --------  --------  --------  --------  --------
      0O   1s             95.7       2.7       0.3       0.0       0.0       0.2
      ...
    <blank>
                          6         7 ...
    ...
    <blank>
    SPIN DOWN
    ...

Each column block is exactly four header lines - MO numbers, energies in
Hartree, occupations, a dashed rule - followed by one row per basis function
that clears the threshold somewhere in the block, then a blank line. Rows
below threshold in every column of a block are omitted entirely, which is
why the frame is assembled with zeros rather than by concatenating blocks.

Other per-MO tables (MULLIKEN ORBITAL POPULATIONS PER MO, LOEWDIN REDUCED
ORBITAL POPULATIONS PER MO) and the ORBITAL ENERGIES block use the same
SPIN UP / SPIN DOWN markers, so those markers are only honoured inside the
Loewdin table itself. If the file contains the table more than once (an
optimisation with populations at every step) the last one wins.
"""

import numpy as np

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

TABLE_HEADER = "LOEWDIN ORBITAL POPULATIONS PER MO"


class PushbackIterator:
    """Iterator wrapper allowing one-line pushback."""
    def __init__(self, iterator):
        self.iterator = iterator
        self.pushback_line = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.pushback_line is not None:
            line = self.pushback_line
            self.pushback_line = None
            return line
        return next(self.iterator)

    def push(self, line):
        self.pushback_line = line


class LoewdinParseError(ValueError):
    """The Loewdin table was found but a column block did not have the
    expected layout. Raised instead of guessing, because a mis-bound header
    row produces a plausible-looking table with the wrong HOMO/LUMO."""


def parse_orca_loewdin_populations_streaming(filename):
    """Parse ORCA's printed Loewdin MO populations.

    Returns {} if the table is absent, otherwise a dict with 'spin_up' and,
    for unrestricted output, 'spin_down'. Restricted output is returned under
    'spin_up' alone with ``df.attrs['restricted'] = True``.

    Each value is a DataFrame: rows = basis-function labels ('0Cu_3dxy'),
    columns = 'MO_n'. Basis functions ORCA left out of a block (below the
    print threshold) are 0.0, not NaN. DataFrame.attrs carries mo_numbers,
    mo_energies (Hartree), mo_occupations, restricted, exact (False) and
    skipped_rows - the count of data lines that could not be read, which is
    0 for well-formed output and a format-mismatch alarm otherwise.
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for Loewdin parsing")
    results = {}
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        it = PushbackIterator(iter(f))
        for line in it:
            if TABLE_HEADER in line and line.strip() == TABLE_HEADER:
                results = _parse_table(it)
    return results


def _parse_table(it):
    """Parse one LOEWDIN ORBITAL POPULATIONS PER MO table; ``it`` is
    positioned just after the header line."""
    # Skip the closing dashes and the THRESHOLD line, then look at what the
    # first real line is: a SPIN marker (unrestricted) or a block header.
    for line in it:
        s = line.strip()
        if not s or _is_dashes(s) or s.startswith("THRESHOLD"):
            continue
        it.push(line)
        break

    results = {}
    while True:
        line = _next_nonblank(it)
        if line is None:
            break
        s = line.strip()
        if s == "SPIN UP":
            results["spin_up"] = _parse_section(it, restricted=False)
        elif s == "SPIN DOWN":
            results["spin_down"] = _parse_section(it, restricted=False)
        elif _int_row(s) is not None and not results:
            it.push(line)
            results["spin_up"] = _parse_section(it, restricted=True)
        else:
            # Anything else is the next section of the output file.
            it.push(line)
            break
    return {k: v for k, v in results.items() if v is not None}


def _parse_section(it, restricted):
    """Parse column blocks until something that is not a block header."""
    mo_numbers, mo_energies, mo_occupations = [], [], []
    row_index = {}        # label -> row position, in order of first appearance
    blocks = []           # (first_col, {row_pos: [values]})
    skipped = 0

    while True:
        line = _next_nonblank(it)
        if line is None:
            break
        nums = _int_row(line.strip())
        if nums is None:
            it.push(line)
            break

        # Header block. Positional: MO numbers, energies, occupations, dashes.
        # Validate each one; a wrong binding here silently moves the HOMO.
        expected_first = mo_numbers[-1] + 1 if mo_numbers else None
        if any(b - a != 1 for a, b in zip(nums, nums[1:])) or \
                (expected_first is not None and nums[0] != expected_first):
            raise LoewdinParseError(
                f"MO numbers are not consecutive at block starting {nums[0]}: {nums}")
        energies = _float_row(next(it, ""))
        occs = _float_row(next(it, ""))
        rule = next(it, "").strip()
        if energies is None or len(energies) != len(nums):
            raise LoewdinParseError(
                f"expected an energy row of {len(nums)} values after MO numbers {nums}")
        if occs is None or len(occs) != len(nums) or not all(0.0 <= o <= 2.0 for o in occs):
            raise LoewdinParseError(
                f"expected an occupation row of {len(nums)} values in [0, 2] for MOs {nums}")
        if not _is_dashes(rule):
            raise LoewdinParseError(
                f"expected a dashed rule after the occupation row for MOs {nums}")

        first_col = len(mo_numbers)
        mo_numbers.extend(nums)
        mo_energies.extend(energies)
        mo_occupations.extend(occs)

        # Data rows until the blank line that ends the block.
        block = {}
        n = len(nums)
        for line in it:
            s = line.strip()
            if not s:
                break
            parts = s.split()
            if len(parts) != n + 2:
                skipped += 1
                continue
            try:
                vals = [float(x) for x in parts[2:]]
            except ValueError:
                skipped += 1
                continue
            label = f"{parts[0]}_{parts[1]}"
            pos = row_index.setdefault(label, len(row_index))
            block[pos] = vals
        blocks.append((first_col, block))

    if not mo_numbers:
        return None
    n_mo = len(mo_numbers)
    assert len(mo_energies) == n_mo and len(mo_occupations) == n_mo

    data = np.zeros((len(row_index), n_mo))
    for first_col, block in blocks:
        for pos, vals in block.items():
            data[pos, first_col:first_col + len(vals)] = vals

    df = pd.DataFrame(data, index=list(row_index),
                      columns=[f"MO_{k}" for k in mo_numbers])
    df.attrs["mo_numbers"] = mo_numbers
    df.attrs["mo_energies"] = mo_energies
    df.attrs["mo_occupations"] = mo_occupations
    df.attrs["restricted"] = restricted
    df.attrs["exact"] = False
    df.attrs["skipped_rows"] = skipped
    return df


# ---------- line classifiers ----------

def _next_nonblank(it):
    for line in it:
        if line.strip():
            return line
    return None


def _is_dashes(s):
    t = s.replace(" ", "")
    return bool(t) and set(t) == {"-"}


def _int_row(s):
    """[ints] if the line is nothing but integers, else None."""
    parts = s.split()
    if not parts:
        return None
    try:
        return [int(p) for p in parts]
    except ValueError:
        return None


def _float_row(line):
    parts = line.split()
    if not parts:
        return None
    try:
        return [float(p) for p in parts]
    except ValueError:
        return None
