"""Parser for ORCA's printed LOEWDIN ORBITAL POPULATIONS PER MO table."""

import gc

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


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


def parse_orca_loewdin_populations_streaming(filename, chunk_size=100):
    """Stream-parse ORCA Loewdin MO populations from a .pop.log.

    Returns dict with 'spin_up' and/or 'spin_down' DataFrames.
    Rows = orbital labels (e.g. '0Cu_3dxy'), columns = MO_N.
    DataFrame.attrs carries mo_numbers, mo_energies, mo_occupations.
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for Loewdin parsing")
    results = {}
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        file_iter = PushbackIterator(iter(f))
        for line in file_iter:
            if 'SPIN UP' in line and line.strip().startswith('SPIN'):
                results['spin_up'] = _parse_loewdin_section(file_iter, chunk_size)
            elif 'SPIN DOWN' in line and line.strip().startswith('SPIN'):
                results['spin_down'] = _parse_loewdin_section(file_iter, chunk_size)
    return results


def _parse_loewdin_section(file_handle, chunk_size=100):
    all_chunks, all_mo_nums, all_energies, all_occs = [], [], [], []
    while True:
        chunk_data, spin_line = _parse_column_block(file_handle)
        if chunk_data is None:
            if spin_line is not None:
                file_handle.push(spin_line)
            break
        df_chunk, mo_nums, energies, occs = chunk_data
        if df_chunk is not None and not df_chunk.empty:
            all_chunks.append(df_chunk)
            all_mo_nums.extend(mo_nums)
            all_energies.extend(energies)
            all_occs.extend(occs)
        gc.collect()
        if spin_line is not None:
            file_handle.push(spin_line)
            break
        if df_chunk is None or df_chunk.empty:
            break
    if not all_chunks:
        return None
    full_df = pd.concat(all_chunks, axis=1)
    full_df.attrs['mo_numbers'] = all_mo_nums
    full_df.attrs['mo_energies'] = all_energies
    full_df.attrs['mo_occupations'] = all_occs
    return full_df


def _parse_column_block(file_handle):
    """Parse one block of MO columns from the Loewdin table.

    Returns ((DataFrame, mo_numbers, mo_energies, mo_occupations), spin_line).
    spin_line is set when parsing hits the next SPIN section header.
    """
    mo_numbers, mo_energies, mo_occupations = [], [], []
    orbital_data = {}
    header_lines = []
    in_data = False
    spin_line = None

    for line in file_handle:
        stripped = line.strip()
        if ('SPIN UP' in line or 'SPIN DOWN' in line) and stripped.startswith('SPIN'):
            spin_line = line
            break

        if not stripped:
            if in_data:
                break
            continue

        if 'LOEWDIN REDUCED ORBITAL' in line or 'THRESHOLD' in line:
            if in_data:
                break
            continue

        if not in_data:
            header_lines.append(line)

        if '--------' in line:
            in_data = True
            if len(header_lines) >= 4:
                for offset in [4, 3, 2]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            nums = [int(x) for x in vals]
                            if all(0 <= n < 10000 for n in nums):
                                mo_numbers = nums
                                break
                        except ValueError:
                            continue
                for offset in [3, 2]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            energies = [float(x) for x in vals]
                            if len(energies) == len(mo_numbers):
                                mo_energies = energies
                                break
                        except ValueError:
                            continue
                for offset in [2, 1]:
                    if len(header_lines) >= offset:
                        try:
                            vals = header_lines[-offset].split()
                            occs = [float(x) for x in vals]
                            if len(occs) == len(mo_numbers) and all(0 <= o <= 2 for o in occs):
                                mo_occupations = occs
                                break
                        except ValueError:
                            continue
            continue

        if in_data:
            parts = line.split()
            if len(parts) < 3 or '---' in line:
                continue
            try:
                label = f"{parts[0]}_{parts[1]}"
                pops = [float(x) for x in parts[2:]]
                if len(pops) == len(mo_numbers):
                    orbital_data[label] = pops
            except (ValueError, IndexError):
                continue

    if orbital_data and mo_numbers:
        df = pd.DataFrame.from_dict(
            orbital_data, orient='index',
            columns=[f"MO_{n}" for n in mo_numbers]
        )
        return (df, mo_numbers, mo_energies, mo_occupations), spin_line
    return None, spin_line
