"""Regression tests for the printed-table parser against the failure modes
in issues #3, #4, #6, #8 and #15. Variants are built from the UKS fixture."""
import re

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from whalewatcher import parse_orca_loewdin_populations_streaming as parse
from whalewatcher.loewdin_table import LoewdinParseError, TABLE_HEADER


def _table_text(text):
    """The Loewdin table from the fixture, header through the last block."""
    start = text.index(TABLE_HEADER)
    start = text.rfind("\n", 0, text.rfind("\n", 0, start))   # include the dashes above
    end = text.index("\n\n\n", start) if "\n\n\n" in text[start:] else len(text)
    # The table ends at the next blank-blank or section; find the Mayer
    # header that follows it in the fixture.
    end = text.index("MAYER POPULATION ANALYSIS", start)
    end = text.rfind("\n", 0, text.rfind("*", 0, end))
    return text[start:end]


@pytest.fixture
def uks_text(data_dir):
    return (data_dir / "orca611_oh_uks.out").read_text()


def _write(tmp_path, text):
    p = tmp_path / "variant.out"
    p.write_text(text)
    return p


def test_restricted_table_is_parsed(data_dir):
    # Issue 3: closed-shell output has no SPIN UP / SPIN DOWN lines at all.
    r = parse(data_dir / "orca611_h2o_rks.out")
    assert list(r) == ["spin_up"]
    df = r["spin_up"]
    assert df.attrs["restricted"] is True
    assert df.shape == (24, 24)
    assert df.attrs["mo_occupations"][:6] == [2.0] * 5 + [0.0]
    assert df.attrs["skipped_rows"] == 0


def test_unrestricted_is_flagged(data_dir):
    r = parse(data_dir / "orca611_oh_uks.out")
    assert r["spin_up"].attrs["restricted"] is False
    assert r["spin_down"].attrs["restricted"] is False


def test_orbital_energies_block_is_ignored(data_dir):
    # Issue 15: the ORBITAL ENERGIES block prints 'SPIN UP ORBITALS' long
    # before the population table. The fixture contains it.
    text = (data_dir / "orca611_oh_uks.out").read_text()
    assert "SPIN UP ORBITALS" in text
    r = parse(data_dir / "orca611_oh_uks.out")
    assert r["spin_up"].attrs["mo_energies"][0] == pytest.approx(-18.83109)


def test_other_per_mo_tables_do_not_clobber(uks_text, tmp_path):
    # Real lab output has MULLIKEN ORBITAL POPULATIONS PER MO before the
    # Loewdin table and LOEWDIN REDUCED ORBITAL POPULATIONS PER MO after it,
    # each with their own SPIN UP / SPIN DOWN. The old parser took whichever
    # came last.
    table = _table_text(uks_text)
    mulliken = table.replace(TABLE_HEADER, "MULLIKEN ORBITAL POPULATIONS PER MO")
    mulliken = re.sub(r"(\d+)\.(\d)", "0.0", mulliken.split("--------", 1)[1])
    mulliken = table.split("--------", 1)[0].replace(TABLE_HEADER, "MULLIKEN ORBITAL POPULATIONS PER MO") + "--------" + mulliken
    reduced = table.replace(TABLE_HEADER, "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO")
    reduced = re.sub(r"^(\s*)(\d+)([A-Z][a-z]?)\s+(\d)([a-z]\S*)", r"\1\2 \3 \5", reduced, flags=re.M)

    good = parse(_write(tmp_path, uks_text))
    variant = _write(tmp_path, uks_text.replace(table, mulliken + "\n\n" + table + "\n\n" + reduced))
    r = parse(variant)
    assert set(r) == {"spin_up", "spin_down"}
    for k in r:
        pd.testing.assert_frame_equal(r[k], good[k])
        assert r[k].attrs["mo_energies"] == good[k].attrs["mo_energies"]


def test_last_table_wins_when_repeated(uks_text, tmp_path):
    table = _table_text(uks_text)
    shifted = table.replace("-18.83109", "-99.00000")
    r = parse(_write(tmp_path, uks_text.replace(table, shifted + "\n\n" + table)))
    assert r["spin_up"].attrs["mo_energies"][0] == pytest.approx(-18.83109)


def test_missing_energy_row_raises(uks_text, tmp_path):
    # Issues 4 and 8: a block whose header is not the four-line layout must
    # raise, not silently bind occupations to energies or desync the lists.
    table = _table_text(uks_text)
    bad = uks_text.replace(table, re.sub(r"^\s*-18\.83109 .*\n", "", table, count=1, flags=re.M))
    assert bad != uks_text
    with pytest.raises(LoewdinParseError):
        parse(_write(tmp_path, bad))


def test_energies_that_look_like_occupations_raise(uks_text, tmp_path):
    # Two float rows in [0, 2] of the right length; the second is the real
    # occupation row but the first, positionally, is where energies go.
    # Dropping the occupation row then leaves energies in the occupation
    # slot, which the old range-guess would have accepted for virtuals.
    occ_row = re.compile(r"^\s*1\.00000 +1\.00000 +1\.00000 +1\.00000 +1\.00000 +0\.00000 *\n", re.M)
    table = _table_text(uks_text)
    assert occ_row.search(table)
    bad = uks_text.replace(table, occ_row.sub("", table, count=1))
    with pytest.raises(LoewdinParseError):
        parse(_write(tmp_path, bad))


def test_skipped_rows_are_counted(uks_text, tmp_path):
    good = parse(_write(tmp_path, uks_text))["spin_up"]
    bad = uks_text.replace("  0O   1s             95.7       2.7       0.3       0.0       0.0       0.2\n",
                           "  0O   1s             95.7       2.7       xx        0.0       0.0       0.2\n"
                           "  0O   1s   extra     95.7       2.7       0.3       0.0       0.0       0.2\n", 1)
    df = parse(_write(tmp_path, bad))["spin_up"]
    assert df.attrs["skipped_rows"] == 2
    assert good.attrs["skipped_rows"] == 0
    # Everything else is untouched; the missing row reads as zero.
    assert df.loc["0O_1s", "MO_0"] == 0.0
    assert df.loc["0O_2s", "MO_1"] == good.loc["0O_2s", "MO_1"]


def test_rows_absent_from_a_block_are_zero_not_nan(data_dir):
    # Rows below threshold in every column of a block are omitted by ORCA.
    df = parse(data_dir / "orca611_oh_uks.out")["spin_up"]
    assert not df.isna().any().any()
    assert (df == 0.0).any().any()


def test_no_table_returns_empty(data_dir):
    assert parse(data_dir / "orca611_h2o_optfreq.out") == {}
