import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from whalewatcher import parse_orca_loewdin_populations_streaming, parse_orca_exact_loewdin


def test_unrestricted(data_dir):
    r = parse_orca_loewdin_populations_streaming(data_dir / "orca611_oh_uks.out")
    assert set(r) == {"spin_up", "spin_down"}
    for df in r.values():
        assert df.shape == (19, 19)
        a = df.attrs
        assert len(a["mo_numbers"]) == len(a["mo_energies"]) == len(a["mo_occupations"]) == 19
        assert a["mo_numbers"] == list(range(19))
    assert sum(r["spin_up"].attrs["mo_occupations"]) == 5
    assert sum(r["spin_down"].attrs["mo_occupations"]) == 4


def test_printed_table_agrees_with_exact(data_dir):
    # On a small basis nothing falls under the 0.1% print threshold for the
    # occupied MOs, so the printed table should match the exact populations
    # to the printed precision (one decimal).
    path = data_dir / "orca611_oh_uks.out"
    printed = parse_orca_loewdin_populations_streaming(path)["spin_up"]
    exact = parse_orca_exact_loewdin(path)["spin_up"]
    cols = [f"MO_{i}" for i in range(5)]
    e = exact.loc[printed.index, cols]
    assert np.abs(printed[cols].values - e.values).max() < 0.15
    assert np.allclose(printed.attrs["mo_energies"][:5], exact.attrs["mo_energies"][:5], atol=1e-4)
