import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from whalewatcher import parse_orca_exact_loewdin


def _check_channel(df, nbas):
    assert df.shape == (nbas, nbas)
    a = df.attrs
    assert len(a["mo_numbers"]) == len(a["mo_energies"]) == len(a["mo_occupations"]) == nbas
    assert a["mo_numbers"] == list(range(nbas))
    assert np.all(np.diff(a["mo_energies"]) >= -1e-6), "energies not in aufbau order"
    # Every column sums to 100% by construction.
    sums = df.sum(axis=0).values
    assert np.allclose(sums, 100.0, atol=0.05)
    assert a["exact"] is True


def test_restricted(data_dir):
    r = parse_orca_exact_loewdin(data_dir / "orca611_h2o_rks.out")
    assert list(r) == ["spin_up"]
    df = r["spin_up"]
    _check_channel(df, 24)
    assert df.attrs["restricted"] is True
    assert df.attrs["mo_occupations"][:6] == [2.0] * 5 + [0.0]


def test_unrestricted(data_dir):
    r = parse_orca_exact_loewdin(data_dir / "orca611_oh_uks.out")
    assert list(r) == ["spin_up", "spin_down"]
    for df in r.values():
        _check_channel(df, 19)
        assert df.attrs["restricted"] is False
    assert sum(r["spin_up"].attrs["mo_occupations"]) == 5
    assert sum(r["spin_down"].attrs["mo_occupations"]) == 4
    # Genuinely spin polarised: the channels differ.
    assert not np.allclose(r["spin_up"].values, r["spin_down"].values)


def test_labels_and_core_orbital(data_dir):
    df = parse_orca_exact_loewdin(data_dir / "orca611_h2o_rks.out")["spin_up"]
    assert "0O_1s" in df.index and "1H_1s" in df.index
    # MO 0 is the oxygen 1s core orbital.
    assert df.loc["0O_1s", "MO_0"] > 90.0


def test_progress_callback(data_dir):
    msgs = []
    parse_orca_exact_loewdin(data_dir / "orca611_h2o_rks.out", progress=msgs.append)
    assert msgs and msgs[-1] == "done"


def test_missing_matrices_raise(data_dir):
    with pytest.raises(Exception):
        parse_orca_exact_loewdin(data_dir / "orca611_h2o_optfreq.out")
