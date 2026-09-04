import numpy as np
import pytest

from whalewatcher import parse_orca_output, detect_bonds


def test_opt_freq_round_trip(data_dir):
    atoms, coords, freqs, modes = parse_orca_output(data_dir / "orca611_h2o_optfreq.out")
    assert atoms == ["O", "H", "H"]
    assert coords.shape == (3, 3)
    assert len(freqs) == 9
    assert modes.shape == (9, 3, 3)


def test_takes_last_geometry_block(data_dir):
    # Opt Freq prints one geometry per optimisation step. The final one is
    # the optimised structure; the initial O-H distance was 0.958 A.
    atoms, coords, _, _ = parse_orca_output(data_dir / "orca611_h2o_optfreq.out")
    d_oh = np.linalg.norm(coords[1] - coords[0])
    assert 0.96 < d_oh < 0.99


def test_frequencies_are_sane(data_dir):
    _, _, freqs, modes = parse_orca_output(data_dir / "orca611_h2o_optfreq.out")
    assert all(abs(f) < 5 for f in freqs[:6])          # translations/rotations
    bend, sym, asym = freqs[6:]
    assert 1500 < bend < 1700
    assert 3500 < sym < asym < 3900
    # Real modes carry displacement; the zero modes are all-zero in ORCA output.
    assert np.abs(modes[6:]).max() > 0.1
    assert np.abs(modes[:6]).max() == pytest.approx(0.0)


def test_bond_detection(data_dir):
    atoms, coords, _, _ = parse_orca_output(data_dir / "orca611_h2o_optfreq.out")
    assert sorted(detect_bonds(atoms, coords)) == [(0, 1), (0, 2)]


def test_missing_blocks_raise(tmp_path):
    p = tmp_path / "empty.out"
    p.write_text("nothing here\n")
    with pytest.raises(ValueError):
        parse_orca_output(p)
