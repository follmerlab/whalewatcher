"""Drives the real Tk widgets through a load / group / plot / table cycle.

Skipped where there is no display (CI). Run locally with ``pytest -q``.
"""
import time
from unittest import mock

import pytest

pytest.importorskip("pandas")
tk = pytest.importorskip("tkinter")


@pytest.fixture
def app():
    import orca_vib_viewer as ovv
    try:
        a = ovv.OrcaVibViewer()
    except tk.TclError as e:          # no display
        pytest.skip(f"no display: {e}")
    a.withdraw()
    yield a
    a.destroy()


def _load_and_wait(app, path, source="printed table", timeout=30):
    app.pop_mode_var.set(source)
    app._load_pop_file(str(path))
    t0 = time.time()
    while app.pop_file_label.cget("text") in ("Loading…",) or app.loewdin_data is None:
        app.update()
        time.sleep(0.02)
        if app.orb_status_label.cget("text").startswith(("Error", "No Loewdin")):
            raise AssertionError(app.orb_status_label.cget("text"))
        if time.time() - t0 > timeout:
            raise AssertionError("timed out loading")


@pytest.mark.parametrize("fixture,source,expect_spin", [
    ("orca611_h2o_rks.out", "printed table", "closed-shell"),
    ("orca611_h2o_rks.out", "exact (S,C)", "closed-shell"),
    ("orca611_oh_uks.out", "printed table", "up"),
    ("orca611_oh_uks.out", "exact (S,C)", "up"),
])
def test_load_group_plot_table(app, data_dir, fixture, source, expect_spin):
    import orca_vib_viewer as ovv
    with mock.patch.object(ovv.messagebox, "showwarning") as warn, \
            mock.patch.object(ovv.messagebox, "showinfo"), \
            mock.patch.object(ovv.messagebox, "showerror") as err:
        _load_and_wait(app, data_dir / fixture, source)
        assert app.orb_spin_var.get() == expect_spin
        assert not warn.called and not err.called

        app.orb_groups = {
            "O": [l for l in app._avail_orbitals if l.startswith("0O")],
            "H": [l for l in app._avail_orbitals if "H_" in l],
        }
        app.group_colors = {"O": "#e41a1c", "H": "#377eb8"}
        app._update_orbital_plot()
        app.update()
        assert not err.called
        assert "eV" in app.gap_label.cget("text")
        header, *rows = app._table_rows
        assert header[:4] == ("Orbital", "MO#", "Energy (eV)", "Occ")
        assert "O" in header and "H" in header and "Total" in header
        assert any(r[0] == "HOMO" for r in rows) and any(r[0] == "LUMO" for r in rows)


def test_freq_tab_loads(app, data_dir):
    app._load_file(str(data_dir / "orca611_h2o_optfreq.out"))
    app.update()
    assert app.atoms == ["O", "H", "H"]
    assert len(app.freqs) == 9
    assert app.listbox.size() == 9
