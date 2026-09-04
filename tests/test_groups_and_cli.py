import json

import pytest

from whalewatcher import save_groups, load_groups
from whalewatcher.cli import parse_args


def test_group_round_trip(tmp_path):
    groups = {"Metal d": ["0Cu_3dxy", "0Cu_3dz2"], "Ligand": ["2N_2px"], "Empty": []}
    colors = {"Metal d": "#e41a1c", "Ligand": "#377eb8"}
    p = tmp_path / "g.json"
    save_groups(p, groups, colors)
    g2, c2 = load_groups(p)
    assert g2 == groups
    assert list(g2) == list(groups)          # order kept
    assert c2 == colors                      # no color for Empty, none invented


def test_load_dedupes_within_group(tmp_path):
    p = tmp_path / "g.json"
    p.write_text(json.dumps({"whalewatcher_groups": 1,
                             "groups": [{"name": "A", "orbitals": ["x", "y", "x"]}]}))
    assert load_groups(p) == ({"A": ["x", "y"]}, {})


@pytest.mark.parametrize("doc,msg", [
    ("not json", "not valid JSON"),
    ({}, "not a whalewatcher group file"),
    ({"whalewatcher_groups": 2, "groups": []}, "unsupported"),
    ({"whalewatcher_groups": 1, "groups": {}}, "must be a list"),
    ({"whalewatcher_groups": 1, "groups": [{"orbitals": []}]}, "has no name"),
    ({"whalewatcher_groups": 1, "groups": [{"name": "A"}, {"name": "A"}]}, "appears twice"),
    ({"whalewatcher_groups": 1, "groups": [{"name": "A", "orbitals": "x"}]}, "list of strings"),
    ({"whalewatcher_groups": 1, "groups": [{"name": "A", "color": "red"}]}, "#rrggbb"),
])
def test_load_rejects_bad_files(tmp_path, doc, msg):
    p = tmp_path / "g.json"
    p.write_text(doc if isinstance(doc, str) else json.dumps(doc))
    with pytest.raises(ValueError, match=msg):
        load_groups(p)


def test_cli_defaults():
    a = parse_args([])
    assert (a.freq_file, a.pop, a.groups, a.tab) == (None, None, None, "modes")


def test_cli_positional_keeps_working():
    a = parse_args(["run.out"])
    assert a.freq_file == "run.out" and a.tab == "modes"


def test_cli_pop_alone_opens_orbital_tab():
    assert parse_args(["--pop", "x.log"]).tab == "orbital"
    assert parse_args(["--groups", "g.json"]).tab == "orbital"
    assert parse_args(["run.out", "--pop", "x.log"]).tab == "modes"
    assert parse_args(["run.out", "--pop", "x.log", "--tab", "orbital"]).tab == "orbital"


def test_cli_rejects_unknown_tab():
    with pytest.raises(SystemExit):
        parse_args(["--tab", "plots"])
