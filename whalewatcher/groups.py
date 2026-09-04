"""Save and load basis-function group definitions as JSON.

A group file looks like::

    {
      "whalewatcher_groups": 1,
      "groups": [
        {"name": "Metal d", "color": "#e41a1c",
         "orbitals": ["0Cu_3dxy", "0Cu_3dxz", "0Cu_3dyz", "0Cu_3dx2y2", "0Cu_3dz2"]},
        {"name": "Ligand N", "color": "#377eb8", "orbitals": ["2N_2px", "2N_2py"]}
      ]
    }

Orbital labels are the parser's ``<atomindex><element>_<shell>`` strings, so
a file written for one calculation applies to any other with the same atom
ordering and basis. Labels that do not exist in a loaded file are kept, so
a definition can be built once and reused; the GUI reports how many were
not matched.
"""

import json

FORMAT_KEY = "whalewatcher_groups"
FORMAT_VERSION = 1


def save_groups(path, groups, colors=None):
    """Write ``groups`` ({name: [labels]}) and optional ``colors``
    ({name: '#rrggbb'}) to ``path``. Group order is preserved."""
    colors = colors or {}
    doc = {
        FORMAT_KEY: FORMAT_VERSION,
        "groups": [
            {"name": name, "color": colors.get(name), "orbitals": list(labels)}
            for name, labels in groups.items()
        ],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")


def load_groups(path):
    """Read a group file. Returns (groups, colors) with the same shapes
    :func:`save_groups` takes. Raises ValueError with a readable message on
    anything that is not a group file."""
    with open(path, "r", encoding="utf-8") as f:
        try:
            doc = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"not valid JSON: {e}") from None

    if not isinstance(doc, dict) or FORMAT_KEY not in doc:
        raise ValueError(f"not a whalewatcher group file (no '{FORMAT_KEY}' key)")
    if doc[FORMAT_KEY] != FORMAT_VERSION:
        raise ValueError(f"unsupported group file version {doc[FORMAT_KEY]!r}")
    entries = doc.get("groups")
    if not isinstance(entries, list):
        raise ValueError("'groups' must be a list")

    groups, colors = {}, {}
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            raise ValueError(f"group {i} has no name")
        name = entry["name"].strip()
        if not name:
            raise ValueError(f"group {i} has an empty name")
        if name in groups:
            raise ValueError(f"group {name!r} appears twice")
        orbitals = entry.get("orbitals", [])
        if not isinstance(orbitals, list) or not all(isinstance(o, str) for o in orbitals):
            raise ValueError(f"group {name!r}: 'orbitals' must be a list of strings")
        # Dedupe within the group, keeping order, the way the GUI does.
        seen = set()
        groups[name] = [o for o in orbitals if not (o in seen or seen.add(o))]
        color = entry.get("color")
        if color is not None:
            if not (isinstance(color, str) and len(color) == 7 and color[0] == "#"):
                raise ValueError(f"group {name!r}: color must look like '#rrggbb'")
            colors[name] = color
    return groups, colors
