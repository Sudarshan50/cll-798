"""Where the large input files live.

The 2 GB dataset and the 67 MB sample are too big to ship in the repository, so their
location is resolved at run time rather than hard-coded to one machine. Order:

    $ATLAS_DATA/<name>          explicit override
    <repo>/data/<name>          if the file was placed in-repo
    <repo>/../chemont_project/data/<name>   the layout used during development
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ZENODO = "https://zenodo.org/records/20472700  (doi:10.5281/zenodo.20472700)"


def candidates(name):
    env = os.environ.get("ATLAS_DATA")
    out = [Path(env) / name] if env else []
    return out + [ROOT / "data" / name, ROOT.parent / "chemont_project" / "data" / name]


def find(name):
    """First existing candidate, or None."""
    return next((p for p in candidates(name) if p.exists()), None)


def require(name):
    p = find(name)
    if p is None:
        raise SystemExit(
            f"missing input: {name}\n"
            f"  looked in: {', '.join(str(c.parent) for c in candidates(name))}\n"
            f"  set ATLAS_DATA to the directory holding it, or download from\n"
            f"  {ZENODO}"
        )
    return p
