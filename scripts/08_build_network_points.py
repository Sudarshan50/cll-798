#!/usr/bin/env python3
"""
Compound points for the radial-network view (mode 3).

Each class gets a halo of its own compounds around its node. Points are written in the
same order as the per-class [off, cnt] ranges already recorded in network_layout.json,
so the renderer can issue one draw call per VISIBLE class and never upload or draw the
rest. That is what keeps density at a few points per pixel instead of 1,500.
"""
import json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "frontend" / "public"
OUT.mkdir(parents=True, exist_ok=True)

L = json.loads((RES / "network_layout.json").read_text())
nodes = L["nodes"]
TOTAL = L["meta"]["points"]

order = sorted(nodes, key=lambda n: n["off"])
assert order[0]["off"] == 0
run = 0
for n in order:
    assert n["off"] == run, f"gap in point ranges at {n['n']}"
    run += n["cnt"]
assert run == TOTAL, f"{run} != {TOTAL}"

# K chosen so density is ~30 points/px at zoom 1500 (about 3x the default view):
#   K = sqrt(1 / (pi * 30 * 1500^2))
K_HALO = math.sqrt(1.0 / (math.pi * 30.0 * 1500.0 ** 2))
print(f"halo constant K = {K_HALO:.3e}  -> uniform density ~30 pts/px at zoom 1500")

# extent must come from the halo the loop below actually emits; 4*nr understated it
# (max hr 0.094 vs max 4*nr 0.038), leaving only a 0.35% margin before the int16 guard
def _halo(n):
    return max(K_HALO * math.sqrt(n["cnt"]), n["nr"] * 1.05) if n["cnt"] else n["nr"] * 1.05
ext = max(max(abs(n["px"]), abs(n["py"])) + _halo(n) for n in nodes)
S = 32600.0 / ext
print(f"extent {ext:.3f}  int16 scale {S:.1f}")

GA = math.pi * (3 - math.sqrt(5))
CH = 2_000_000
written = 0
with (OUT / "points.bin").open("wb") as fh:
    for n in order:
        c = n["cnt"]
        if c == 0:
            continue
        # Halo radius proportional to sqrt(compounds) makes on-screen density the SAME
        # for every class: density = c / (pi*(K*sqrt(c)*zoom)^2) = 1/(pi*K^2*zoom^2),
        # which is independent of c. Tying it to the node radius (the previous version)
        # left the biggest class at 5,813 points/pixel while small ones were empty.
        rr = max(K_HALO * math.sqrt(c), n["nr"] * 1.05)
        n["hr"] = round(rr, 6)          # frontend culls against this exact value
        done = 0
        while done < c:
            m = min(CH, c - done)
            j = np.arange(done, done + m, dtype=np.float64)
            th = j * GA
            rq = rr * np.sqrt((j + 0.5) / c)
            xf = (n["px"] + np.cos(th) * rq) * S
            yf = (n["py"] + np.sin(th) * rq) * S
            if np.abs(xf).max() > 32767 or np.abs(yf).max() > 32767:
                raise SystemExit("int16 overflow")
            np.stack([xf.astype(np.int16), yf.astype(np.int16)], 1).tofile(fh)
            done += m
            written += m
assert written == TOTAL, f"{written} != {TOTAL}"
p = OUT / "points.bin"
print(f"wrote {p}  {p.stat().st_size/1e6:.1f} MB  ({written:,} points)")

chk = np.fromfile(p, dtype=np.int16, count=4_000_000).reshape(-1, 2)
assert (np.abs(chk).sum(1) == 0).mean() < 0.02, "points collapsed"
print(f"check: x {chk[:,0].min()}..{chk[:,0].max()}  y {chk[:,1].min()}..{chk[:,1].max()}")

for n in nodes:
    n.setdefault("hr", round(max(n["nr"] * 1.05, 0.0), 6))
import hashlib
h = hashlib.sha256()
with p.open("rb") as fh:
    for blk in iter(lambda: fh.read(1 << 22), b""):
        h.update(blk)
# the frontend requests points.bin?v=<version>, so an immutable year-long cache is
# safe: regenerating the cloud changes the version and therefore the URL
L["meta"]["version"] = h.hexdigest()[:12]
print(f"content version: {L['meta']['version']}")
L["meta"]["scale"] = S
L["meta"]["halo_k"] = K_HALO
(RES / "network_layout.json").write_text(json.dumps(L, separators=(",", ":")))
(OUT / "layout.json").write_text(json.dumps(L, separators=(",", ":")))
print("layout.json updated with the matching int16 scale")
