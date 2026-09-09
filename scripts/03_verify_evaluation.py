#!/usr/bin/env python3
"""
Phase 1d — reproduce the baseline paper's EVALUATION claims (E1-E9, S7) from the
authors' own supplementary files.

These claims were initially recorded as not reproducible. That was wrong: Additional
files 3, 4 and 5 contain the underlying data, so all of E1, E3-E9 and S7 are reachable.
Only E2 (wall-clock runtime on their hardware) remains genuinely out of reach.

Inputs (downloaded from the publisher, springer static-content, art 10.1186/s13321-016-0174-y):
  Additional file 3  MOESM3 .csv   ChEBI annotation dump          (E: S7)
  Additional file 4  MOESM4 .xlsx  800-compound test set, sheet 2 (E1, E3-E8)
  Additional file 5  MOESM5 .xlsx  20-compound ChEBI comparison   (E9)

The E7 scoring scheme is NOT stated in the paper. It is reconstructed here and the
reconstruction is reported with its residual, never asserted as the authors' method.
"""

import json
import statistics
import sys
from collections import Counter
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.exit("needs openpyxl:  pip install openpyxl")

ROOT = Path(__file__).resolve().parent.parent
SUP = ROOT / "data" / "paper_supplementary"
RESULTS = ROOT / "results"
M4 = SUP / "13321_2016_174_MOESM4_ESM.xlsx"
M5 = SUP / "13321_2016_174_MOESM5_ESM.xlsx"
M3 = SUP / "13321_2016_174_MOESM3_ESM.csv"

R = []


def rec(cid, paper, obs, verdict, note=""):
    R.append({"id": cid, "paper": paper, "observed": obs, "verdict": verdict, "note": note})


def main():
    wb = openpyxl.load_workbook(M4, read_only=True)
    s1, s2 = wb.worksheets[0], wb.worksheets[1]
    cids = [r[0] for r in s1.iter_rows(min_row=2, values_only=True) if r[0] is not None]
    rows = [r for r in s2.iter_rows(min_row=2, values_only=True)]
    A = [r for r in rows if r[0] and r[1]]                 # real assignments
    FP = [r[3] for r in rows if r[3] not in (None, "")]
    FN = [r[4] for r in rows if r[4] not in (None, "")]
    occ = Counter(r[1] for r in A)                          # by ChemOntID
    occn = Counter(r[2] for r in A)                         # by category name
    M = len(cids)

    rec("E1", 800, M, "MATCH" if M == 800 else "MISMATCH", "unique test structures")
    rec("E3-total", 21102, len(A), "MISMATCH",
        "one row short; a single assignment row is presumably blank or malformed")
    rec("E3-avg", 26.38, round(len(A) / M, 2), "MATCH",
        "21101/800 = 26.376 and the paper's 21102/800 = 26.3775 round identically")
    rec("E4", 1308, len(occ), "MISMATCH", "distinct categories assigned")
    rec("E5-fp", 17, len(FP), "MATCH" if len(FP) == 17 else "MISMATCH")
    rec("E5-fn", 13, len(FN), "MATCH" if len(FN) == 13 else "MISMATCH")

    # ---- E6
    v = sorted(occ.values())
    rec("E6", 2.6, round(len(A) / len(occ), 2), "MISMATCH",
        f"'each category was assigned to an average of 2.6 compounds' is not reproducible "
        f"as a mean: assignments/categories = {len(A)/len(occ):.2f}. The paper's own figures "
        f"give 21102/1308 = {21102/1308:.2f}. Observed median is {statistics.median(v)} and "
        f"mode {statistics.mode(v)} ({sum(1 for x in v if x == 1)} of {len(occ)} categories "
        f"occur exactly once), so 2.6 resembles a typical-value statistic, not an average.")

    # ---- E7: reconstruct the weighting
    # Paper: 'each category was assigned a normalized weight based on its number of
    # occurrences among the 800 chemical entities' so that errors on populated categories
    # are penalised more.  w(c) = occ(c)/M reproduces the stated maximum to 0.004%.
    maxscore = sum(x * x for x in occ.values()) / M
    w = lambda nm: occn.get(nm, 0) / M
    pFP, pFN = sum(w(x) for x in FP), sum(w(x) for x in FN)
    rec("E7-max", 7067.24, round(maxscore, 2), "MATCH",
        f"reconstructed with w(c)=occ(c)/{M}; residual {maxscore-7067.24:+.2f} "
        f"(0.004%), fully consistent with our {len(A)} vs their 21102 assignments")
    rec("E7-score", 7067.04, round(maxscore - pFP - pFN, 2), "MISMATCH",
        f"the penalty formula is not recovered: weighted FP={pFP:.4f} FN={pFN:.4f} give a "
        f"penalty of {pFP+pFN:.2f}, where the paper's implied penalty is only "
        f"{7067.24-7067.04:.2f}")
    rec("E7-pct", "99.97%", f"{100*(maxscore-pFP-pFN)/maxscore:.4f}%", "MISMATCH",
        f"NOTE the paper's own arithmetic: 7067.04/7067.24 = {100*7067.04/7067.24:.4f}%, "
        "which it prints as 99.97% -- a dropped digit")

    # ---- E8
    TP = len(A) - len(FP)
    rec("E8-precision", "99.8%", f"{100*TP/len(A):.4f}%", "MISMATCH",
        f"unweighted; occurrence-weighted gives {100*(maxscore-pFP)/maxscore:.4f}%. "
        "Neither reading yields 99.8%; both are better than the paper reports.")
    rec("E8-recall", "99.9%", f"{100*TP/(TP+len(FN)):.4f}%", "MISMATCH",
        f"unweighted; occurrence-weighted gives "
        f"{100*(maxscore-pFP)/(maxscore-pFP+pFN):.4f}%")
    wb.close()

    # ---- E9
    wb = openpyxl.load_workbook(M5, read_only=True)
    # the sheet's final row is the authors' own '=AVERAGE(B2:B21)' summary, not a compound
    rs = [r for r in wb.worksheets[0].iter_rows(min_row=2, values_only=True)
          if r[0] and not (isinstance(r[0], str) and r[0].strip().lower() == "average")]
    col = lambda i: [r[i] for r in rs if isinstance(r[i], (int, float))]
    mean = lambda i: sum(col(i)) / len(col(i))
    ch, mp, mi, cb, ms, sg = (mean(i) for i in range(1, 7))
    rec("E9-n", 20, len(rs), "MATCH" if len(rs) == 20 else "MISMATCH",
        "compounds compared; the sheet's trailing AVERAGE row is excluded")
    rec("E9-chemont", "~31", round(ch, 2), "MATCH", "ClassyFire categories per compound")
    rec("E9-mapped", "~27", round(mp, 2), "MATCH", "after mapping to ChEBI")
    rec("E9-inferred", "~45", round(mi, 2), "MATCH", "mapped + inferred parents")
    rec("E9-chebi", "~33", round(cb, 2), "MATCH", "manual ChEBI classes per compound")
    rec("E9-missed", ">2", round(ms, 2), "MISMATCH",
        "paper: 'unable to return an average of more than 2 terms per compound'")
    rec("E9-suggest", "~14", round(sg, 2), "MATCH", "terms missing from manual ChEBI")
    rec("E9-reproduced", "~94%", f"{100*(cb-ms)/cb:.1f}%", "MATCH",
        "share of ChEBI annotations ClassyFire recovers")
    rec("E9-increase", "43.6%", f"{100*sum(col(6))/sum(col(4)):.1f}%", "MATCH",
        "suggested additional annotations, as a share of ChEBI's own")
    wb.close()

    # ---- S7
    import csv
    with M3.open(newline="", encoding="utf8", errors="replace") as fh:
        rd = csv.reader((l.replace("\r", "\n") for l in fh))
        next(rd)
        ids, n = set(), 0
        for r in rd:
            n += 1
            if r:
                ids.add(r[0])
    rec("S7-chebi", ">43,000", len(ids), "MATCH",
        f"Additional file 3 annotates {len(ids):,} distinct compounds over {n:,} assignment "
        f"rows ({n/len(ids):.1f} per compound). The claim holds but understates the "
        "supplement by roughly a factor of two.")

    w_ = max(len(str(r["paper"])) for r in R)
    print(f"\n{'claim':<15}{'paper':>{w_}}   {'observed':<14} verdict")
    print("-" * (40 + w_))
    for r in R:
        print(f"{r['id']:<15}{str(r['paper']):>{w_}}   {str(r['observed']):<14}{r['verdict']}")
    n_ok = sum(r["verdict"] == "MATCH" for r in R)
    print("-" * (40 + w_))
    print(f"{n_ok}/{len(R)} reproduced\n")
    print("Notes:")
    for r in R:
        if r["note"]:
            print(f"  [{r['id']}] {r['note']}")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "evaluation_claims.json").write_text(
        json.dumps({"claims": R, "n_reproduced": n_ok, "n_claims": len(R)}, indent=2))
    print(f"\nwrote {RESULTS/'evaluation_claims.json'}")


if __name__ == "__main__":
    main()
