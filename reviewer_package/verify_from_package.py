#!/usr/bin/env python3
"""
verify_from_package.py

Rebuilds the manuscript's headline numbers from THIS FOLDER ALONE, using only the stored
scores. It does not import the analysis code, does not touch the network, and makes no model
query. Run it first: if these numbers match the manuscript, the rest of the package is worth
reading.

    python3 verify_from_package.py

Requires pandas and numpy only.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def auc(x: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUC, ties counted as one half."""
    return float(((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum())
                 / (len(x) * len(y)))


def main() -> None:
    d = pd.read_csv(HERE / "analysis_set_membership.tsv", sep="\t")
    st = json.loads((HERE / "stratified.json").read_text())
    rm = json.loads((HERE / "ranking_manifest.json").read_text())
    ok = True

    def check(label: str, got, want) -> None:
        nonlocal ok
        good = got == want
        ok &= good
        print(f"{'OK  ' if good else 'FAIL'}  {label:46s} {got}  (manuscript: {want})")

    sp = d[d.in_primary_splicing_set]
    prim = st["primary_splicing_class"]
    check("splicing class positives", int((sp.label == "POS").sum()), prim["n_pos"])
    check("splicing class negatives", int((sp.label == "NEG").sum()), prim["n_neg"])

    tools = {"SpliceAI": "spliceai_max", "AlphaGenome splicing composite": "ag_splicing_composite_anygene",
             "AVI": "avi", "CADD": "cadd_raw"}
    for tool, col in tools.items():
        g = sp.dropna(subset=[col])
        got = round(auc(g[g.label == "POS"][col].to_numpy(float),
                        g[g.label == "NEG"][col].to_numpy(float)), 3)
        check(f"AUC, {tool}", got, prim["auc"][tool]["auc"])

    cm = sp.dropna(subset=list(tools.values()))
    check("common four-tool subset, positives", int((cm.label == "POS").sum()),
          st["common_subset_all_tools"]["n_pos"])

    g = sp.dropna(subset=["spliceai_max"])
    for t in (0.2, 0.5):
        th = st["thresholds_splicing_class"][f"SpliceAI >= {t}"]
        check(f"SpliceAI >= {t}, true positives",
              int(((g.label == "POS") & (g.spliceai_max >= t)).sum()), th["TP"])
        check(f"SpliceAI >= {t}, no-effect variants flagged",
              int(((g.label == "NEG") & (g.spliceai_max >= t)).sum()), th["FP"])

    p1 = d[d.in_PKD1_only_subset].dropna(subset=["avi"])
    check("AUC, AVI, PKD1 only",
          round(auc(p1[p1.label == "POS"].avi.to_numpy(float),
                    p1[p1.label == "NEG"].avi.to_numpy(float)), 3),
          st["composite_scores_PKD1_only"]["all_tools"]["AVI"]["auc"])

    prov = sp.groupby(["origin", "label"]).size().to_dict()
    check("patient-derived negatives (must be zero)", prov.get(("patient", "NEG"), 0), 0)
    check("patient-derived positives", prov.get(("patient", "POS"), 0),
          st["provenance_by_outcome"]["POS"]["patient"])

    print("\nRanking comparisons, both matched on the PROPORTION of their own background flagged:")
    for name, v in rm["backgrounds"]["PKD1"].items():
        s_, a_ = v["spliceai"], v["avi"]
        e, r = v["eligible_benchmark_variants"], v["recovery"]
        same = s_["flag_proportion_pct"] == a_["flag_proportion_pct"]
        ok &= same
        print(f"  {'OK  ' if same else 'FAIL'}  {name:28s} SpliceAI {r['spliceai_pos']}/{e['n_pos']} at "
              f"{s_['flag_proportion_pct']}% of {s_['alleles_scored']:,} alleles; "
              f"AVI {r['avi_pos']}/{e['n_pos']} at {a_['flag_proportion_pct']}% of "
              f"{a_['alleles_scored']:,}")

    print("\nAll headline numbers reproduce from stored scores." if ok
          else "\nMISMATCH: at least one number does not reproduce.")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
