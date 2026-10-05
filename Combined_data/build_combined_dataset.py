#!/usr/bin/env python3
"""Combine the Double/Triple/Four/Five module-level data into one tidy CSV.

Grain: one row per (pheno_score_method, l_max, network_name, module_abbr).

Sources, all under the repo root:
  <N>_EMT_Mec?Sens/correlations_<tag>[_stateps]/Team_strength_pearson_equal/
      lmax<L>_{dominant,top30,top60}_module_data*.csv   module-level metrics
      lmax<L>_panel_*.csv                               directed team strengths
  <N>_EMT_Mec?Sens/<Tag>_PCA_variance.csv               per-network PCA spectrum

The `_stateps` folders are a separate phenotypic-score calculation and are kept
as a distinct `pheno_score_method` value rather than merged.

The `_noPS` folders are deliberately skipped: they are a strict subset of the
base folders (networks containing the PS module removed), so the base folders
already carry every row. Use the `has_PS` column to reproduce a noPS view.
"""

import csv
import glob
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/.."
ARROW = "→"
VARIANTS = ("dominant", "top30", "top60")
# tag -> (module count, directory holding that tag)
TAGS = {
    "double": (2, "Double_EMT_MecnSens"),
    "triple": (3, "Triple_EMT_MechSens"),
    "four": (4, "Four_EMT_MechSens"),
    "five": (5, "Five_EMT_MechSens"),
}
PCA_FILE = {
    "double": "Double_EMT_MecnSens/Double_PCA_variance.csv",
    "triple": "Triple_EMT_MechSens/Triple_PCA_variance.csv",
    "four": "Four_EMT_MechSens/Four_PCA_variance.csv",
    "five": "Five_EMT_MechSens/Five_PCA_variance.csv",
}


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return [{k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
                for row in csv.DictReader(fh)]


def num(x):
    """Parse a float, treating '', 'NA' and 'NaN' as missing."""
    if x is None or x == "" or x.upper() in ("NA", "NAN"):
        return None
    return float(x)


def pca_metrics(row):
    """PC1_var / NumPC_90 / PCA_entropy from a *_PCA_variance.csv row.

    Definitions verified to reproduce the published correlations_*.csv
    coefficients to 4 decimals across all 4 PCA metrics.
    """
    spectrum = []
    for key in row:
        m = re.fullmatch(r"PC(\d+)_variance_explained", key)
        if m:
            v = num(row[key])
            if v is not None:
                spectrum.append((int(m.group(1)), v))
    spectrum = [v for _, v in sorted(spectrum)]
    total = sum(spectrum)
    if not spectrum or total <= 0:
        return {}
    cum, n90 = 0.0, 0
    for v in spectrum:
        cum += v / total
        n90 += 1
        if cum >= 0.90:
            break
    p = [v / total for v in spectrum if v > 0]
    entropy = -sum(x * math.log(x) for x in p)
    return {
        "PC1_var": spectrum[0],
        "NumPC_90": n90,
        "PCA_entropy": entropy,
        "PCA_entropy_norm": entropy / math.log(len(spectrum)) if len(spectrum) > 1 else None,
        "n_PCs": len(spectrum),
        "n_nodes": row.get("n_nodes"),
        "n_samples": row.get("n_samples"),
    }


def panel_derived(panel_row, modules):
    """Per-network and per-module sums of the directed team-strength panel."""
    intra_sum = 0.0
    inter_sum = 0.0
    out_sum = {m: 0.0 for m in modules}
    in_sum = {m: 0.0 for m in modules}
    for m in modules:
        v = num(panel_row.get(f"T_{m}_intra", ""))
        if v is not None:
            intra_sum += v
    for a in modules:
        for b in modules:
            if a == b:
                continue
            v = num(panel_row.get(f"T_{a}{ARROW}{b}", ""))
            if v is not None:
                inter_sum += v
                out_sum[a] += v
                in_sum[b] += v
    return intra_sum, inter_sum, out_sum, in_sum


FIELDS = [
    "pheno_score_method", "l_max", "n_modules", "network_name", "modules_str",
    "modules_sorted", "has_PS", "module_abbr", "point_label",
    "network_T_s", "module_T_s",
    "P_s_dominant", "P_s_top30", "P_s_top60",
    "impurity", "density", "intra_density", "inter_density",
    "T_module_out_sum", "T_module_in_sum", "T_network_intra_sum", "T_network_inter_sum",
    "n_nodes", "n_samples", "n_PCs",
    "PC1_var", "NumPC_90", "PCA_entropy", "PCA_entropy_norm",
]


def main():
    out_rows = []
    warnings = []

    for tag, (n_modules, parent) in TAGS.items():
        pca = {r["topology"]: pca_metrics(r)
               for r in read_csv(os.path.join(ROOT, PCA_FILE[tag]))}

        for method, suffix in (("newNorm", ""), ("stateps", "_stateps")):
            base = os.path.join(ROOT, parent, f"correlations_{tag}{suffix}",
                                "Team_strength_pearson_equal")
            if not os.path.isdir(base):
                warnings.append(f"missing directory: {base}")
                continue

            lmaxes = sorted({int(re.search(r"lmax(\d+)", os.path.basename(f)).group(1))
                             for f in glob.glob(base + "/lmax*_dominant_module_data*.csv")})
            for lmax in lmaxes:
                # module-level metrics, one file per P_s variant
                variant_rows = {}
                for v in VARIANTS:
                    hits = glob.glob(f"{base}/lmax{lmax}_{v}_module_data*.csv")
                    if len(hits) != 1:
                        warnings.append(f"{tag}/{method}/lmax{lmax}: {len(hits)} files for '{v}'")
                        continue
                    variant_rows[v] = {(r["network_name"], r["module_abbr"]): r
                                       for r in read_csv(hits[0])}
                if "dominant" not in variant_rows:
                    continue

                # directed team-strength panel, one file per lmax
                hits = [f for f in glob.glob(f"{base}/lmax{lmax}_panel*.csv")]
                panel = {}
                if len(hits) == 1:
                    panel = {r["network_name"]: r for r in read_csv(hits[0])}
                else:
                    warnings.append(f"{tag}/{method}/lmax{lmax}: {len(hits)} panel files")

                for (net, mod), r in sorted(variant_rows["dominant"].items()):
                    prow = panel.get(net, {})
                    modules_str = prow.get("modules_str", "")
                    modules = [m.strip() for m in modules_str.split(",") if m.strip()]
                    if not modules:
                        modules = sorted({m for (n, m) in variant_rows["dominant"] if n == net})
                    intra_sum, inter_sum, out_sum, in_sum = panel_derived(prow, modules) \
                        if prow else (None, None, {}, {})

                    row = {
                        "pheno_score_method": method,
                        "l_max": lmax,
                        "n_modules": n_modules,
                        "network_name": net,
                        "modules_str": modules_str,
                        "modules_sorted": ",".join(sorted(modules)),
                        "has_PS": int("PS" in modules),
                        "module_abbr": mod,
                        "point_label": r["point_label"],
                        "network_T_s": r["network_T_s"],
                        "module_T_s": r["module_T_s"],
                        "impurity": r["impurity"],
                        "density": r["density"],
                        "intra_density": r["intra_density"],
                        "inter_density": r["inter_density"],
                        "T_module_out_sum": out_sum.get(mod),
                        "T_module_in_sum": in_sum.get(mod),
                        "T_network_intra_sum": intra_sum,
                        "T_network_inter_sum": inter_sum,
                    }
                    for v in VARIANTS:
                        src = variant_rows.get(v, {}).get((net, mod))
                        row[f"P_s_{v}"] = src["module_pheno_score"] if src else None
                    row.update({k: v for k, v in pca.get(net, {}).items()})
                    if net not in pca:
                        warnings.append(f"no PCA row for {net}")
                    out_rows.append(row)

    out_rows.sort(key=lambda r: (r["pheno_score_method"], r["n_modules"], r["l_max"],
                                 r["network_name"], r["module_abbr"]))
    dest = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "combined_module_level_data.csv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    print(f"wrote {dest}: {len(out_rows)} rows x {len(FIELDS)} cols")
    for w_ in sorted(set(warnings)):
        print("  WARN:", w_, file=sys.stderr)


if __name__ == "__main__":
    main()
