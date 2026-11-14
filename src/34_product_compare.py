"""Does the product or metric change the trend? Same units and years, four measures:
  duration     daily MOD10A1 snow-cover duration (src/13), our main measure
  sca_clear    8-day MOD10A2 snow-covered area, cloud left out (src/26)
  sca_all      8-day MOD10A2 snow-covered area, cloud counted as no snow (src/26)
  moydgl06     improved Terra-Aqua 8-day snow-covered area, MOYDGL06* (src/31), 2003-2018 only
For each unit: over 2001-2025 (2003-2018 for MOYDGL06 and for a like-for-like run of the others) and over each
published study's years, Sen slope as % of the measure's mean per decade, Mann-Kendall p, and the year-to-year
correlation between measures. Outputs: outputs/product_compare.csv, outputs/product_correlation.csv,
figures/product_compare.png/.pdf
"""
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("det", R / "src" / "22_detectability.py")
det = importlib.util.module_from_spec(spec)
spec.loader.exec_module(det)
O = R / "outputs"
MEASURES = ["duration", "sca_clear", "sca_all", "moydgl06"]


def load():
    u = pd.read_csv(O / "unit_year.csv")[["unit", "year", "duration"]]
    e = pd.read_csv(O / "eightday_year.csv")[["unit", "year", "sca_clear", "sca_all", "cloud"]]
    m = pd.read_csv(O / "moydgl06_year.csv")[["unit", "year", "sca"]].rename(columns={"sca": "moydgl06"}) \
        if (O / "moydgl06_year.csv").exists() else pd.DataFrame(columns=["unit", "year", "moydgl06"])
    return e.merge(u, on=["unit", "year"], how="left").merge(m, on=["unit", "year"], how="left")


def trend(g, col, a, b):
    g = g[g.year.between(a, b)].dropna(subset=[col]).sort_values("year")
    if len(g) < 10:
        return None
    x, y = g.year.to_numpy(float), g[col].to_numpy(float)
    s = det.sen(x, y)
    return dict(n=len(y), pct_per_decade=100 * 10 * s / y.mean(), mk_p=2 * stats.norm.sf(abs(det.mk_z(y[None])[0])))


def main():
    d = load()
    pub = pd.read_csv(R / "references" / "published_trends.csv").dropna(subset=["our_unit"])
    periods = [("2001-2025", 2001, 2025, None), ("2003-2018", 2003, 2018, None)] + \
              [(f"{p.study}", max(int(p.start), 2001), min(int(p.end), 2025), p.our_unit) for p in pub.itertuples()]
    rows = []
    for unit, g in d.groupby("unit"):
        for label, a, b, only in periods:
            if only and only != unit:
                continue
            for col in MEASURES:
                r = trend(g, col, a, b)
                if r:
                    rows.append(dict(unit=unit, period=label, years=f"{a}-{b}", measure=col, **r))
    t = pd.DataFrame(rows)
    t.to_csv(O / "product_compare.csv", index=False)
    wide = t[t.period == "2003-2018"].pivot(index="unit", columns="measure", values="pct_per_decade")
    print("Trend, % of mean per decade, 2003-2018\n", wide.round(1).to_string())
    pv = t[t.period == "2003-2018"].pivot(index="unit", columns="measure", values="mk_p")
    print("p\n", pv.round(3).to_string())

    corr = []
    for unit, g in d.groupby("unit"):
        g = g[g.year.between(2003, 2018)]
        for a in MEASURES:
            for b in MEASURES:
                if a < b and g[[a, b]].dropna().shape[0] >= 10:
                    corr.append(dict(unit=unit, a=a, b=b, r=g[a].corr(g[b])))
    c = pd.DataFrame(corr)
    c.to_csv(O / "product_correlation.csv", index=False)
    print(c.groupby(["a", "b"]).r.median().round(2))

    units = list(wide.index)
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(units) + 1.5), dpi=200)
    style = {"duration": ("o", "#000000"), "sca_clear": ("s", "#0072B2"), "sca_all": ("^", "#E69F00"), "moydgl06": ("D", "#009E73")}
    names = {"duration": "daily duration (ours)", "sca_clear": "8-day area, cloud excluded", "sca_all": "8-day area, cloud as no snow",
             "moydgl06": "MOYDGL06 improved 8-day area"}
    for k, col in enumerate([m for m in MEASURES if m in wide.columns]):
        y = np.arange(len(units)) + (k - 1.5) * 0.15
        sig = pv[col].reindex(units) < 0.05
        ax.scatter(wide[col].reindex(units), y, marker=style[col][0], color=style[col][1], s=22, label=names[col],
                   facecolors=[style[col][1] if s else "white" for s in sig], linewidths=1)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(range(len(units)))
    ax.set_yticklabels(units, fontsize=8)
    ax.set_xlabel("Trend 2003-2018 (% of mean per decade); filled = p < 0.05")
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(R / "figures" / f"product_compare.{ext}", bbox_inches="tight")


if __name__ == "__main__":
    main()
