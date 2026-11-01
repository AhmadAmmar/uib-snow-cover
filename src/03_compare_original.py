"""Compare the new series with the original, and summarise cloud and trends. Reads data/uib_sca_8day.csv.

Saves outputs/comparison.csv and figures/original_vs_new.png and .pdf.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
name = sys.argv[1] if len(sys.argv) > 1 else "uib_sca_8day"
new = pd.read_csv(R / "data" / f"{name}.csv", parse_dates=["date"])
old = pd.read_csv(R / "outputs" / "original_8day.csv", parse_dates=["date"])
new["valid"] = new[["snow", "land", "water", "ice"]].sum(axis=1)
new["total"] = new[["snow", "cloud", "land", "water", "ice", "night", "nodecision", "fill"]].sum(axis=1)
new["cloud_pct"] = 100 * new.cloud / new.total
new["snow_pct_valid"] = 100 * new.snow / new.valid
m = old.merge(new, on="date")
r = stats.pearsonr(m.area, m.snow).statistic
s = stats.linregress(m.snow, m.area)
print(f"matched {len(m)} of {len(old)}; r={r:.3f}; original = {s.slope:.2f} x new + {s.intercept:.0f}; median ratio {(m.area / m.snow).median():.2f}")
print("cloud % by month:", new.groupby(new.date.dt.month).cloud_pct.median().round(1).to_dict())
print("composites with cloud over 10%:", int((new.cloud_pct > 10).sum()), "of", len(new))
new["year"] = new.date.dt.year
y = new[new.groupby("year").year.transform("size") >= 44]
for col in ("snow", "snow_pct_valid"):
    a = y.groupby("year")[col].mean()
    t = stats.theilslopes(a.values, a.index.values)
    k = stats.kendalltau(a.index.values, a.values)
    print(f"{col} annual mean {a.index.min()}-{a.index.max()}: Sen {t.slope:.2f} [{t.low_slope:.2f}, {t.high_slope:.2f}] per year; tau {k.statistic:.2f}, p={k.pvalue:.3f}")
(R / "outputs").mkdir(exist_ok=True)
(R / "figures").mkdir(exist_ok=True)
m[["date", "area", "snow", "cloud_pct"]].to_csv(R / "outputs" / "comparison.csv", index=False)
fig, ax = plt.subplots(figsize=(9, 3.2))
ax.plot(m.date, m.area / 1e3, color="#0072B2", lw=0.9, label="Original")
ax.plot(m.date, m.snow / 1e3, color="#D55E00", lw=0.9, ls="--", label="New, Upper Indus above Tarbela")
ax.set_ylabel("Snow cover area (thousand km²)")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(R / "figures" / "original_vs_new.png", dpi=300)
fig.savefig(R / "figures" / "original_vs_new.pdf")
