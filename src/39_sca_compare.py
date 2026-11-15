"""Snow-covered area beside snow-cover days, for every range, from outputs/ only.

Snow-covered area (SCA) here is the mean over each hydrological year of the snow share of clear pixels in the 8-day
MOD10A2 maximum snow extent (outputs/eightday_year.csv, sca_clear), times the range area. Snow-cover days come from
daily MOD10A1 (outputs/unit_year.csv). The two are related: a range's mean snow days divided by 365 is its mean snow
share. They differ in product (8-day maximum extent shows more snow than daily), in cloud handling and in what they
weight. Agreement between them is a check; disagreement is a finding to explain.

Writes figures/atlas_sca.{png,pdf} and outputs/sca_summary.csv.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymannkendall as mk
from scipy import stats

R = Path(__file__).resolve().parent.parent
F, O = R / "figures", R / "outputs"
ORDER = ["Karakoram", "Hindu Kush", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
COL = {"Karakoram": "#000000", "Hindu Kush": "#E69F00", "Himalaya": "#D55E00", "Pamir Mountains": "#0072B2",
       "Tian Shan": "#009E73", "Tibetan Plateau": "#CC79A7"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})

e = pd.read_csv(O / "eightday_year.csv"); e = e[e.level == "range"]
u = pd.read_csv(O / "unit_year.csv"); u = u[u.level == "range"]
d = e.merge(u[["unit", "year", "km2", "duration"]], on=["unit", "year"])
d["sca_km2"] = d.sca_clear * d.km2
d["days_share"] = d.duration / 365

rows = []
for r in ORDER:
    x = d[d.unit == r].sort_values("year")
    s, _, lo, hi = stats.theilslopes(x.sca_km2, x.year, 0.95)
    t = max(mk.original_test(x.sca_km2.values).p, mk.trend_free_pre_whitening_modification_test(x.sca_km2.values).p)  # project rule: the larger p
    rd = stats.pearsonr(x.sca_clear, x.days_share)[0]
    rd_diff = stats.pearsonr(np.diff(x.sca_clear), np.diff(x.days_share))[0]
    ds, *_ = stats.theilslopes(x.duration, x.year)
    rows.append({"range": r, "km2": x.km2.iloc[0], "mean_sca_pct": 100 * x.sca_clear.mean(), "mean_sca_km2": x.sca_km2.mean(),
                 "sen_km2_decade": 10 * s, "ci_lo": 10 * lo, "ci_hi": 10 * hi, "pct_decade": 1000 * s / x.sca_km2.mean(),
                 "mk_p": t, "days_pct_decade": 1000 * ds / x.duration.mean(), "mean_days_share_pct": 100 * x.days_share.mean(),
                 "r_sca_days": rd, "r_sca_days_yearly_change": rd_diff, "mean_cloud_pct": 100 * x.cloud.mean()})
S = pd.DataFrame(rows).set_index("range")
from statsmodels.stats.multitest import multipletests
S["q"] = multipletests(S.mk_p, method="fdr_bh")[1]
S.round(4).to_csv(O / "sca_summary.csv")
print(S[["mean_sca_pct", "mean_days_share_pct", "mean_sca_km2", "sen_km2_decade", "ci_lo", "ci_hi", "pct_decade",
         "days_pct_decade", "mk_p", "q", "r_sca_days", "r_sca_days_yearly_change", "mean_cloud_pct"]].round(3).to_string())

fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
a = ax[0, 0]
for r in ORDER:
    x = d[d.unit == r].sort_values("year")
    a.plot(x.year, 100 * x.sca_km2 / x.sca_km2.mean(), color=COL[r], lw=2.2 if r == "Karakoram" else 1, label=r)
a.axhline(100, color="grey", lw=0.5); a.set_ylabel("Snow-covered area (% of range mean)"); a.set_xlabel("Hydrological year")
a.set_title("Snow-covered area, every range on one scale"); a.legend(fontsize=6.5, ncol=2, frameon=False)

a = ax[0, 1]
y = np.arange(len(ORDER))[::-1]
for k, r in zip(y, ORDER):
    v = S.loc[r]
    a.barh(k + 0.18, v.pct_decade, height=0.34, color=COL[r])
    a.barh(k - 0.18, v.days_pct_decade, height=0.34, color=COL[r], alpha=0.35, hatch="///", edgecolor="white")
a.set_yticks(y, ORDER); a.axvline(0, color="k", lw=0.6); a.set_xlabel("Trend (% per decade)")
a.set_title("Trend: area (solid) and snow-cover days (hatched)")

a = ax[1, 0]
for r in ORDER:
    v = S.loc[r]
    a.barh(ORDER.index(r), v.mean_sca_km2 / 1e3, color=COL[r])
    a.text(v.mean_sca_km2 / 1e3, ORDER.index(r), f"  {v.mean_sca_pct:.0f} % of range", va="center", fontsize=6.5)
a.invert_yaxis(); a.set_yticks(range(len(ORDER)), ORDER); a.set_xlabel("Mean snow-covered area (thousand km²)")
a.set_title("How much snow: mean area and share of each range")

a = ax[1, 1]
for r in ORDER:
    x = d[d.unit == r]
    a.scatter(100 * x.days_share, 100 * x.sca_clear, s=8, color=COL[r], label=r)
lim = [0, 100 * max(d.sca_clear.max(), d.days_share.max()) * 1.05]
a.plot(lim, lim, color="grey", lw=0.6, ls=":"); a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("Daily product: mean snow days / 365 (%)"); a.set_ylabel("8-day product: snow share of clear pixels (%)")
a.set_title("The two measures, each range-year")
plt.tight_layout()
fig.savefig(F / "atlas_sca.png", bbox_inches="tight"); fig.savefig(F / "atlas_sca.pdf", bbox_inches="tight")
print("written figures/atlas_sca")
