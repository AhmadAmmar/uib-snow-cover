"""Rerun our snow-duration trend over each published study's own period and region, and set it against what the study reported.

For each study in references/published_trends.csv that maps to one of our units: Sen slope and Mann-Kendall p of our
yearly duration over the study's years (clipped to 2001-2025), and the smallest trend that many years could detect
(src/22, using the unit's noise over 2001-2025). Answers: does our data agree in direction, is it significant, and could
a trend of that size have been detected at all?

Outputs: outputs/published_rerun.csv, figures/published_rerun.png and .pdf.
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

pub = pd.read_csv(R / "references" / "published_trends.csv")
uy = pd.read_csv(R / "outputs" / "unit_year.csv")
noise = pd.read_csv(R / "outputs" / "detectability.csv").set_index("unit")
unit_noise = pd.read_csv(R / "outputs" / "detectability_unit_noise.csv").set_index("lag1")


def mdt_years(n, sd, phi):
    """Detectable slope (days per year) for n years: interpolate the unit-noise table in lag-1, then in log n."""
    ns = np.array([int(c[1:]) for c in unit_noise.columns])
    per_n = [np.interp(np.clip(phi, 0, 0.9), unit_noise.index, unit_noise[c]) for c in unit_noise.columns]
    return sd * np.exp(np.interp(np.log(n), np.log(ns), np.log(per_n)))


rows = []
for _, p in pub.dropna(subset=["our_unit"]).iterrows():
    a, b = max(int(p.start), 2001), min(int(p.end), 2025)
    g = uy[(uy.unit == p.our_unit) & uy.year.between(a, b)].sort_values("year")
    x, y = g.year.to_numpy(float), g.duration.to_numpy(float)
    slope = det.sen(x, y)
    z = det.mk_z(y[None])[0]
    nz = noise.loc[p.our_unit]
    m = mdt_years(len(y), nz.resid_sd, nz.lag1)
    ours = "none" if 2 * stats.norm.sf(abs(z)) >= 0.05 else ("decrease" if slope < 0 else "increase")
    rows.append(dict(study=p.study, doi=p.doi, our_unit=p.our_unit, years=f"{a}-{b}", n=len(y),
                     reported_direction=p.reported_direction, reported_significance=p.reported_significance,
                     reported_metric=p.metric, our_sen_days_per_decade=10 * slope,
                     our_mk_p=2 * stats.norm.sf(abs(z)), our_result=ours,
                     our_direction="decrease" if slope < 0 else "increase",
                     detectable_days_per_decade=10 * m, our_over_detectable=abs(slope) / m))
t = pd.DataFrame(rows)
t.to_csv(R / "outputs" / "published_rerun.csv", index=False)
print(t.drop(columns=["doi", "reported_metric"]).round(3).to_string(index=False))

# Figure: one row per study; our trend over its period with the band a trend must exceed to be detected.
colour = {"decrease": "#D55E00", "increase": "#0072B2", "none": "#555555"}
t = t.iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(8, 0.42 * len(t) + 1.4), dpi=200)
for i, r in t.iterrows():
    ax.fill_betweenx([i - 0.3, i + 0.3], -r.detectable_days_per_decade, r.detectable_days_per_decade, color="#dddddd", lw=0)
    ax.plot(r.our_sen_days_per_decade, i, "o", color="black", ms=5)
    ax.text(1.02, i, r.reported_direction + (" (sig.)" if str(r.reported_significance).startswith("signif") else ""),
            transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=colour[r.reported_direction])
ax.axvline(0, color="black", lw=0.6)
ax.set_yticks(range(len(t)))
ax.set_yticklabels([f"{s} ({u}, {y})" for s, u, y in zip(t.study, t.our_unit, t.years)], fontsize=8)
ax.set_xlabel("Our trend in snow-cover duration over the study's years (days per decade)")
ax.set_title("Dots: our trend. Grey: smaller than this cannot be detected reliably (80% power). Right: what the study reported.",
             fontsize=8, loc="left")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(R / "figures" / f"published_rerun.{ext}", bbox_inches="tight")
