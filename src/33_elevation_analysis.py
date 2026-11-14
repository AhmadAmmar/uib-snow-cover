"""Elevation-dependent change in snow duration, from outputs/band_year.csv (src/24).

Per range and 500 m band (bands with at least 200 pixels in every year): Sen slope, Mann-Kendall p, false discovery
control across all bands, and the detectable trend for 25 years from the band's own noise (src/22).
Snowline per range and year: the elevation where duration first reaches half a year (182.5 days), interpolated
linearly between band midpoints; its Sen slope and p.
Outputs: outputs/band_trends.csv, outputs/snowline_year.csv, outputs/snowline_trends.csv, figures/elevation_bands.png/.pdf
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
unit_noise = pd.read_csv(R / "outputs" / "detectability_unit_noise.csv").set_index("lag1")
HALF = 182.5
MIN_PIX = 200


def trend(x, y):
    b = det.sen(x, y)
    z = det.mk_z(y[None])[0]
    return b, 2 * stats.norm.sf(abs(z))


def mdt25(y, x, b):
    res = y - b * x
    res -= res.mean()
    sd, phi = res.std(ddof=1), np.corrcoef(res[:-1], res[1:])[0, 1]
    return sd * np.interp(np.clip(phi, 0, 0.9), unit_noise.index, unit_noise["n25"])


def snowline(g):
    g = g.sort_values("band_low_m")
    mid, dur = g.band_low_m.to_numpy() + 250, g.duration.to_numpy()
    above = np.nonzero(dur >= HALF)[0]
    if not len(above) or above[0] == 0:
        return np.nan
    k = above[0]
    return mid[k - 1] + (HALF - dur[k - 1]) * (mid[k] - mid[k - 1]) / (dur[k] - dur[k - 1])


def main():
    b = pd.read_csv(R / "outputs" / "band_year.csv")
    full = b.groupby(["range", "band_low_m"]).filter(lambda g: len(g) == 25 and g.pixels.min() >= MIN_PIX)
    rows = []
    for (r, lo), g in full.groupby(["range", "band_low_m"]):
        g = g.sort_values("year")
        x, y = g.year.to_numpy(float), g.duration.to_numpy(float)
        s, p = trend(x, y)
        rows.append(dict(range=r, band_low_m=lo, mean_days=y.mean(), pixels=int(g.pixels.median()),
                         sen_days_per_decade=10 * s, mk_p=p, mdt25_days_per_decade=10 * mdt25(y, x, s)))
    t = pd.DataFrame(rows)
    t["q"] = stats.false_discovery_control(t.mk_p)
    t["observed_over_mdt25"] = t.sen_days_per_decade.abs() / t.mdt25_days_per_decade
    t.to_csv(R / "outputs" / "band_trends.csv", index=False)
    print(len(t), "bands;", (t.mk_p < 0.05).sum(), "p<0.05;", (t.q < 0.05).sum(), "q<0.05")

    sl = b.groupby(["range", "year"]).apply(snowline, include_groups=False).rename("snowline_m").reset_index()
    sl.to_csv(R / "outputs" / "snowline_year.csv", index=False)
    st = []
    for r, g in sl.dropna().groupby("range"):
        if len(g) < 20:
            continue
        s, p = trend(g.year.to_numpy(float), g.snowline_m.to_numpy(float))
        st.append(dict(range=r, years=len(g), mean_m=g.snowline_m.mean(), sen_m_per_decade=10 * s, mk_p=p))
    st = pd.DataFrame(st)
    st.to_csv(R / "outputs" / "snowline_trends.csv", index=False)
    print(st.round(2).to_string(index=False))

    ranges = sorted(t.range.unique())
    fig, axes = plt.subplots(1, len(ranges), figsize=(2.1 * len(ranges), 3.6), dpi=200, sharey=True, squeeze=False)
    for ax, r in zip(axes[0], ranges):
        g = t[t.range == r].sort_values("band_low_m")
        mid = g.band_low_m + 250
        ax.fill_betweenx(mid, -g.mdt25_days_per_decade, g.mdt25_days_per_decade, color="#dddddd", lw=0, label="not detectable")
        ax.plot(g.sen_days_per_decade, mid, "-o", color="#0072B2", ms=3, lw=1, label="observed trend")
        sig = g[g.q < 0.05]
        ax.plot(sig.sen_days_per_decade, sig.band_low_m + 250, "o", color="#D55E00", ms=5, label="significant after correction")
        ax.axvline(0, color="black", lw=0.5)
        ax.set_title(r, fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0][0].set_ylabel("Elevation (m, band midpoint)", fontsize=8)
    fig.supxlabel("Trend in snow-cover duration, 2001-2025 (days per decade)", fontsize=8)
    axes[0][0].legend(fontsize=6, frameon=False, loc="lower left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(R / "figures" / f"elevation_bands.{ext}")


if __name__ == "__main__":
    main()
