"""The Karakoram and the other ranges side by side: one summary page and a set of single charts, from outputs/ only.

Figures (figures/atlas_*.{png,pdf}):
  atlas_page          one page: map of trends, all ranges on one chart, trend bars with detectability, climate drivers,
                      snowline, elevation bands
  atlas_all_ranges    every range's yearly snow duration as anomalies on one chart, with Sen trend lines
  atlas_small_mult    each range on its own axis with trend line and 95 % band of the slope
  atlas_trend_bars    Sen slope with 95 % confidence interval and the minimum detectable trend for 25 years
  atlas_climate       snow trend against warm-season temperature and cold-season snowfall trends
  atlas_cumulative    cumulative anomaly: when each range ran above or below its own mean
  atlas_bands         trend by 500 m elevation band for each range
  atlas_snowline      snowline altitude per year for each range
  atlas_correlation   year-to-year correlation between ranges: does the Karakoram move with its neighbours?
Writes outputs/atlas_summary.csv with every number the figures show.
"""
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
F, O = R / "figures", R / "outputs"
ORDER = ["Karakoram", "Hindu Kush", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
COL = {"Karakoram": "#000000", "Hindu Kush": "#E69F00", "Himalaya": "#D55E00", "Pamir Mountains": "#0072B2",
       "Tian Shan": "#009E73", "Tibetan Plateau": "#CC79A7"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})

uy = pd.read_csv(O / "unit_year.csv")
uy = uy[uy.level == "range"].copy()
det = pd.read_csv(O / "detectability.csv").set_index("unit")
tr = pd.read_csv(O / "unit_trends.csv").set_index("unit")
bands = pd.read_csv(O / "band_trends.csv")
snl = pd.read_csv(O / "snowline_year.csv")


def sen(y, x):
    s, i, lo, hi = stats.theilslopes(y, x, 0.95)
    return s, i, lo, hi


rows = []
for r in ORDER:
    d = uy[uy.unit == r].sort_values("year")
    s, i, lo, hi = sen(d.duration, d.year)
    ts, *_ = sen(d.t2m_warm, d.year)
    tc, *_ = sen(d.t2m_cold, d.year)
    sf, *_ = sen(d.snowfall_cold, d.year)
    rows.append({"range": r, "km2": d.km2.iloc[0], "mean_days": d.duration.mean(), "sen_days_decade": 10 * s,
                 "ci_lo": 10 * lo, "ci_hi": 10 * hi, "mk_p": tr.loc[r, "p"], "q": tr.loc[r, "q"],
                 "mdt25": det.loc[r, "mdt_25y_days_per_decade"] if r in det.index else np.nan,
                 "warm_t_C_decade": 10 * ts, "cold_t_C_decade": 10 * tc, "snowfall_pct_decade": 1000 * sf / d.snowfall_cold.mean(),
                 "sd_days": d.duration.std()})
S = pd.DataFrame(rows).set_index("range")
S["pct_decade"] = 100 * S.sen_days_decade / S.mean_days
S.round(4).to_csv(O / "atlas_summary.csv")
print(S[["mean_days", "sen_days_decade", "ci_lo", "ci_hi", "pct_decade", "q", "mdt25", "warm_t_C_decade", "snowfall_pct_decade"]].round(2).to_string())


def z(d):
    return (d - d.mean()) / d.std()


def save(fig, name):
    fig.savefig(F / f"{name}.png", bbox_inches="tight"); fig.savefig(F / f"{name}.pdf", bbox_inches="tight"); plt.close(fig)


# --- single charts -------------------------------------------------------------------------------------------
def all_ranges(ax):
    for r in ORDER:
        d = uy[uy.unit == r].sort_values("year")
        zz = z(d.duration)
        s, i, *_ = sen(zz, d.year)
        lw = 2.2 if r == "Karakoram" else 1.0
        ax.plot(d.year, zz, color=COL[r], lw=lw, alpha=0.9 if r == "Karakoram" else 0.55, label=r)
        ax.plot(d.year, i + s * d.year, color=COL[r], lw=lw, ls="--")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("Snow duration anomaly (standard deviations)"); ax.set_xlabel("Hydrological year")
    ax.set_title("All ranges on one scale: yearly anomaly and trend")
    ax.legend(ncol=3, fontsize=6.5, frameon=False, loc="lower left")


def trend_bars(ax):
    y = np.arange(len(ORDER))[::-1]
    for k, r in zip(y, ORDER):
        v = S.loc[r]
        ax.barh(k, v.sen_days_decade, color=COL[r], height=0.6)
        ax.plot([v.ci_lo, v.ci_hi], [k, k], color="k", lw=1)
        ax.plot([-v.mdt25, v.mdt25], [k - 0.38, k - 0.38], color="grey", lw=3, alpha=0.4, solid_capstyle="butt")
        ax.text(max(v.ci_hi, v.mdt25) + 0.5, k, f"q = {v.q:.2f}", va="center", fontsize=6.5)
    ax.set_yticks(y, ORDER); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Trend in snow duration (days per decade)")
    ax.set_title("Trend, 95 % interval (black) and smallest\ntrend detectable in 25 years (grey)")


def climate(ax):
    for r in ORDER:
        v = S.loc[r]
        ax.scatter(v.warm_t_C_decade, v.pct_decade, s=20 + v.km2 / 2e4, color=COL[r], edgecolor="k", lw=0.4, zorder=3)
        ax.annotate(r, (v.warm_t_C_decade, v.pct_decade), fontsize=6.5, xytext=(4, 3), textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("Warm-season air temperature trend (°C per decade)"); ax.set_ylabel("Snow duration trend (% per decade)")
    ax.set_title("Warm-season warming is similar everywhere;\nonly the Karakoram's snow season lengthens")


def bands_heat(ax):
    p = bands.pivot_table(index="range", columns="band_low_m", values="sen_days_per_decade").reindex(ORDER)
    m = np.nanmax(np.abs(p.values))
    im = ax.imshow(p.values, aspect="auto", cmap="RdBu", vmin=-m, vmax=m)
    ax.set_yticks(range(len(p)), p.index); ax.set_xticks(range(len(p.columns)), [f"{int(c / 1000)}" if c % 1000 == 0 else "" for c in p.columns])
    ax.set_xlabel("Elevation band, lower edge (km)")
    plt.colorbar(im, ax=ax, label="days per decade", shrink=0.8)
    done = p.notna().any(axis=1).sum()
    ax.set_title(f"Trend by elevation band ({done} of {len(p)} ranges computed so far)")


def snowline(ax):
    for r in ORDER:
        d = snl[snl.range == r]
        if len(d):
            ax.plot(d.year, d.snowline_m, color=COL[r], lw=2 if r == "Karakoram" else 1, label=r)
    ax.set_ylabel("Snowline altitude (m)"); ax.set_xlabel("Hydrological year"); ax.set_title("Late-summer snowline")


def cumulative(ax):
    for r in ORDER:
        d = uy[uy.unit == r].sort_values("year")
        ax.plot(d.year, z(d.duration).cumsum(), color=COL[r], lw=2.2 if r == "Karakoram" else 1, label=r)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("Cumulative anomaly (standard deviations)"); ax.set_xlabel("Hydrological year")
    ax.set_title("Cumulative anomaly: rising = run of snowy years")


def correlation(ax):
    W = uy.pivot(index="year", columns="unit", values="duration")[ORDER]
    c = W.diff().corr()                                   # first differences: year-to-year co-movement, trend removed
    im = ax.imshow(c.values, cmap="PuOr", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ORDER)), [o.split()[0] for o in ORDER], rotation=45, ha="right"); ax.set_yticks(range(len(ORDER)), ORDER)
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            ax.text(j, i, f"{c.values[i, j]:.2f}", ha="center", va="center", fontsize=6)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title("Year-to-year correlation between ranges")
    return c


def trend_map(ax):
    g = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    g = g[g.level == "range"].copy()
    key = "name" if "name" in g else [c for c in g.columns if g[c].isin(ORDER).any()][0]
    g = g[g[key].isin(ORDER)]
    g["pct"] = g[key].map(S.pct_decade)
    m = np.abs(g.pct).max()
    g.plot(column="pct", cmap="RdBu", vmin=-m, vmax=m, ax=ax, edgecolor="white", lw=0.4, legend=True,
           legend_kwds={"label": "Snow duration trend (% per decade)", "shrink": 0.6})
    g[g[key] == "Karakoram"].boundary.plot(ax=ax, color="k", lw=1.5)
    for _, row in g.iterrows():
        p = row.geometry.representative_point()
        ax.annotate(row[key].replace(" Mountains", ""), (p.x, p.y), fontsize=6.5, ha="center")
    x0, y0, x1, y1 = g.total_bounds; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))
    ax.set_axis_off(); ax.set_title("Snow duration trend 2001–2025 by range (Karakoram outlined)")


for name, fn, size in [("atlas_all_ranges", all_ranges, (7, 4)), ("atlas_trend_bars", trend_bars, (6, 3.5)),
                       ("atlas_climate", climate, (5, 4)), ("atlas_bands", bands_heat, (7, 3)), ("atlas_snowline", snowline, (6, 3.5)),
                       ("atlas_cumulative", cumulative, (7, 3.5)), ("atlas_correlation", correlation, (5, 4.2)), ("atlas_map", trend_map, (8, 5))]:
    fig, ax = plt.subplots(figsize=size)
    out = fn(ax)
    if name == "atlas_snowline":
        ax.legend(fontsize=6.5, frameon=False, ncol=2)
    save(fig, name)
    if name == "atlas_correlation":
        print("\nyear-to-year correlation (first differences):\n", out.round(2).to_string())

# small multiples with the same y-scale in % of each range's mean, so shapes compare
fig, axes = plt.subplots(2, 3, figsize=(9, 5), sharex=True, sharey=True)
for ax, r in zip(axes.flat, ORDER):
    d = uy[uy.unit == r].sort_values("year")
    pct = 100 * d.duration / d.duration.mean()
    s, i, lo, hi = sen(pct, d.year)
    ax.plot(d.year, pct, color=COL[r], marker="o", ms=2.5, lw=1)
    yrs = d.year.values; mid = np.median(yrs); c = np.median(pct)
    ax.fill_between(yrs, c + lo * (yrs - mid), c + hi * (yrs - mid), color=COL[r], alpha=0.15)
    ax.plot(yrs, i + s * yrs, color="k", ls="--", lw=0.8)
    ax.axhline(100, color="grey", lw=0.4)
    ax.set_title(f"{r}: {10 * s:+.1f} % per decade", fontsize=8)
for ax in axes[:, 0]:
    ax.set_ylabel("Snow duration (% of range mean)")
save(fig, "atlas_small_mult")

# --- the one-page summary ----------------------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 14))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1, 1, 1], hspace=0.45, wspace=0.3)
trend_map(fig.add_subplot(gs[0, :]))
all_ranges(fig.add_subplot(gs[1, :]))
trend_bars(fig.add_subplot(gs[2, 0])); climate(fig.add_subplot(gs[2, 1]))
cumulative(fig.add_subplot(gs[3, 0])); bands_heat(fig.add_subplot(gs[3, 1]))
fig.suptitle("Is the Karakoram different? Snow cover duration across the Third Pole, 2001–2025 (MODIS MOD10A1)", fontsize=11, y=0.995)
save(fig, "atlas_page")
print("written", sorted(p.name for p in F.glob("atlas_*.png")))
