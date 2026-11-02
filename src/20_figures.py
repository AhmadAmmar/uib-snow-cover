"""Every figure, table and quoted number for the Third Pole report, from outputs/ only. Rerun after each analysis step.

Writes figures/*.{pdf,png,svg}, docs/generated/values.tex and docs/generated/table_*.tex.
"""
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymannkendall as mk
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import Resampling, calculate_default_transform, reproject
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy import stats
from statsmodels.stats.multitest import multipletests

import cartopy.crs as ccrs
from research_helpers import carto
from research_helpers.autotex import write_values

R = Path(__file__).resolve().parent.parent
F, G, O = R / "figures", R / "docs/generated", R / "outputs"
MONSOON = Path.home() / "Downloads/monsoon/_derived/master_model_tables/master_annual_selected_predictors.csv"
PROJ = carto.albers(85, 25, 45)
EXTENT = [60, 108, 22, 47]
CLASS = {0: ("Too little data", "#bdbdbd"), 1: ("Rarely snow", "#f2efe6"), 2: ("No detectable trend", "#a6cee3"), 3: ("Gain", "#0072B2"), 4: ("Loss", "#D55E00"), 5: ("Persistent snow or ice", "#222222")}
RANGE_ORDER = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
CAPTION_NOTE = "Preliminary."


def regions():
    g = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    fig, ax = plt.subplots(1, 2, figsize=carto.size("double", 0.42), subplot_kw={"projection": PROJ})
    for a, level, col, title in ((ax[0], "division", "parent", "Ranges and divisions"), (ax[1], "subbasin", "parent", "Sub-basins")):
        carto.map_axes(PROJ, EXTENT, ax=a)
        u = g[g.level == level]
        cats = sorted(u[col].unique())
        pal = {c: carto.QUALITATIVE[i % 8] if level == "division" else plt.get_cmap("tab20")(i % 20) for i, c in enumerate(cats)}
        for c in cats:
            a.add_geometries(u[u[col] == c].geometry, ccrs.PlateCarree(), facecolor=pal[c], edgecolor="white", linewidth=0.15, zorder=2)
        a.set_title(title)
        if level == "division":
            handles = [plt.Rectangle((0, 0), 1, 1, color=pal[c]) for c in cats]
            a.legend(handles, cats, loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=6, ncol=3, handlelength=1)
        else:
            carto.scale_bar(a, "lower right")
    carto.panel_labels(ax)
    carto.save(fig, F / "fig_regions")
    plt.close(fig)


def change_maps(values):
    """Warp each class raster to the map projection at 1 km, masking outside the range, and draw it without further resampling."""
    tifs = sorted((R / "data").glob("change_*.tif"))
    if not tifs:
        return
    g = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    g = g[g.level == "range"]
    cmap = ListedColormap([CLASS[k][1] for k in CLASS])
    norm = BoundaryNorm(np.arange(-0.5, 6.5), 6)
    fig, ax = carto.map_axes(PROJ, None, width="double", ratio=0.62, land=False, step=(2, 2))
    boxes = []
    for p in tifs:
        name = p.stem[7:].replace("_", " ")
        shape = g[g.unit.str.lower() == name].to_crs(PROJ.proj4_init)
        with rasterio.open(p) as s:
            transform, w, h = calculate_default_transform(s.crs, PROJ.proj4_init, s.width, s.height, *s.bounds, resolution=1000)
            out = np.full((h, w), 255, dtype="uint8")
            reproject(s.read(1), out, src_transform=s.transform, src_crs=s.crs, dst_transform=transform, dst_crs=PROJ.proj4_init, resampling=Resampling.nearest)
        inside = geometry_mask(shape.geometry, (h, w), transform, invert=True)
        arr = np.ma.masked_where(~inside | (out == 255), out.astype(float))
        box = (transform.c, transform.c + w * transform.a, transform.f + h * transform.e, transform.f)
        ax.imshow(arr, extent=box, origin="upper", cmap=cmap, norm=norm, transform=PROJ, interpolation="nearest")
        boxes.append(box)
    ax.set_extent([min(b[0] for b in boxes), max(b[1] for b in boxes), min(b[2] for b in boxes), max(b[3] for b in boxes)], PROJ)
    ax.legend([plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="#555555", linewidth=0.3) for _, c in CLASS.values()], [n for n, _ in CLASS.values()], loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=6.5)
    carto.scale_bar(ax, "lower left")
    carto.save(fig, F / "fig_change_classes")
    plt.close(fig)


def trends(t, values):
    rows = []
    for (unit, level, parent), d in t.groupby(["unit", "level", "parent"]):
        d = d.sort_values("year")
        v = d.duration.values
        if np.isnan(v).sum() > 5:
            continue
        clean = v[~np.isnan(v)]
        # the larger p of the plain test and trend-free pre-whitening: the Hamed-Rao correction is unstable on 25 values
        plain, tfpw = mk.original_test(clean), mk.trend_free_pre_whitening_modification_test(clean)
        rows.append(dict(unit=unit, level=level, parent=parent, km2=d.km2.iloc[0], mean_days=np.nanmean(v), slope_days_decade=plain.slope * 10, tau=plain.Tau, p=max(plain.p, tfpw.p), p_plain=plain.p, p_tfpw=tfpw.p))
    x = pd.DataFrame(rows)
    x["q"] = multipletests(x.p, method="fdr_bh")[1]
    x.to_csv(O / "unit_trends.csv", index=False)
    values.update(nUnits=len(x), nUnitsSig=int((x.q < 0.05).sum()), nUnitsRaw=int((x.p < 0.05).sum()))
    return x


def series(t):
    d = t[t.level == "range"]
    units = [u for u in RANGE_ORDER if u in set(d.unit)]
    fig, axes = plt.subplots(len(units), 1, figsize=carto.size("single", 0.55 * max(len(units), 1) + 0.25), sharex=True, squeeze=False)
    for ax, u in zip(axes[:, 0], units):
        s = d[d.unit == u].sort_values("year")
        ax.plot(s.year, s.duration, color=carto.QUALITATIVE[0], marker="o", ms=2)
        slope, intercept, lo, hi = stats.theilslopes(s.duration, s.year)
        ax.plot(s.year, intercept + slope * s.year, color=carto.QUALITATIVE[5], ls="--", lw=0.8)
        ax.set_title(u)
        ax.set_ylabel("Days per year")
    axes[-1, 0].set_xlabel("Hydrological year (October to September)")
    carto.save(fig, F / "fig_series")
    plt.close(fig)


def heatmap(t):
    d = t[t.level == "division"].copy()
    if d.empty:
        return
    d["anom"] = d.duration - d.groupby("unit").duration.transform("mean")
    p = d.pivot_table(index=["parent", "unit"], columns="year", values="anom")
    fig, ax = plt.subplots(figsize=carto.size("double", max(0.3, 0.021 * len(p) + 0.12)))
    lim = np.nanpercentile(np.abs(p.values), 98)
    m = ax.pcolormesh(p.columns, range(len(p)), p.values, cmap=carto.DIVERGING.reversed(), vmin=-lim, vmax=lim, shading="nearest", rasterized=True)
    ax.set_yticks(range(len(p)), [u for _, u in p.index], fontsize=6)
    ax.set_xlabel("Hydrological year")
    fig.colorbar(m, ax=ax, label="Anomaly from unit mean (days)", shrink=0.8)
    carto.save(fig, F / "fig_anomaly_heatmap")
    plt.close(fig)


def attribution(t, values):
    """Spearman correlation of detrended series, per range, with false discovery control. These are associations."""
    m = pd.read_csv(MONSOON)[["year", "rain_jjas_mm", "nino34_djf", "dmi_son"]] if MONSOON.exists() else None
    names = {"t2m_cold": "Cold-season temperature", "t2m_warm": "Warm-season temperature", "snowfall_cold": "Cold-season snowfall", "precip_wy": "Total precipitation"}
    if m is not None:
        names.update(rain_jjas_mm="South Asian JJAS rainfall", nino34_djf="Nino 3.4, DJF", dmi_son="Dipole mode index, SON")
    d = t[t.level == "range"]
    if m is not None:
        d = d.merge(m, on="year", how="left")
    rows = []
    for u, s in d.groupby("unit"):
        for c in names:
            z = s[["year", "duration", c]].dropna()
            if len(z) < 15:
                continue
            a, b = [stats.linregress(z.year, z[k]) for k in ("duration", c)]
            res = stats.spearmanr(z.duration - (a.intercept + a.slope * z.year), z[c] - (b.intercept + b.slope * z.year))
            rows.append(dict(unit=u, driver=names[c], r=res.statistic, p=res.pvalue, n=len(z)))
    a = pd.DataFrame(rows)
    a["q"] = multipletests(a.p, method="fdr_bh")[1]
    a.to_csv(O / "attribution.csv", index=False)
    p = a.pivot(index="driver", columns="unit", values="r").reindex(index=list(names.values()))
    q = a.pivot(index="driver", columns="unit", values="q").reindex(index=list(names.values()))
    fig, ax = plt.subplots(figsize=carto.size("single", 0.26 * len(p) / max(p.shape[1], 1) * 2.2 + 0.25))
    im = ax.imshow(p.values, cmap=carto.DIVERGING.reversed(), vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(p.shape[1]), p.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(p)), p.index)
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            if not np.isnan(p.values[i, j]):
                ax.text(j, i, f"{p.values[i, j]:.2f}" + ("*" if q.values[i, j] < 0.05 else ""), ha="center", va="center", fontsize=6, color="white" if abs(p.values[i, j]) > 0.55 else "#222222")
    fig.colorbar(im, ax=ax, label="Spearman correlation", shrink=0.8)
    carto.save(fig, F / "fig_attribution")
    plt.close(fig)
    values["nDriverSig"] = int((a.q < 0.05).sum())
    return a


def main():
    carto.use_style()
    values = {}
    regions()
    g = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    values.update(nRanges=int((g.level == "range").sum()), nDivisions=int((g.level == "division").sum()), nSubbasins=int((g.level == "subbasin").sum()), domainKm=int(g[g.level == "range"].km2.sum()))
    t = pd.read_csv(O / "unit_year.csv")
    values["rangesDone"] = ", ".join(u for u in RANGE_ORDER if u in set(t.unit))
    values["nRangesDone"] = t[t.level == "range"].unit.nunique()
    change_maps(values)
    x = trends(t, values)
    series(t)
    heatmap(t)
    a = attribution(t, values)
    G.mkdir(parents=True, exist_ok=True)
    for r in x[x.level == "range"].itertuples():
        k = r.unit.replace(" ", "")
        values.update({f"slope {k}": (r.slope_days_decade, 1), f"p {k}": (r.p, 2), f"mean {k}": (r.mean_days, 0)})
    c = pd.concat(pd.read_csv(p).assign(rng=p.stem[7:]) for p in O.glob("change_*.csv"))
    c = c[c.level == "range"]
    total = c.km2.sum()
    for name, tag in (("no detectable trend", "stable"), ("gain", "gain"), ("loss", "loss"), ("persistent", "persistent"), ("too little data", "nodata"), ("rarely snow", "rare")):
        lab = "stable seasonal" if name == "no detectable trend" else name
        values[f"km {tag}"] = int(c[c.cls_name == lab].km2.sum())
    values["km classified"] = int(total)
    write_values(G / "values.tex", values)
    show = x.sort_values(["level", "parent", "km2"], ascending=[False, True, False])
    head = ["Unit", "Parent", "Area (km$^2$)", "Mean (days)", "Sen slope (days per decade)", "$p$", "$q$"]
    body = [" & ".join([u.replace("&", "\\&"), p.replace("Third Pole", "--"), f"{k:,.0f}".replace(",", "{,}"), f"{m:.0f}", f"{s:+.1f}", f"{pp:.3f}", f"{q:.3f}"]) + r" \\" for u, p, k, m, s, pp, q in show[["unit", "parent", "km2", "mean_days", "slope_days_decade", "p", "q"]].values]
    (G / "table_trends.tex").write_text("\\begin{tabular}{llrrrrr}\n\\toprule\n" + " & ".join(head) + " \\\\\n\\midrule\n" + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n", encoding="utf8")
    print(values)
    print(x.sort_values("p").head(8).round(3).to_string())
    print(a.sort_values("p").head(8).round(3).to_string())


if __name__ == "__main__":
    main()
