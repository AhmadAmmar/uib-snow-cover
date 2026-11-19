"""Batch jobs for what is too large to compute interactively, kept as table assets and then read back.

    python src/16_full_resolution_and_subbasins.py start   Tibetan Plateau p histogram at 500 m; sub-basin yearly duration at 1 km
    python src/16_full_resolution_and_subbasins.py read    writes outputs/tibet_fdr_500m.csv, outputs/subbasin_year.csv, outputs/subbasin_trends.csv
"""
import sys
from pathlib import Path

import ee
import numpy as np
import pandas as pd
import pymannkendall as mk
from statsmodels.stats.multitest import multipletests

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
A = "projects/ee-ahmadammargeo/assets/"
ROOT = A + "tp_snow_v1"
RANGES = ["hindu_kush", "karakoram", "himalaya", "pamir_mountains", "tian_shan", "tibetan_plateau"]
YEARS = list(range(2001, 2026))
TOP, BINS, Q, SCALE = 0.002, 2000, 0.05, 463.3127
MIN_VALID = 120
TIBET = A + "tp_table_tibet_hist"
SUB = A + "tp_table_subbasin_year"


def start(which="all"):
    """Start exports of Tibetan Plateau snow duration at full resolution and by sub-basin."""
    rng = ee.FeatureCollection(A + "tp_regions_v2_range").filter(ee.Filter.eq("unit", "Tibetan Plateau")).geometry()
    m = ee.Image(ROOT + "_mk/tibetan_plateau")
    p = m.select("p").updateMask(m.select("n").gte(20)).rename("p")
    # native grid of the stored image
    props = ee.Dictionary({"count": p.reduceRegion(ee.Reducer.count(), rng, SCALE, maxPixels=1e11, tileScale=16).get("p"), "hist": ee.List(p.updateMask(p.lt(TOP)).reduceRegion(ee.Reducer.fixedHistogram(0, TOP, BINS), rng, SCALE, maxPixels=1e11, tileScale=16).getArray("p").slice(1, 1, 2).project([0]).toList().map(lambda v: ee.Number(v).format("%.2f"))).join(",")})
    if which in ("all", "tibet"):
        ee.batch.Export.table.toAsset(collection=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([85, 30]), props)]), description="tp_table_tibet_hist", assetId=TIBET).start()

    if which == "tibet":
        return
    subs = ee.FeatureCollection(A + "tp_regions_v2_subbasin").map(lambda f: f.select(["unit", "parent", "km2"]))

    def year(y):
        """Snow duration of one year over the plateau ranges."""
        i = ee.ImageCollection([ee.Image(f"{ROOT}/{r}_{y}") for r in RANGES]).mosaic()
        dur = i.select("snow_days").multiply(365).divide(i.select("valid_days")).updateMask(i.select("valid_days").gte(MIN_VALID)).rename("duration")
        return dur.addBands(i.select("valid_days").toFloat()).reduceRegions(subs, ee.Reducer.mean(), 1000, tileScale=16).map(lambda f: ee.Feature(ee.Geometry.Point([85, 30]), f.toDictionary()).set("year", y))

    fc = ee.FeatureCollection([year(y) for y in YEARS]).flatten()
    ee.batch.Export.table.toAsset(collection=fc, description="tp_table_subbasin_year", assetId=SUB).start()
    print("started")


def read():
    """Read the exported histogram and totals."""
    f = ee.FeatureCollection(TIBET).first().getInfo()["properties"]
    h, total = [float(x) for x in f["hist"].split(",")], f["count"]
    cum, thr = 0, 0.0
    for k, c in enumerate(h):
        cum += c
        edge = (k + 1) * TOP / BINS
        if cum and edge <= Q * cum / total:
            thr = edge
    pd.DataFrame([dict(pixels=total, below_002=sum(h), threshold=thr)]).to_csv(R / "outputs/tibet_fdr_500m.csv", index=False)
    print("Tibet at 500 m: pixels", round(total), "threshold", thr, "share below 0.002", round(sum(h) / total, 5))

    rows = [x["properties"] for x in ee.FeatureCollection(SUB).getInfo()["features"]]
    d = pd.DataFrame(rows)
    d = d.rename(columns={"mean": "duration"}) if "mean" in d and "duration" not in d else d
    d.to_csv(R / "outputs/subbasin_year.csv", index=False)
    out = []
    for (u, par), g in d.groupby(["unit", "parent"]):
        v = g.sort_values("year").duration.dropna().values
        if len(v) < 20:
            out.append(dict(unit=u, parent=par, km2=g.km2.iloc[0], n=len(v)))
            continue
        a, b = mk.original_test(v), mk.trend_free_pre_whitening_modification_test(v)
        out.append(dict(unit=u, parent=par, km2=g.km2.iloc[0], n=len(v), mean_days=v.mean(), slope_days_decade=a.slope * 10, tau=a.Tau, p=max(a.p, b.p)))
    t = pd.DataFrame(out)
    ok = t.p.notna()
    t.loc[ok, "q"] = multipletests(t.p[ok], method="fdr_bh")[1]
    t.to_csv(R / "outputs/subbasin_trends.csv", index=False)
    print(len(t), "sub-basins;", int(ok.sum()), "tested;", int((t.p < 0.05).sum()), "p<0.05;", int((t.q < 0.05).sum()), "q<0.05")
    print(t[t.p < 0.05].sort_values("p")[["unit", "km2", "mean_days", "slope_days_decade", "p", "q"]].round(3).head(25).to_string())
    print(t.groupby("parent").agg(n=("unit", "size"), km2=("km2", "sum"), median_slope=("slope_days_decade", "median")).round(1).to_string())


if __name__ == "__main__":
    {"start": start, "start-subbasins": lambda: start("sub"), "start-tibet": lambda: start("tibet"), "read": read}[sys.argv[1]]()
