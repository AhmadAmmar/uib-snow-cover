"""Fit sensor harmonisation from near-same-time image pairs, instead of medians with different dates.

Pairs lie within MAX_DAYS. Each pixel is clear in both scenes, and half of the pairs, chosen by scene, are held back for testing.
The fit is reduced major axis, as neither sensor is free of error, and it is compared with no change and with Roy et al. (2016).
Writes outputs/harmonisation_fit.csv.
"""
import sys
from pathlib import Path

import ee
import numpy as np
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
from research_helpers import landsat as L

R = Path(__file__).resolve().parent.parent
AREAS = {"Muzaffarabad": [73.3, 34.25, 73.65, 34.5], "Hunza": [74.3, 36.0, 75.5, 37.0], "Khumbu": [86.6, 27.8, 87.0, 28.1], "Hengduan": [99.0, 28.0, 100.0, 29.0], "Sindh": [67.3, 26.0, 68.6, 27.4], "TianShan": [78.0, 41.5, 79.5, 42.5]}
PAIRS = {"ETM_to_OLI": ("LANDSAT/LE07/C02/T1_L2", "LANDSAT/LC08/C02/T1_L2", "2013-04-01", "2018-01-01"), "TM_to_ETM": ("LANDSAT/LT05/C02/T1_L2", "LANDSAT/LE07/C02/T1_L2", "2000-01-01", "2011-11-01")}
MAX_DAYS, MAX_PAIRS, PIXELS = 1, 80, 300
ROOT = "projects/ee-ahmadammargeo/assets/tp_harm_"


def pair_table(a_asset, b_asset, start, end, region, tag):
    a = L.prepare(a_asset, region, start, end, harmonise=False)
    b = L.prepare(b_asset, region, start, end, harmonise=False)
    window = ee.Filter.And(ee.Filter.maxDifference(MAX_DAYS * 86400000, leftField="system:time_start", rightField="system:time_start"), ee.Filter.intersects(leftField=".geo", rightField=".geo"))
    pairs = ee.ImageCollection(ee.Join.saveFirst("match").apply(a, b, window)).limit(MAX_PAIRS)

    def one(im):
        im = ee.Image(im)
        both = im.rename([f"a_{x}" for x in L.COMMON]).addBands(ee.Image(im.get("match")).rename([f"b_{x}" for x in L.COMMON]))
        return both.sample(region=region, scale=60, numPixels=PIXELS, seed=1, dropNulls=True, tileScale=4, geometries=True).map(lambda f: f.set("pair", im.get("system:index"), "area", tag))
    return pairs.map(one).flatten()


def start():
    for label, (a, b, s, e) in PAIRS.items():
        fc = None
        for area, box in AREAS.items():
            x = pair_table(a, b, s, e, ee.Geometry.Rectangle(box), area)
            fc = x if fc is None else fc.merge(x)
        ee.batch.Export.table.toAsset(collection=fc, description=f"tp_harm_{label}", assetId=ROOT + label).start()
    print("started")


def rma(x, y):
    s = np.sign(np.corrcoef(x, y)[0, 1]) * y.std() / x.std()
    return s, y.mean() - s * x.mean()


def clean(d):
    """Valid reflectance, and no snow or thin cloud left, as both change quickly between dates."""
    keep = (d.filter(like="a_").min(axis=1) > 0) & (d.filter(like="b_").min(axis=1) > 0) & (d.a_blue < 0.2) & (d.b_blue < 0.2)
    ndsi = lambda s: (d[f"{s}_green"] - d[f"{s}_swir1"]) / (d[f"{s}_green"] + d[f"{s}_swir1"])
    keep &= (ndsi("a") < 0.2) & (ndsi("b") < 0.2)
    # scenes that disagree strongly overall have changed, or hold undetected cloud
    d = d[keep].copy()
    gap = (d.a_red - d.b_red).abs().groupby(d.pair).transform("median")
    return d[gap < 0.05]


def main():
    rows = []
    for label in PAIRS:
        d = pd.DataFrame([f["properties"] for f in ee.FeatureCollection(ROOT + label).limit(40000).getInfo()["features"]])
        raw = len(d)
        d = clean(d)
        keys = sorted(d.pair.unique())
        test = set(keys[::2])
        tr, te = d[~d.pair.isin(test)], d[d.pair.isin(test)]
        print(label, "pixels", raw, "after cleaning", len(d), "pairs", len(keys), d.area.value_counts().to_dict(), flush=True)
        for i, b_ in enumerate(L.COMMON):
            s, c = rma(tr[f"a_{b_}"].values, tr[f"b_{b_}"].values)
            tx, ty = te[f"a_{b_}"].values, te[f"b_{b_}"].values
            rmse = lambda p: float(np.sqrt(np.mean((p - ty) ** 2)))
            roy = tx * L.ROY_SLOPE[i] + L.ROY_INTERCEPT[i]
            rows.append(dict(pair=label, band=b_, n_train=len(tr), n_test=len(te), slope=s, intercept=c, r=np.corrcoef(tx, ty)[0, 1], rmse_none=rmse(tx), rmse_roy=rmse(roy), rmse_fit=rmse(tx * s + c), bias_none=float((tx - ty).mean()), bias_roy=float((roy - ty).mean()), bias_fit=float((tx * s + c - ty).mean())))
    t = pd.DataFrame(rows)
    t.to_csv(R / "outputs/harmonisation_fit.csv", index=False)
    print(t.round(4).to_string())


if __name__ == "__main__":
    {"start": start, "read": main}[sys.argv[1]]()
