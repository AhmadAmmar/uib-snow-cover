"""Is the Himalayan loss real? Follow the pixels classed as loss, by year, with other sensors, thresholds and clear-day counts.

    python src/15_cloud_tests.py mask    batch export of the loss pixels, from the Terra trend test
    python src/15_cloud_tests.py run     yearly table over those pixels; writes outputs/loss_tests.csv and outputs/loss_tests_summary.csv

Aqua starts in mid-2002, so all comparisons use 2003-2025. The pixels were chosen for loss in Terra at a threshold of 40, so
Terra at 40 is expected to fall; the tests are whether Aqua and other thresholds fall too.
"""
import importlib.util
import sys
import time
from pathlib import Path

import ee
import numpy as np
import pandas as pd
import pymannkendall as mk
from scipy import stats

R = Path(__file__).resolve().parent.parent
sp = importlib.util.spec_from_file_location("c", R / "src/12_change_classes.py")
C = importlib.util.module_from_spec(sp)
sp.loader.exec_module(C)
MASK = "projects/ee-ahmadammargeo/assets/tp_snow_v1_loss_himalaya"
THRESHOLD_P = 0.00027
YEARS = list(range(2003, 2026))


def mask():
    """Pixels of the Himalaya range used in the cloud tests."""
    rng = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range").filter(ee.Filter.eq("unit", "Himalaya")).geometry()
    dur, late = C.stack("himalaya")
    s, _z, p, n = C.mann_kendall(dur)
    m = p.lte(THRESHOLD_P).And(s.lt(0)).And(n.gte(C.MIN_YEARS)).And(late.lt(C.PERSISTENT)).selfMask().rename("loss").toByte()
    proj = ee.Image(f"{C.ROOT}/himalaya_2001").select(0).projection()
    t = ee.batch.Export.image.toAsset(image=m.clip(rng), description="tp_loss_himalaya", assetId=MASK, crs=proj.crs().getInfo(), crsTransform=proj.getInfo()["transform"], region=rng.bounds(), maxPixels=1e11)
    t.start()
    print("started", t.id)


def counts(product, start, end, thr):
    """Snow and valid days of a MODIS product over a period at an NDSI threshold."""
    col = ee.ImageCollection(product).filterDate(start, end).select("NDSI_Snow_Cover")
    snow = col.map(lambda i: i.gte(thr).unmask(0).toFloat()).sum()
    valid = col.map(lambda i: i.mask().gt(0).toFloat().rename("NDSI_Snow_Cover")).sum()
    return snow, valid


def run():
    """Compare snow duration under different cloud treatments on the same pixels (the
    snow-loss pixels of the mask)."""
    m = ee.Image(MASK)
    geom = m.geometry()
    px = m.reduceRegion(ee.Reducer.count(), geom, C.SCALE, maxPixels=1e9).getInfo()["loss"]
    print("loss pixels", px, flush=True)
    rows = []
    for y in YEARS:
        a, b = f"{y - 1}-10-01", f"{y}-10-01"
        bands = []
        for name, prod, thr in (("terra40", "MODIS/061/MOD10A1", 40), ("terra30", "MODIS/061/MOD10A1", 30), ("terra50", "MODIS/061/MOD10A1", 50), ("aqua40", "MODIS/061/MYD10A1", 40), ("aqua30", "MODIS/061/MYD10A1", 30)):
            sn, va = counts(prod, a, b, thr)
            bands += [sn.rename(f"{name}_snow"), va.rename(f"{name}_valid")]
        # gap filled: a day is clear if either sensor saw the ground, and snow if either saw snow
        t = ee.ImageCollection("MODIS/061/MOD10A1").filterDate(a, b).select("NDSI_Snow_Cover")
        q = ee.ImageCollection("MODIS/061/MYD10A1").filterDate(a, b).select("NDSI_Snow_Cover")
        joined = ee.ImageCollection(ee.Join.inner().apply(t, q, ee.Filter.equals(leftField="system:index", rightField="system:index"))).map(lambda f: ee.Image(f.get("primary")).rename("t").addBands(ee.Image(f.get("secondary")).rename("q")))
        both_snow = joined.map(lambda i: i.select("t").gte(40).unmask(0).Or(i.select("q").gte(40).unmask(0)).toFloat().rename("x")).sum().rename("merged_snow")
        both_valid = joined.map(lambda i: i.select("t").mask().gt(0).Or(i.select("q").mask().gt(0)).toFloat().rename("x")).sum().rename("merged_valid")
        im = ee.Image.cat(bands + [both_snow, both_valid]).updateMask(m)
        for attempt in range(6):  # the connection here drops at times
            try:
                r = im.reduceRegion(ee.Reducer.mean(), geom, C.SCALE, maxPixels=1e9, tileScale=8).getInfo()
                break
            except Exception as e:
                print("retry", y, str(e)[:80], flush=True)
                time.sleep(60)
        else:
            raise SystemExit(f"failed at {y}")
        r["year"] = y
        rows.append(r)
        print(y, round(r["terra40_snow"]), round(r["terra40_valid"]), round(r["aqua40_snow"]), round(r["aqua40_valid"]), flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(R / "outputs/loss_tests.csv", index=False)
    out = []
    for n in ("terra40", "terra30", "terra50", "aqua40", "aqua30", "merged"):
        dur = 365 * d[f"{n}_snow"] / d[f"{n}_valid"]
        plain = mk.original_test(dur.values)
        raw = mk.original_test(d[f"{n}_snow"].values)
        clear = mk.original_test(d[f"{n}_valid"].values)
        out.append(dict(series=n, mean_days=dur.mean(), slope_days_decade=plain.slope * 10, p=plain.p, raw_snow_days_slope=raw.slope * 10, raw_p=raw.p, clear_days_mean=d[f"{n}_valid"].mean(), clear_slope=clear.slope * 10, clear_p=clear.p, r_with_clear=stats.spearmanr(dur, d[f"{n}_valid"]).statistic))
    o = pd.DataFrame(out)
    o.to_csv(R / "outputs/loss_tests_summary.csv", index=False)
    print(o.round(3).to_string())


if __name__ == "__main__":
    {"mask": mask, "run": run}[sys.argv[1]]()
