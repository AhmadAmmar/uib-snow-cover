"""Check whether the coarse p-value bins in src/12 hide pixels that pass false discovery control, using bins of one millionth.

Writes outputs/fdr_check.csv with the threshold and signed area for each range. Pixels are counted at 1 km to stay within limits.
"""
import importlib.util
from pathlib import Path

import ee
import pandas as pd

R = Path(__file__).resolve().parent.parent
sp = importlib.util.spec_from_file_location("c", R / "src/12_change_classes.py")
C = importlib.util.module_from_spec(sp)
sp.loader.exec_module(C)
SCALE, TOP, BINS, Q = 1000, 0.002, 2000, 0.05
rows = []
for name in ["Karakoram", "Pamir Mountains", "Hindu Kush", "Himalaya", "Tian Shan", "Tibetan Plateau"]:
    key = name.lower().replace(" ", "_")
    rng = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range").filter(ee.Filter.eq("unit", name)).geometry()
    stored = f"{C.ROOT}_mk/{key}"
    try:
        ee.data.getAsset(stored)
        m = ee.Image(stored)
        s, p, n = m.select("s"), m.select("p"), m.select("n")
    except ee.EEException:
        dur, _ = C.stack(key)
        s, _z, p, n = C.mann_kendall(dur)
    p = p.updateMask(n.gte(C.MIN_YEARS)).rename("p")
    total = p.reduceRegion(ee.Reducer.count(), rng, SCALE, maxPixels=1e10, tileScale=8).getInfo()["p"]
    hist = p.updateMask(p.lt(TOP)).reduceRegion(ee.Reducer.fixedHistogram(0, TOP, BINS), rng, SCALE, maxPixels=1e10, tileScale=8).getInfo()["p"]
    cum, thr = 0, 0.0
    for lo, c in hist:
        cum += c
        edge = lo + TOP / BINS
        if cum and edge <= Q * cum / total:
            thr = edge
    area = ee.Image.pixelArea().divide(1e6)
    sig = p.lte(thr) if thr else p.lt(0)
    gain = area.updateMask(sig.And(s.gt(0))).reduceRegion(ee.Reducer.sum(), rng, SCALE, maxPixels=1e10, tileScale=8).getInfo()["area"] or 0
    loss = area.updateMask(sig.And(s.lt(0))).reduceRegion(ee.Reducer.sum(), rng, SCALE, maxPixels=1e10, tileScale=8).getInfo()["area"] or 0
    low = sum(c for lo, c in hist)
    rows.append(dict(range=name, pixels=total, share_p_below_002=low / total, threshold=thr, gain_km2=gain, loss_km2=loss))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(R / "outputs/fdr_check.csv", index=False)
