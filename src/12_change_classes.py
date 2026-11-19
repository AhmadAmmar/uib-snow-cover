"""Classify snow change per pixel for one range from the yearly assets, and tabulate the area of each class per division.

Duration = 365 x snow_days / valid_days, only for years with at least MIN_VALID clear days. A trend needs MIN_YEARS such
years. Significance uses a Mann-Kendall normal approximation, with Benjamini-Hochberg control across pixels at Q.

Classes: 0 too little data, 1 rarely snow, 2 stable seasonal, 3 gain, 4 loss, 5 persistent snow or ice.
Glacier outlines (GLIMS) give each class an on-glacier or off-glacier area. Persistence means snow or ice is seen at the
end of the melt season; it is not a measure of glacier health.

    python src/12_change_classes.py <range> [--map]      saves outputs/change_<range>.csv, and with --map a quick-look PNG
    python src/12_change_classes.py <range> --export-mk  batch export of the trend test, for ranges too large to compute interactively
"""
import sys
from pathlib import Path

import ee
import pandas as pd
import requests

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
ROOT = "projects/ee-ahmadammargeo/assets/tp_snow_v1"
YEARS = list(range(2001, 2026))
MIN_VALID, MIN_YEARS, Q = 120, 20, 0.05
RARE_DAYS, PERSISTENT = 15, 0.8
BINS, TOP = 2000, 0.002
SCALE = 463.3127


def stack(key):
    """Bands d_<year> (duration, masked where clear days are too few), and the overall late-season snow fraction."""
    imgs = [ee.Image(f"{ROOT}/{key}_{y}") for y in YEARS]
    dur = ee.Image.cat([i.select("snow_days").multiply(365).divide(i.select("valid_days")).updateMask(i.select("valid_days").gte(MIN_VALID)).rename(f"d_{y}").toFloat() for i, y in zip(imgs, YEARS)])
    late = ee.ImageCollection(imgs)
    frac = late.select("late_snow").sum().divide(late.select("late_valid").sum()).rename("late_frac")
    return dur, frac


def mann_kendall(dur):
    """S, tau-style z and two-sided p per pixel, computed from all year pairs (no tie correction beyond a variance on n)."""
    names = [f"d_{y}" for y in YEARS]
    s = ee.Image(0).toFloat()
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            diff = dur.select(b).subtract(dur.select(a))
            s = s.add(diff.gt(0).subtract(diff.lt(0)).unmask(0))
    n = dur.reduce(ee.Reducer.count()).toFloat()
    var = n.multiply(n.subtract(1)).multiply(n.multiply(2).add(5)).divide(18)
    z = s.subtract(s.gt(0)).add(s.lt(0)).divide(var.sqrt())
    # erfc gives the two-sided normal p
    p = z.abs().divide(ee.Number(2).sqrt()).erfc()
    return s.rename("s"), z.rename("z"), p.rename("p"), n.rename("n")


def bh_threshold(p, geom, scale=SCALE):
    """Largest p threshold meeting Benjamini-Hochberg. Bins of one millionth up to TOP: a coarser histogram returned zero wrongly."""
    p = p.rename("p")
    total = p.reduceRegion(ee.Reducer.count(), geom, scale, maxPixels=1e10, tileScale=8).getInfo()["p"]
    h = p.updateMask(p.lt(TOP)).reduceRegion(ee.Reducer.fixedHistogram(0, TOP, BINS), geom, scale, maxPixels=1e10, tileScale=8).getInfo()["p"]
    cum, best = 0, 0.0
    for lo, c in h:
        cum += c
        edge = lo + TOP / BINS
        if cum and edge <= Q * cum / total:
            best = edge
    return best, total


def main(name, *flags):
    """Per-pixel change classes of snow duration in a range from a per-pixel Mann-Kendall test; --export-mk stores the test
    result as an asset."""
    key = name.lower().replace(" ", "_")
    rng = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range").filter(ee.Filter.eq("unit", name)).geometry()
    divs = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_division").filter(ee.Filter.eq("parent", name))
    dur, late = stack(key)
    s, z, p, n = mann_kendall(dur)
    stored = f"{ROOT}_mk/{key}"
    if "--export-mk" in flags:
        try:
            ee.data.createAsset({"type": "FOLDER"}, f"{ROOT}_mk")
        except ee.EEException:
            pass
        proj = ee.Image(f"{ROOT}/{key}_2001").select(0).projection()
        task = ee.batch.Export.image.toAsset(image=s.addBands(p).addBands(n).clip(rng).toFloat(), description=f"tp_mk_{key}", assetId=stored, crs=proj.crs().getInfo(), crsTransform=proj.getInfo()["transform"], region=rng.bounds(), maxPixels=1e11)
        task.start()
        print("started", task.id)
        return
    scale = SCALE
    try:
        ee.data.getAsset(stored)
        m = ee.Image(stored)
        s, p, n = m.select("s"), m.select("p"), m.select("n")
        scale = 1000  # the largest range: tabulated at 1 km to stay within interactive limits
    except ee.EEException:
        pass
    enough = n.gte(MIN_YEARS)
    p = p.updateMask(enough)
    thr, tested = bh_threshold(p, rng, scale)
    print(name, "tested pixels", int(tested), "BH p threshold", round(thr, 5))
    mean = dur.reduce(ee.Reducer.mean())
    sig = (p.lte(thr) if thr else p.lt(0)).unmask(0)
    cls = (ee.Image(0).where(enough, 2).where(enough.And(mean.lt(RARE_DAYS)), 1)
           .where(enough.And(sig).And(s.gt(0)), 3).where(enough.And(sig).And(s.lt(0)), 4)
           .where(enough.And(late.gte(PERSISTENT)), 5).rename("cls").clip(rng).toByte())
    glacier = ee.Image(0).paint(ee.FeatureCollection("GLIMS/current").filterBounds(rng), 1).rename("g")
    code = cls.multiply(10).add(glacier).rename("code")
    area = ee.Image.pixelArea().divide(1e6).addBands(code)
    reducer = ee.Reducer.sum().group(1, "code")
    rows = []
    units = [ee.Feature(rng, {"unit": name, "level": "range"})] + divs.toList(100).getInfo()
    for u in units:
        info = u if isinstance(u, dict) else None
        unit = info["properties"]["unit"] if info else name
        geom = ee.Feature(ee.Dictionary(info)).geometry() if False else (divs.filter(ee.Filter.eq("unit", unit)).geometry() if info else rng)
        got = area.reduceRegion(reducer, geom, scale, maxPixels=1e10, tileScale=8).getInfo()["groups"]
        for g in got:
            rows.append(dict(unit=unit, level="division" if info else "range", cls=int(g["code"]) // 10, on_glacier=int(g["code"]) % 10, km2=g["sum"]))
    t = pd.DataFrame(rows)
    t["cls_name"] = t.cls.map({0: "too little data", 1: "rarely snow", 2: "stable seasonal", 3: "gain", 4: "loss", 5: "persistent"})
    (R / "outputs").mkdir(exist_ok=True)
    t.to_csv(R / "outputs" / f"change_{key}.csv", index=False)
    print(t[t.level == "range"].pivot_table(index="cls_name", columns="on_glacier", values="km2", aggfunc="sum").round(0).to_string())
    # georeferenced 1 km copy for the figures, by nearest pixel
    (R / "data").mkdir(exist_ok=True)
    for res in (1000, 2000, 4000):  # the service refuses large downloads and returns an error message instead of a file
        url = cls.unmask(255).clip(rng).getDownloadURL(dict(scale=res, crs="EPSG:4326", region=rng.bounds(), format="GEO_TIFF", bands=["cls"]))
        body = requests.get(url, timeout=600).content
        if body[:2] in (b"II", b"MM"):
            (R / "data" / f"change_{key}.tif").write_bytes(body)
            print("map saved at", res, "m")
            break
    else:
        raise SystemExit("no valid map download")
    if "--map" in flags:
        (R / "figures").mkdir(exist_ok=True)
        url = cls.getThumbURL(dict(region=rng.bounds(), dimensions=1100, min=0, max=5, palette=["#bdbdbd", "#f2efe6", "#a6cee3", "#2166ac", "#b2182b", "#111111"], format="png"))
        (R / "figures" / f"change_{key}.png").write_bytes(requests.get(url, timeout=300).content)


if __name__ == "__main__":
    main(*sys.argv[1:])
