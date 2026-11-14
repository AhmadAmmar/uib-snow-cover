"""Check daily MOD10A1 snow against same-day Landsat 8/9, by elevation band and season.

For each range and season, up to SCENES Landsat Collection 2 Level 2 scenes (2013-2025, under 20% cloud) are drawn at
random (fixed seed). Landsat snow: NDSI = (green - SWIR1)/(green + SWIR1) > 0.4 and green reflectance > 0.1, with cloud,
cloud shadow and fill removed by QA_PIXEL. The 30 m snow mask is averaged to the 500 m MODIS grid (fraction of snow in
each MODIS cell, cells at least 90% clear). MODIS snow: NDSI_Snow_Cover >= 40, as in src/11, same date, clear only.
POINTS cells per scene are sampled with their elevation (Copernicus GLO-30, 2024_1).

Output: outputs/landsat_validation_points.csv (one row per sampled cell); summary by src/30.
"""
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
OUT = R / "outputs" / "landsat_validation_points.csv"
REG = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range")
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
SEASONS = {"winter": [12, 1, 2], "spring": [3, 4, 5], "summer": [6, 7, 8], "autumn": [9, 10, 11]}
SCENES, POINTS = 8, 300
MODIS = ee.ImageCollection("MODIS/061/MOD10A1").select("NDSI_Snow_Cover")
DEM = ee.ImageCollection("COPERNICUS/DEM/GLO30_2024_1").select("DEM").mosaic().setDefaultProjection("EPSG:4326", None, 30)
PROJ = MODIS.first().projection()
LS = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))


def landsat_fraction(img):
    qa = img.select("QA_PIXEL")
    clear = qa.bitwiseAnd(0b11111).eq(0)  # fill, dilated cloud, cirrus, cloud, shadow
    g = img.select("SR_B3").multiply(0.0000275).add(-0.2)
    s = img.select("SR_B6").multiply(0.0000275).add(-0.2)
    snow = g.subtract(s).divide(g.add(s)).gt(0.4).And(g.gt(0.1)).updateMask(clear)
    frac = snow.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(PROJ).rename("ls_frac")
    cov = clear.unmask(0).reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(PROJ).rename("ls_clear")
    return frac.addBands(cov)


def scene_rows(img, r, season):
    date = ee.Date(img.get("system:time_start"))
    mod = MODIS.filterDate(date.update(hour=0, minute=0, second=0), date.update(hour=0, minute=0, second=0).advance(1, "day")).first()
    stack = landsat_fraction(img).addBands(ee.Image(mod).rename("modis")).addBands(
        DEM.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(PROJ).rename("elev"))
    stack = stack.updateMask(stack.select("ls_clear").gte(0.9))
    pts = stack.sample(region=img.geometry(), projection=PROJ, numPixels=POINTS, seed=1, dropNulls=True, tileScale=4)
    d = pd.DataFrame([f["properties"] for f in pts.getInfo()["features"]])
    if d.empty:
        return d
    return d.assign(range=r, season=season, scene=img.get("system:index").getInfo(), date=date.format("YYYY-MM-dd").getInfo())


def main():
    have = set(pd.read_csv(OUT)[["range", "season"]].apply(tuple, axis=1)) if OUT.exists() else set()
    for r in RANGES:
        geom = REG.filter(ee.Filter.eq("unit", r)).filter(ee.Filter.eq("level", "range")).geometry()
        for season, months in SEASONS.items():
            if (r, season) in have:
                continue
            col = (LS.filterBounds(geom).filterDate("2013-04-01", "2025-10-01").filter(ee.Filter.lt("CLOUD_COVER", 20))
                   .filter(ee.Filter.calendarRange(months[0], months[-1], "month") if months[0] < months[-1]
                           else ee.Filter.Or(ee.Filter.calendarRange(12, 12, "month"), ee.Filter.calendarRange(1, 2, "month")))
                   .randomColumn("rnd", 7).sort("rnd").limit(SCENES))
            parts = []
            for k in range(col.size().getInfo()):
                img = ee.Image(col.toList(SCENES).get(k))
                try:
                    parts.append(scene_rows(img, r, season))
                except ee.EEException as e:
                    print("skip", r, season, k, str(e)[:80], flush=True)
            d = pd.concat(parts) if parts else pd.DataFrame()
            if not d.empty:
                d.to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
            print(r, season, len(d), flush=True)


if __name__ == "__main__":
    main()
