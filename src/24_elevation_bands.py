"""Snow-cover duration by 500 m elevation band, per range and hydrological year, from the yearly assets of src/11.

Duration is defined as in src/13: snow days scaled to 365 by clear days, pixels with at least 120 clear days.
Elevation is the Copernicus GLO-30 DEM (2024_1 release), averaged to the 500 m MODIS grid.
Output: outputs/band_year.csv (range, band_low_m, year, duration, pixels).
Runs at 1 km scale with tileScale 8, like src/13. Resumes: years already in the file are skipped.
"""
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
OUT = R / "outputs" / "band_year.csv"
ROOT = "projects/ee-ahmadammargeo/assets/tp_snow_v1"
REG = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range")
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
YEARS = range(2001, 2026)
STEP = 500
MIN_VALID = 120
DEM = ee.ImageCollection("COPERNICUS/DEM/GLO30_2024_1").select("DEM").mosaic().setDefaultProjection("EPSG:4326", None, 30)


def key(n):
    return n.lower().replace(" ", "_")


def rows(r, y, geom):
    i = ee.Image(f"{ROOT}/{key(r)}_{y}")
    v = i.select("valid_days")
    dur = i.select("snow_days").multiply(365).divide(v).updateMask(v.gte(MIN_VALID)).rename("duration")
    band = DEM.reduceResolution(ee.Reducer.mean(), maxPixels=1024).reproject(i.projection()).divide(STEP).floor().int().rename("band")
    red = ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True).group(groupField=1, groupName="band")
    g = dur.addBands(band).reduceRegion(red, geom, 1000, maxPixels=1e12, tileScale=8).get("groups").getInfo()
    return [dict(range=r, band_low_m=int(x["band"]) * STEP, year=y, duration=x.get("mean"), pixels=x.get("count", 0)) for x in g if x.get("mean") is not None]


def main():
    done = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["range", "year"])
    have = set(zip(done["range"], done["year"]))
    for r in RANGES:
        geom = REG.filter(ee.Filter.eq("unit", r)).filter(ee.Filter.eq("level", "range")).geometry()
        for y in YEARS:
            if (r, y) in have:
                continue
            new = pd.DataFrame(rows(r, y, geom))
            new.to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
            print(r, y, len(new), flush=True)


if __name__ == "__main__":
    main()
