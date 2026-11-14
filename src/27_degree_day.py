"""Degree-day snow model driven by daily ERA5-Land, as a check on what temperature and snowfall alone explain.

Each 0.1 degree cell stores snow water equivalent (mm):
  swe(t) = max(0, swe(t-1) + snowfall(t) - DDF * max(T(t), 0))
and a day counts as snow-covered while swe >= 5 mm. DDF = 4 mm per degree-day, within the 3-5 mm range fitted for
the Upper Indus (Azizi et al. 2019, doi:10.2208/jscejer.75.i_125). The run starts with no snow on 1 August before
each hydrological year, so October already holds any early snow.
Daily grids are fetched from Earth Engine a month at a time (computePixels) and the model runs locally, one year in
memory at a time. Yearly modelled snow days are averaged over each range and division, like MODIS duration in src/13.
Compare year-to-year variation, not absolute days: the 9 km grid smooths out the high ground.

Output: outputs/degree_day_year.csv (unit, level, year, model_days). Resumes by year.
"""
from pathlib import Path

import ee
import numpy as np
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
OUT = R / "outputs" / "degree_day_year.csv"
REG = "projects/ee-ahmadammargeo/assets/tp_regions_v2_"
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
DAILY = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
DDF, SWE_MIN, RES = 4.0, 5.0, 0.1
UNITS = ee.FeatureCollection(REG + "range").filter(ee.Filter.inList("unit", RANGES)).merge(
    ee.FeatureCollection(REG + "division").filter(ee.Filter.inList("parent", RANGES)))
W, S, E, N = [round(v, 1) for v in np.array(UNITS.geometry().bounds(1000).coordinates().getInfo()[0]).min(0).tolist()
              + np.array(UNITS.geometry().bounds(1000).coordinates().getInfo()[0]).max(0).tolist()]
GRID = dict(dimensions={"width": int(round((E - W) / RES)) + 1, "height": int(round((N - S) / RES)) + 1},
            affineTransform={"scaleX": RES, "shearX": 0, "translateX": W, "shearY": 0, "scaleY": -RES, "translateY": N},
            crsCode="EPSG:4326")


def pixels(img):
    a = ee.data.computePixels({"expression": img, "fileFormat": "NUMPY_NDARRAY", "grid": GRID})
    return np.stack([a[n] for n in a.dtype.names]).astype("float32")


def unit_masks():
    feats = UNITS.getInfo()["features"]
    masks = {}
    for k, f in enumerate(feats):
        fc = ee.FeatureCollection([ee.Feature(ee.Geometry(f["geometry"]))])
        m = pixels(ee.Image(0).paint(fc, 1).rename("m"))[0] > 0
        masks[(f["properties"]["unit"], f["properties"]["level"])] = m
    return masks


def year_days(y):
    swe = None
    days = None
    for start in pd.date_range(f"{y - 1}-08-01", f"{y}-09-01", freq="MS"):
        end = start + pd.offsets.MonthBegin(1)
        col = DAILY.filterDate(str(start.date()), str(end.date()))
        t = pixels(col.select("temperature_2m").toBands()) - 273.15
        sf = pixels(col.select("snowfall_sum").toBands()) * 1000
        if swe is None:
            swe = np.zeros(t.shape[1:], "float32")
            days = np.zeros(t.shape[1:], "float32")
        for d in range(t.shape[0]):
            swe = np.maximum(0, swe + np.nan_to_num(sf[d]) - DDF * np.maximum(np.nan_to_num(t[d]), 0))
            if start >= pd.Timestamp(f"{y - 1}-10-01"):
                days += swe >= SWE_MIN
    return days


def main():
    done = set(pd.read_csv(OUT).year) if OUT.exists() else set()
    masks = unit_masks()
    print("grid", GRID["dimensions"], "units", len(masks), flush=True)
    for y in range(2001, 2026):
        if y in done:
            continue
        d = year_days(y)
        rows = [dict(unit=u, level=l, year=y, model_days=float(d[m].mean())) for (u, l), m in masks.items() if m.any()]
        pd.DataFrame(rows).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
        print(y, round(np.mean([r["model_days"] for r in rows]), 1), flush=True)


if __name__ == "__main__":
    main()
