"""Yearly snow duration and climate per range and division, for every range whose yearly assets are complete.

Saves outputs/unit_year.csv (duration, late-season fraction, clear days, ERA5-Land climate) and data/change_<range>.tif comes from src/12.
"""
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
ROOT = "projects/ee-ahmadammargeo/assets/tp_snow_v1"
REG = "projects/ee-ahmadammargeo/assets/tp_regions_v2_"
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
YEARS = list(range(2001, 2026))
MIN_VALID = 120
ERA = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")


def key(n):
    """Asset key of a range name."""
    return n.lower().replace(" ", "_")


def complete():
    """Ranges whose yearly assets all exist."""
    have = {a["name"].split("/")[-1] for a in ee.data.listAssets({"parent": ROOT})["assets"]}
    return [r for r in RANGES if all(f"{key(r)}_{y}" in have for y in YEARS)]


def snow_image(r, y):
    """Snow duration (days per year, scaled to 365 from valid days, masked below MIN_VALID) late-summer snow fraction and valid days of a range-year."""
    i = ee.Image(f"{ROOT}/{key(r)}_{y}")
    dur = i.select("snow_days").multiply(365).divide(i.select("valid_days")).updateMask(i.select("valid_days").gte(MIN_VALID)).rename("duration")
    late = i.select("late_snow").divide(i.select("late_valid")).rename("late_frac")
    return dur.addBands(late).addBands(i.select("valid_days").toFloat())


def climate_image(y):
    """ERA5-Land climate of a hydrological year: mean 2 m temperature over the year, cold
    season (October-April) and warm season (May-September), precipitation and cold-season snowfall (mm), and more."""
    wy = ERA.filterDate(f"{y - 1}-10-01", f"{y}-10-01")
    cold = ERA.filterDate(f"{y - 1}-10-01", f"{y}-05-01")
    warm = ERA.filterDate(f"{y}-05-01", f"{y}-10-01")
    return ee.Image.cat(
        wy.select("temperature_2m").mean().subtract(273.15).rename("t2m_wy"),
        cold.select("temperature_2m").mean().subtract(273.15).rename("t2m_cold"),
        warm.select("temperature_2m").mean().subtract(273.15).rename("t2m_warm"),
        wy.select("total_precipitation_sum").sum().multiply(1000).rename("precip_wy"),
        cold.select("snowfall_sum").sum().multiply(1000).rename("snowfall_cold"),
    )


def table(img, units, scale, y):
    """Mean of an image over each unit as a table, with the year."""
    fc = img.reduceRegions(units, ee.Reducer.mean(), scale, tileScale=8).map(lambda f: f.setGeometry(None).set("year", y))
    return pd.DataFrame([f["properties"] for f in fc.getInfo()["features"]])


def main():
    """Unit-year table of snow and climate for every complete range."""
    done = complete()
    print("complete ranges:", done)
    rows = []
    for r in done:
        ra = ee.FeatureCollection(REG + "range").filter(ee.Filter.eq("unit", r))
        di = ee.FeatureCollection(REG + "division").filter(ee.Filter.eq("parent", r))
        units = ra.merge(di)
        for y in YEARS:
            a = table(snow_image(r, y), units, 1000, y)
            b = table(climate_image(y), units, 11132, y)
            rows.append(a.merge(b[["unit", "level", "t2m_wy", "t2m_cold", "t2m_warm", "precip_wy", "snowfall_cold"]], on=["unit", "level"]))
            print(r, y, flush=True)
    t = pd.concat(rows)[["unit", "level", "parent", "km2", "year", "duration", "late_frac", "valid_days", "t2m_wy", "t2m_cold", "t2m_warm", "precip_wy", "snowfall_cold"]]
    (R / "outputs").mkdir(exist_ok=True)
    t.to_csv(R / "outputs/unit_year.csv", index=False)
    print(len(t), "rows")


if __name__ == "__main__":
    main()
