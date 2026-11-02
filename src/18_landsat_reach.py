"""How far back can Landsat support the work? Scenes per year, and harmonisation checks, for test areas in the Third Pole and Kashmir.

Writes outputs/landsat_availability.csv and outputs/landsat_cross_sensor.csv.
"""
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
from research_helpers import landsat

R = Path(__file__).resolve().parent.parent
AREAS = {"Muzaffarabad": [73.3, 34.25, 73.65, 34.5], "Hunza": [74.3, 36.0, 75.5, 37.0], "Baltoro": [75.8, 35.5, 76.8, 36.0], "Khumbu": [86.6, 27.8, 87.0, 28.1], "Hengduan": [99.0, 28.0, 100.0, 29.0], "Sindh": [67.3, 26.0, 68.6, 27.4]}
rows, checks = [], []
for name, box in AREAS.items():
    region = ee.Geometry.Rectangle(box)
    for r in landsat.availability(region):
        rows.append(dict(area=name, **r))
    c = landsat.cross_sensor_check(region)
    for band in landsat.COMMON:
        checks.append(dict(area=name, band=band, before=c["before"][band], after=c["after"][band]))
    print(name, flush=True)
d = pd.DataFrame(rows)
d.to_csv(R / "outputs/landsat_availability.csv", index=False)
k = pd.DataFrame(checks)
k.to_csv(R / "outputs/landsat_cross_sensor.csv", index=False)
d["group"] = d.sensor.str[:2].map({"LM": "MSS", "LT": "TM", "LE": "ETM+", "LC": "OLI"})
d["period"] = pd.cut(d.year, [1971, 1979, 1984, 1989, 1994, 1999, 2004, 2012, 2025], labels=["1972-79", "80-84", "85-89", "90-94", "95-99", "2000-04", "05-12", "13-25"])
print(d.pivot_table(index="area", columns="period", values="scenes", aggfunc="sum", observed=False).fillna(0).astype(int).to_string())
print(d.pivot_table(index="group", columns="period", values="scenes", aggfunc="sum", observed=False).fillna(0).astype(int).to_string())
print(k.groupby("band")[["before", "after"]].agg(lambda s: round(s.abs().mean(), 4)).to_string())
