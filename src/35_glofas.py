"""Daily river discharge for the Upper Indus from the GloFAS v4 historical reanalysis (Copernicus Emergency Management
Service, via the Early Warning Data Store, ewds.climate.copernicus.eu: dataset cems-glofas-historical; licence
accepted on the dataset page first). The same token as in ~/.cdsapirc is used.

Needs ~/.cdsapirc with the user's own token (the CDS and EWDS share one ECMWF account).
Downloads one year per request over a small box around the gauges, to data/glofas/glofas_<year>.nc, then picks for each
gauge the cell within 0.15 degrees with the largest mean discharge (the main channel; GloFAS cells are 0.05 degrees).
Gauges (approximate published positions; check against the GloFAS upstream-area map before relying on one):
  Indus at Besham Qila (inflow to Tarbela)  34.92 N, 72.87 E
  Hunza at Dainyor bridge                   35.93 N, 74.37 E
  Gilgit at Gilgit                          35.92 N, 74.30 E
  Astore at Doyian                          35.55 N, 74.70 E
  Shyok at Yugo                             35.18 N, 76.10 E
Output: outputs/glofas_daily.csv (gauge, date, discharge m3/s), outputs/glofas_year.csv (hydrological-year means and
warm-season, May-September, means).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent.parent
D = R / "data" / "glofas"
# Only the main stem: "largest mean flow near the gauge" finds the Indus correctly, but near a tributary it finds the
# Indus instead (tested 2026-11-14: Astore came out at 1,469 m3/s, about ten times its real flow). Tributary gauges need
# the GloFAS upstream-area map matched to each gauge's published catchment area before they are added back.
GAUGES = {"Indus at Besham Qila": (34.92, 72.87)}
AREA = [36.5, 72.4, 34.6, 76.5]  # N, W, S, E


def download(first=2000, last=2025):
    import cdsapi
    D.mkdir(parents=True, exist_ok=True)
    key = [ln.split(":", 1)[1].strip() for ln in (Path.home() / ".cdsapirc").read_text().splitlines() if ln.startswith("key")][0]
    c = cdsapi.Client(url="https://ewds.climate.copernicus.eu/api", key=key)  # GloFAS moved to the Early Warning Data Store
    for y in range(first, last + 1):
        out = D / f"glofas_{y}.nc"
        if out.exists():
            continue
        c.retrieve("cems-glofas-historical", {
            "system_version": ["version_4_0"], "hydrological_model": ["lisflood"], "product_type": ["consolidated"],
            "variable": ["average_river_discharge_in_the_last_24_hours"], "timespan": ["time_mean"], "year": [str(y)],
            "month": [f"{m:02d}" for m in range(1, 13)], "day": [f"{d:02d}" for d in range(1, 32)],
            "data_format": "netcdf", "download_format": "unarchived", "area": AREA}, str(out))
        print("downloaded", y, flush=True)


def extract():
    import xarray as xr
    # one year at a time, joined in memory: the files are small, and this needs no dask
    parts = []
    for f in sorted(D.glob("glofas_*.nc")):
        with xr.open_dataset(f) as one:
            parts.append(one[[v for v in one.data_vars if one[v].ndim == 3][0]].load())
    q = xr.concat(parts, dim=[c for c in parts[0].dims if "time" in c][0])
    tname = [c for c in q.dims if "time" in c][0]
    mean = q.mean(tname).compute()
    rows = []
    for g, (la, lo) in GAUGES.items():
        box = mean.sel(latitude=slice(la + 0.15, la - 0.15), longitude=slice(lo - 0.15, lo + 0.15))
        k = np.unravel_index(int(np.nanargmax(box.values)), box.shape)
        cla, clo = float(box.latitude[k[0]]), float(box.longitude[k[1]])
        s = q.sel(latitude=cla, longitude=clo).to_series()
        # GloFAS stamps each 24-hour mean at the end of the day it covers: date it to that day
        rows.append(pd.DataFrame({"gauge": g, "date": s.index - pd.Timedelta(days=1), "discharge": s.values, "cell_lat": cla, "cell_lon": clo}))
    d = pd.concat(rows, ignore_index=True)
    d.to_csv(R / "outputs/glofas_daily.csv", index=False)
    d["date"] = pd.to_datetime(d.date)
    d["hy"] = d.date.dt.year + (d.date.dt.month >= 10)
    warm = d[d.date.dt.month.between(5, 9)].groupby(["gauge", d.date.dt.year]).discharge.mean().rename("warm_mean")
    year = d.groupby(["gauge", "hy"]).discharge.mean().rename("hy_mean")
    t = pd.concat([year, warm.rename_axis(["gauge", "hy"])], axis=1).reset_index().rename(columns={"hy": "year"})
    t.to_csv(R / "outputs/glofas_year.csv", index=False)
    print(t.groupby("gauge")[["hy_mean", "warm_mean"]].mean().round(0))


if __name__ == "__main__":
    (download if sys.argv[1:] == ["download"] else extract)()
