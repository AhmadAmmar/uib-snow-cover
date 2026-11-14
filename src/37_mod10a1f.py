"""Cloud-gap-filled daily snow, MOD10A1F v6.1 (NSIDC; doi:10.5067/MODIS/MOD10A1F.061), for one MODIS tile, h24v05,
which covers the Karakoram and most of the Upper Indus. Tests the cloud treatment: our MOD10A1 duration leaves cloudy
days out and scales by clear days; MOD10A1F fills them from the last clear view.

Needs NASA Earthdata login stored once with earthaccess (persist=True writes ~/_netrc). About 3 MB a day, ~28 GB for
2001-2025, saved to data/mod10a1f (outside git). Each file is read once and reduced to yearly per-pixel counts, saved as
data/mod10a1f/counts_<year>.npz (snow days at CGF_NDSI_Snow_Cover >= 40, valid days), so the HDF files can be deleted
after a year is reduced if space is short; nothing else in the project is touched.

  python src/37_mod10a1f.py download [first last]
  python src/37_mod10a1f.py reduce   [first last]
  python src/37_mod10a1f.py compare     -> outputs/mod10a1f_year.csv: Karakoram duration, gap-filled vs our MOD10A1
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent.parent
D = R / "data" / "mod10a1f"
TILE = "h24v05"


def hydro_days(y):
    return pd.Timestamp(f"{y - 1}-10-01"), pd.Timestamp(f"{y}-09-30")


def download(first=2001, last=2025):
    import earthaccess
    earthaccess.login(strategy="netrc")
    D.mkdir(parents=True, exist_ok=True)
    for y in range(first, last + 1):
        a, b = hydro_days(y)
        g = earthaccess.search_data(short_name="MOD10A1F", version="61", granule_name=f"*{TILE}*",
                                    temporal=(str(a.date()), str(b.date())))
        earthaccess.download(g, str(D / str(y)))
        print(y, len(g), flush=True)


def reduce(first=2001, last=2025):
    from pyhdf.SD import SD, SDC
    for y in range(first, last + 1):
        out = D / f"counts_{y}.npz"
        files = sorted((D / str(y)).glob("*.hdf"))
        if out.exists() or not files:
            continue
        snow = valid = None
        for f in files:
            a = SD(str(f), SDC.READ).select("CGF_NDSI_Snow_Cover")[:].astype(np.int16)
            if snow is None:
                snow = np.zeros(a.shape, np.int16)
                valid = np.zeros(a.shape, np.int16)
            ok = a <= 100
            snow += ok & (a >= 40)
            valid += ok
        np.savez_compressed(out, snow=snow, valid=valid, days=len(files))
        print(y, len(files), flush=True)


def tile_mask(unit="Karakoram"):
    """Mask of the tile's 2400 x 2400 sinusoidal cells whose centres fall in the unit."""
    import geopandas as gpd
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    size = 463.312716528
    h, v = 24, 5
    x0 = -20015109.354 + h * 2400 * size
    y0 = 10007554.677 - v * 2400 * size
    sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
    reg = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    geom = reg[(reg.unit == unit) & (reg.level == "range")].to_crs(sinu).geometry
    return rasterize([(g, 1) for g in geom], out_shape=(2400, 2400), transform=from_origin(x0, y0, size, size)).astype(bool)


def compare():
    m = tile_mask()
    rows = []
    for f in sorted(D.glob("counts_*.npz")):
        z = np.load(f)
        y = int(f.stem.split("_")[1])
        ok = m & (z["valid"] >= 300)
        rows.append(dict(year=y, duration_cgf=float((z["snow"][ok] * 365 / z["valid"][ok]).mean()), cells=int(ok.sum()),
                         share_of_unit_in_tile=float(m.sum())))
    t = pd.DataFrame(rows)
    u = pd.read_csv(R / "outputs/unit_year.csv").query("unit == 'Karakoram' and level == 'range'")[["year", "duration"]]
    t = t.merge(u, on="year", how="left")
    t.to_csv(R / "outputs/mod10a1f_year.csv", index=False)
    print(t.round(2).to_string(index=False))
    print("r =", round(t.duration_cgf.corr(t.duration), 3))


if __name__ == "__main__":
    cmd, *rest = sys.argv[1:] or ["compare"]
    {"download": download, "reduce": reduce, "compare": compare}[cmd](*map(int, rest))
