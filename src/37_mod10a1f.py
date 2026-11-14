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
  python src/37_mod10a1f.py timing [first last]; python src/37_mod10a1f.py timing_series   -> snow onset and melt-out
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


def reduce_raw(first=2001, last=2025):
    """The same counts from the unfilled MOD10A1 layer inside the same files: same pixels and days, cloud left out."""
    from pyhdf.SD import SD, SDC
    for y in range(first, last + 1):
        out = D / f"counts_raw_{y}.npz"
        files = sorted((D / str(y)).glob("*.hdf"))
        if out.exists() or not files:
            continue
        snow = valid = None
        for f in files:
            a = SD(str(f), SDC.READ).select("MOD10A1_NDSI_Snow_Cover")[:].astype(np.int16)
            if snow is None:
                snow = np.zeros(a.shape, np.int16)
                valid = np.zeros(a.shape, np.int16)
            ok = a <= 100
            snow += ok & (a >= 40)
            valid += ok
        np.savez_compressed(out, snow=snow, valid=valid, days=len(files))
        print("raw", y, len(files), flush=True)


def compare_same_pixels(unit="Karakoram", min_valid=120):
    """Duration on identical pixels: gap-filled versus clear days only (scaled to 365, at least min_valid clear days)."""
    import importlib.util
    from scipy import stats
    spec = importlib.util.spec_from_file_location("det", R / "src" / "22_detectability.py")
    det = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(det)
    m = tile_mask(unit)
    rows = []
    for y in range(2001, 2026):
        f, g = D / f"counts_{y}.npz", D / f"counts_raw_{y}.npz"
        if not (f.exists() and g.exists()):
            continue
        a, b = np.load(f), np.load(g)
        ok = m & (b["valid"] >= min_valid) & (a["valid"] >= 300)
        rows.append(dict(year=y, filled=float((a["snow"][ok].astype("float64") * 365 / a["valid"][ok]).mean()),
                         clear_only=float((b["snow"][ok].astype("float64") * 365 / b["valid"][ok]).mean()),
                         clear_days=float(b["valid"][ok].mean()), cells=int(ok.sum())))
    t = pd.DataFrame(rows)
    t.to_csv(R / "outputs/cloud_treatment_same_pixels.csv", index=False)
    print(t.round(1).to_string(index=False))
    print("r =", round(t.filled.corr(t.clear_only), 3), "| mean difference", round((t.filled - t.clear_only).mean(), 1), "days")
    for col in ("filled", "clear_only"):
        x, v = t.year.to_numpy(float), t[col].to_numpy()
        print(col, "Sen", round(10 * det.sen(x, v), 2), "days/decade, p", round(2 * stats.norm.sf(abs(det.mk_z(v[None])[0])), 3))


def timing(first=2001, last=2025, run=5):
    """Snow onset and melt-out per pixel, in days from 1 October, read in one pass with per-pixel counters.
    Onset: first day of the first spell of at least `run` consecutive snow days. Melt-out: last day of the last such
    spell. Pixels with no such spell get -1. Saved as data/mod10a1f/timing_<year>.npz."""
    from pyhdf.SD import SD, SDC
    for y in range(first, last + 1):
        out = D / f"timing_{y}.npz"
        files = sorted((D / str(y)).glob("*.hdf"))
        if out.exists() or len(files) < 360:
            continue
        streak = onset = melt = None
        for day, f in enumerate(files):
            a = SD(str(f), SDC.READ).select("CGF_NDSI_Snow_Cover")[:]
            snow = (a <= 100) & (a >= 40)
            if streak is None:
                streak = np.zeros(a.shape, np.int16)
                onset = np.full(a.shape, -1, np.int16)
                melt = np.full(a.shape, -1, np.int16)
            streak = np.where(snow, streak + 1, 0).astype(np.int16)
            hit = streak == run
            onset[hit & (onset < 0)] = day - run + 1
            melt[streak >= run] = day
        np.savez_compressed(out, onset=onset, melt=melt, days=len(files))
        print("timing", y, flush=True)


def timing_series(unit="Karakoram"):
    """Mean onset and melt-out over the unit's pixels that have a snow season every year; Sen slope and Mann-Kendall p."""
    import importlib.util
    from scipy import stats
    spec = importlib.util.spec_from_file_location("det", R / "src" / "22_detectability.py")
    det = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(det)
    m = tile_mask(unit)
    z = {int(f.stem.split("_")[1]): np.load(f) for f in sorted(D.glob("timing_*.npz"))}
    every = m.copy()
    for v in z.values():
        every &= (v["onset"] >= 0)
    rows = [dict(year=y, onset_day=float(v["onset"][every].mean()), meltout_day=float(v["melt"][every].mean()),
                 season_days=float((v["melt"][every] - v["onset"][every]).mean()), pixels=int(every.sum()))
            for y, v in z.items()]
    t = pd.DataFrame(rows)
    t.to_csv(R / f"outputs/timing_{unit.lower().replace(' ', '_')}.csv", index=False)
    for col in ("onset_day", "meltout_day", "season_days"):
        x, yv = t.year.to_numpy(float), t[col].to_numpy()
        print(col, "Sen", round(10 * det.sen(x, yv), 2), "days/decade, p", round(2 * stats.norm.sf(abs(det.mk_z(yv[None])[0])), 3))
    print(t.round(1).to_string(index=False))


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
        rows.append(dict(year=y, duration_cgf=float((z["snow"][ok].astype("float64") * 365 / z["valid"][ok]).mean()), cells=int(ok.sum()),
                         share_of_unit_in_tile=float(m.sum())))
    t = pd.DataFrame(rows)
    u = pd.read_csv(R / "outputs/unit_year.csv").query("unit == 'Karakoram' and level == 'range'")[["year", "duration"]]
    t = t.merge(u, on="year", how="left")
    t.to_csv(R / "outputs/mod10a1f_year.csv", index=False)
    print(t.round(2).to_string(index=False))
    print("r =", round(t.duration_cgf.corr(t.duration), 3))


if __name__ == "__main__":
    cmd, *rest = sys.argv[1:] or ["compare"]
    {"download": download, "reduce": reduce, "reduce_raw": reduce_raw, "compare": compare, "compare_same_pixels": compare_same_pixels,
     "timing": timing, "timing_series": timing_series}[cmd](*map(int, rest))
