"""Snow-covered area from MOYDGL06* (Muhammad and Thapa 2020, doi:10.5194/essd-12-345-2020; data doi:10.1594/PANGAEA.901821,
CC BY 4.0), the improved Terra-Aqua 8-day product, for our ranges and divisions, 2002-2018.

Codes as in the paper: 200 and 210 snow; -200, 0, 240, 250 no snow. The 8-day files are read straight from the zip in
data/moydgl06. Yearly value per unit: mean over composites starting in the hydrological year (1 Oct to 30 Sep) of
snow pixels / all pixels in the unit. Only complete hydrological years (2003-2018) are kept.
Output: outputs/moydgl06_year.csv (unit, level, year, sca, composites).
"""
import re
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.windows
from rasterio.features import rasterize

R = Path(__file__).resolve().parent.parent
Z = zipfile.ZipFile(R / "data" / "moydgl06" / "MOYDGL06_2002_2018_HMA.zip")
OUT = R / "outputs" / "moydgl06_year.csv"
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]


def main():
    """Snow cover from the MOYDGL06 product over the ranges and their divisions."""
    names = sorted(n for n in Z.namelist() if n.endswith(".tif"))
    with rasterio.MemoryFile(Z.read(names[0])) as m, m.open() as r:
        shape, transform = r.shape, r.transform
    reg = gpd.read_file(R / "data" / "regions" / "tp_regions_v2.geojson")
    reg = reg[((reg.level == "range") & reg.unit.isin(RANGES)) | ((reg.level == "division") & reg.parent.isin(RANGES))]
    masks = {}
    for u in reg.itertuples():
        win = rasterio.windows.from_bounds(*u.geometry.bounds, transform=transform).round_offsets().round_lengths()
        win = win.intersection(rasterio.windows.Window(0, 0, shape[1], shape[0]))
        mk = rasterize([(u.geometry, 1)], out_shape=(win.height, win.width),
                       transform=rasterio.windows.transform(win, transform), dtype="uint8").astype(bool)
        masks[(u.unit, u.level)] = (win.toslices(), mk)
    rows = []
    for n in names:
        y, doy = map(int, re.search(r"_(\d{4})(\d{3})\.tif", n).groups())
        date = pd.Timestamp(y, 1, 1) + pd.Timedelta(days=doy - 1)
        hy = date.year + (date.month >= 10)
        with rasterio.MemoryFile(Z.read(n)) as m, m.open() as r:
            snow = np.isin(r.read(1), (200, 210))
        for (u, l), (sl, mk) in masks.items():
            rows.append((u, l, hy, snow[sl][mk].mean()))
        del snow
        print(n, flush=True) if doy == 1 else None
    d = pd.DataFrame(rows, columns=["unit", "level", "year", "sca"])
    t = d.groupby(["unit", "level", "year"]).sca.agg(["mean", "size"]).reset_index().rename(columns={"mean": "sca", "size": "composites"})
    t = t[t.year.between(2003, 2018)]
    t.to_csv(OUT, index=False)
    print(t.groupby("year").composites.max().to_dict())


if __name__ == "__main__":
    main()
