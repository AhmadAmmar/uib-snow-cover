"""Long record: NOAA Climate Data Record of Northern Hemisphere snow cover extent (Rutgers), weekly, 1966-2026
(doi:10.7289/V5N014G9; Estilow et al. 2015, doi:10.5194/essd-7-137-2015; no restrictions on use).
File: data/nhsce/nhsce_v01r01_19661004_20261102.nc, 88 x 88 polar stereographic grid, cells of roughly 190 km.

Per range: cells whose centre lies inside the range, area-weighted. Weekly snow share, then per hydrological year
(1 Oct to 30 Sep) the mean share x 365 = snow-cover duration in days, comparable in kind to the MODIS measure.
Weeks with missing cells (fill -127, 1968-1971) are dropped, then years with fewer than 48 weeks. Grid cells are coarse: ranges with fewer than 4 cells are not reported.
Outputs: outputs/nhsce_year.csv, outputs/nhsce_trends.csv (1967-2025 and 2001-2025, with the detectable trend and the
correlation with MODIS duration over 2001-2025), figures/nhsce_long.png/.pdf
"""
import importlib.util
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

R = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("det", R / "src" / "22_detectability.py")
det = importlib.util.module_from_spec(spec)
spec.loader.exec_module(det)
unit_noise = pd.read_csv(R / "outputs" / "detectability_unit_noise.csv").set_index("lag1")
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]


def mdt(n, sd, phi):
    """Minimum detectable trend for length n, noise sd and phi, from the unit noise table."""
    ns = np.array([int(c[1:]) for c in unit_noise.columns])
    per_n = [np.interp(np.clip(phi, 0, 0.9), unit_noise.index, unit_noise[c]) for c in unit_noise.columns]
    if n > ns.max():
        return det.mdt(n, sd, phi)  # direct simulation beyond the table
    return sd * np.exp(np.interp(np.log(n), np.log(ns), np.log(per_n)))


def main():
    """Snow cover from the NOAA Northern Hemisphere snow cover record (weekly snow extent) for the cells within each range."""
    d = xr.open_dataset(R / "data/nhsce/nhsce_v01r01_19661004_20261102.nc")
    lat, lon, area = d.latitude.values, d.longitude.values, d.area.values
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon.ravel(), lat.ravel()), crs="EPSG:4326")
    reg = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
    reg = reg[(reg.level == "range") & reg.unit.isin(RANGES)]
    hit = gpd.sjoin(pts, reg[["unit", "geometry"]], predicate="within")
    sce = d.snow_cover_extent.values.reshape(d.sizes["time"], -1).astype(float)
    sce[(sce < 0) | (sce > 1)] = np.nan  # -127 marks missing weeks (1968-1971)
    t = pd.to_datetime(d.time.values)
    hy = t.year + (t.month >= 10)
    rows, cells = [], {}
    for unit, h in hit.groupby("unit"):
        idx = h.index.to_numpy()
        cells[unit] = len(idx)
        w = area.ravel()[idx]
        ok = ~np.isnan(sce[:, idx]).any(1)
        share = np.nansum(sce[:, idx] * w, 1) / w.sum()
        s = pd.DataFrame({"hy": hy[ok], "share": share[ok]}).groupby("hy").share.agg(["mean", "size"])
        for y, r in s.iterrows():
            if r["size"] >= 48 and 1967 <= y <= 2025:
                rows.append(dict(unit=unit, year=int(y), duration=365 * r["mean"], weeks=int(r["size"]), cells=len(idx)))
    yt = pd.DataFrame(rows)
    yt.to_csv(R / "outputs/nhsce_year.csv", index=False)
    print("cells per range:", cells)

    modis = pd.read_csv(R / "outputs/unit_year.csv").query("level == 'range'")[["unit", "year", "duration"]]
    out = []
    for unit, g in yt.groupby("unit"):
        if cells[unit] < 4:
            continue
        for a, b in ((1967, 2025), (2001, 2025)):
            k = g[g.year.between(a, b)].sort_values("year")
            x, y = k.year.to_numpy(float), k.duration.to_numpy(float)
            s = det.sen(x, y)
            res = y - s * x
            res -= res.mean()
            m = mdt(len(y), res.std(ddof=1), np.corrcoef(res[:-1], res[1:])[0, 1])
            j = k.merge(modis, on=["unit", "year"], suffixes=("", "_modis"))
            out.append(dict(unit=unit, cells=cells[unit], period=f"{a}-{b}", n=len(y), mean_days=y.mean(),
                            sen_days_per_decade=10 * s, mk_p=2 * stats.norm.sf(abs(det.mk_z(y[None])[0])),
                            detectable_days_per_decade=10 * m,
                            r_with_modis=j.duration.corr(j.duration_modis) if len(j) > 10 else np.nan))
    tr = pd.DataFrame(out)
    tr.to_csv(R / "outputs/nhsce_trends.csv", index=False)
    print(tr.round(3).to_string(index=False))

    show = tr.unit.unique()
    fig, axes = plt.subplots(len(show), 1, figsize=(7, 1.6 * len(show) + 0.6), dpi=200, sharex=True, squeeze=False)
    for ax, unit in zip(axes[:, 0], show):
        g = yt[yt.unit == unit]
        ax.plot(g.year, g.duration, color="#555555", lw=1, label="NOAA CDR (weekly, ~190 km)")
        m = modis[modis.unit == unit]
        ax.plot(m.year, m.duration, color="#0072B2", lw=1.2, label="MODIS MOD10A1 (daily, 500 m)")
        ax.set_title(f"{unit} ({cells[unit]} cells)", fontsize=8, loc="left")
        ax.tick_params(labelsize=7)
    axes[0, 0].legend(fontsize=7, frameon=False, ncol=2)
    fig.supylabel("Snow-cover duration (days)", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(R / "figures" / f"nhsce_long.{ext}")


if __name__ == "__main__":
    main()
