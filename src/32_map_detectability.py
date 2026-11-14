"""Map: per sub-basin, (a) the smallest trend 25 years can detect, as % of mean duration per decade, and
(b) the observed trend as a fraction of that size. Input: outputs/detectability.csv (src/22). Output: figures/map_detectability.*
"""
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm
from research_helpers import carto

R = Path(__file__).resolve().parent.parent
PROJ = carto.albers(85, 25, 45)
EXTENT = [60, 108, 22, 47]

g = gpd.read_file(R / "data/regions/tp_regions_v2.geojson")
d = pd.read_csv(R / "outputs/detectability.csv")
sb = g[g.level == "subbasin"].merge(d[d.level == "sub-basin"], on="unit", how="left")
ranges = g[g.level == "range"]
from shapely.geometry import Polygon, MultiPolygon


def outer(geom):
    """Exterior rings only: the plateau outline has thousands of small holes that read as black hatching."""
    parts = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
    big = [Polygon(q.exterior) for q in parts if Polygon(q.exterior).area > 0.5]
    return MultiPolygon(big) if big else geom


outlines = ranges.geometry.apply(outer)

panels = [("mdt25_pct_of_mean_per_decade", "Smallest detectable trend, 25 years\n(% of mean duration per decade)",
           [0, 5, 10, 15, 20, 30, 50, 100], "viridis_r"),
          ("observed_over_mdt25", "Observed trend as a fraction of\nthe smallest detectable",
           [0, 0.1, 0.25, 0.5, 0.75, 1.0], "magma_r")]
fig, axes = plt.subplots(1, 2, figsize=carto.size("double", 0.36), subplot_kw={"projection": PROJ})
for ax, (col, title, bounds, cmap_name) in zip(axes, panels):
    carto.map_axes(PROJ, EXTENT, ax=ax)
    cmap = plt.get_cmap(cmap_name, len(bounds) - 1)
    norm = BoundaryNorm(bounds, cmap.N)
    have = sb[sb[col].notna()]
    ax.add_geometries(sb[sb[col].isna()].geometry, ccrs.PlateCarree(), facecolor="#e0e0e0", edgecolor="white", linewidth=0.1, zorder=2)
    for v, geom in zip(have[col].clip(upper=bounds[-1] - 1e-9), have.geometry):
        ax.add_geometries([geom], ccrs.PlateCarree(), facecolor=cmap(norm(v)), edgecolor="white", linewidth=0.1, zorder=2)
    ax.add_geometries(outlines, ccrs.PlateCarree(), facecolor="none", edgecolor="#333333", linewidth=0.5, zorder=3)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.17, ticks=bounds)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=8)
carto.scale_bar(axes[1], "lower right")
carto.panel_labels(axes)
carto.finish_map(fig, R / "figures" / "map_detectability",
                 title="Detectable trend in snow-cover duration by sub-basin",
                 sources="MODIS MOD10A1 v6.1 (NASA NSIDC DAAC); HydroBASINS level 7; GMBA mountain inventory v2; analysis src/22",
                 crs="Albers equal-area, centre 85E, parallels 25N and 45N",
                 note="Grey: sub-basins with too little snow or too few years. Outlines: mountain ranges.", status="preliminary")
print("saved")
