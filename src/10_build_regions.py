"""Build Third Pole analysis regions, version 2, and upload them to Earth Engine.

Ranges and divisions come from GMBA Mountain Inventory v2 (levels 3 and 4). Sub-basins are HydroBASINS level 5 clipped
to the mountain domain, each named by the major basin that holds most of it.

    python src/10_build_regions.py <GMBA standard .shp> build      writes data/regions/tp_regions_v2.geojson and a preview
    python src/10_build_regions.py <GMBA standard .shp> upload     starts the export to an Earth Engine asset
"""
import json
import sys
from pathlib import Path

import ee
import geopandas as gpd
import matplotlib
import pandas as pd
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union
from shapely.validation import make_valid

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parent.parent
OUT = R / "data/regions"
RANGES = ["Hindu Kush", "Karakoram", "Himalaya", "Pamir Mountains", "Tian Shan", "Tibetan Plateau"]
ASSET = "projects/ee-ahmadammargeo/assets/tp_regions_v2"
SIMPLIFY_DEG = 0.002
MIN_SUBBASIN_KM2 = 500
EQUAL_AREA = "ESRI:54034"
V1 = Path.home() / "dev/gee-backup/ahmadammargeo/assets/data/tp_regions_v1.geojson"


def polygons(geom):
    """Valid polygonal part only; repairs can leave stray lines and points."""
    geom = make_valid(geom)
    parts = [x for x in getattr(geom, "geoms", [geom]) if x.geom_type in ("Polygon", "MultiPolygon")]
    return unary_union(parts) if parts else Polygon()


def clean(g):
    g = g.copy()
    g["geometry"] = [polygons(polygons(x).simplify(SIMPLIFY_DEG, preserve_topology=True)) for x in g.geometry]
    return g


def build(shp):
    g = gpd.read_file(shp, columns=["GMBA_V2_ID", "MapName", "Hier_Lvl", "Area", "Level_03", "Level_04"])
    g = g[g.Level_03.isin(RANGES)]
    rows = []
    for _, r in g[g.Hier_Lvl.astype(str) == "3"].iterrows():
        rows.append(dict(unit=r.MapName, level="range", parent="Third Pole", gmba_id=r.GMBA_V2_ID, source_km2=r.Area, geometry=r.geometry))
    for _, r in g[g.Hier_Lvl.astype(str) == "4"].iterrows():
        rows.append(dict(unit=r.MapName, level="division", parent=r.Level_03, gmba_id=r.GMBA_V2_ID, source_km2=r.Area, geometry=r.geometry))
    mountains = gpd.GeoDataFrame(rows, crs=g.crs)
    domain = polygons(unary_union([polygons(x) for x in mountains[mountains.level == "range"].geometry]))

    ee.Initialize(project="ee-ahmadammargeo")
    basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_5").filterBounds(ee.Geometry.BBox(*domain.bounds))
    basins = basins.map(lambda f: f.simplify(2000).select(["HYBAS_ID", "MAIN_BAS", "PFAF_ID"]))
    info = basins.getInfo()["features"]
    b = gpd.GeoDataFrame([{**f["properties"], "geometry": shape(f["geometry"])} for f in info], crs=4326)
    b["geometry"] = [polygons(polygons(x).intersection(domain)) for x in b.geometry]
    b = b[~b.geometry.is_empty]
    b = b[b.to_crs(EQUAL_AREA).area / 1e6 >= MIN_SUBBASIN_KM2]

    # the major basin name comes from the largest overlap with the named basins of version 1
    v1 = gpd.read_file(V1)
    v1 = v1[v1.level == "basin"][["unit", "geometry"]].copy()
    v1["geometry"] = [polygons(x) for x in v1.geometry]
    names = {}
    for _, r in b.iterrows():
        hit = v1.intersection(r.geometry)
        area = gpd.GeoSeries(hit, crs=4326).to_crs(EQUAL_AREA).area
        names[r.HYBAS_ID] = v1.unit.iloc[int(area.values.argmax())] if area.max() > 0 else "Other"
    subs = gpd.GeoDataFrame([dict(unit=f"{names[r.HYBAS_ID]} {int(r.HYBAS_ID)}", level="subbasin", parent=names[r.HYBAS_ID], gmba_id=int(r.HYBAS_ID), source_km2=None, geometry=r.geometry) for _, r in b.iterrows()], crs=4326)

    all_units = pd.concat([mountains, subs], ignore_index=True)
    all_units = clean(all_units)
    all_units["km2"] = (all_units.to_crs(EQUAL_AREA).area / 1e6).round(0)

    # checks: outlines must stay valid and keep their area after simplifying
    assert all_units.is_valid.all()
    all_units = all_units[all_units.km2 >= MIN_SUBBASIN_KM2].reset_index(drop=True)  # too small to report
    known = all_units[all_units.source_km2.notna()]
    drift = ((known.km2 - known.source_km2).abs() / known.source_km2)
    print("largest area drift against GMBA:", round(float(drift.max()), 4), known.loc[drift.idxmax(), "unit"])
    print(known.assign(d=drift).nlargest(5,"d")[["unit","km2","source_km2","d"]].round(3).to_string())
    assert drift.max() < 0.02
    d = all_units[all_units.level == "division"].groupby("parent").km2.sum()
    r_ = all_units[all_units.level == "range"].set_index("unit").km2
    print("divisions / range area:", (d / r_).round(3).to_dict())

    OUT.mkdir(parents=True, exist_ok=True)
    all_units.drop(columns="source_km2").to_file(OUT / "tp_regions_v2.geojson", driver="GeoJSON")
    print(all_units.level.value_counts().to_dict(), "domain km2", int(all_units[all_units.level == "range"].km2.sum()))
    print(all_units[all_units.level == "subbasin"].groupby("parent").km2.agg(["count", "sum"]).astype(int).to_string())
    print("file MB", round((OUT / "tp_regions_v2.geojson").stat().st_size / 1e6, 1))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), dpi=130)
    for a, lev, col in zip(ax, ("range", "division", "subbasin"), ("unit", "unit", "parent")):
        all_units[all_units.level == lev].plot(ax=a, column=col, cmap="tab20", edgecolor="white", linewidth=.2)
        a.set_title(lev)
        a.set_axis_off()
    fig.tight_layout()
    fig.savefig(OUT / "tp_regions_v2_preview.png")


def upload():
    """One asset per level, because a single request is limited to 10 MB."""
    ee.Initialize(project="ee-ahmadammargeo")
    data = json.loads((OUT / "tp_regions_v2.geojson").read_text(encoding="utf8"))
    for level in ("range", "division", "subbasin"):
        feats = [f for f in data["features"] if f["properties"]["level"] == level]
        print(level, len(feats), "MB", round(len(json.dumps(feats)) / 1e6, 1))
        fc = ee.FeatureCollection([ee.Feature(ee.Geometry(f["geometry"], "EPSG:4326", False), {k: v for k, v in f["properties"].items() if v is not None}) for f in feats])
        t = ee.batch.Export.table.toAsset(collection=fc, description=f"tp_regions_v2_{level}", assetId=f"{ASSET}_{level}")
        t.start()
        print("started", t.id, flush=True)


if __name__ == "__main__":
    build(sys.argv[1]) if sys.argv[2] == "build" else upload()
