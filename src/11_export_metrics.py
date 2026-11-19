"""Yearly snow metrics per pixel from daily MOD10A1, kept as Earth Engine assets so trends can be computed there.

Bands for hydrological year Y (1 Oct Y-1 to 30 Sep Y):
  snow_days, valid_days    days with NDSI snow cover at or above the threshold, and days with a clear view
  late_snow, late_valid    the same for 1 Aug to 15 Sep, the end of the melt season

    python src/11_export_metrics.py check <range>         one year at coarse scale, prints the means; exports nothing
    python src/11_export_metrics.py start <range> [first last]   one asset task per year
"""
import sys

import ee

ee.Initialize(project="ee-ahmadammargeo")
REGIONS = ee.FeatureCollection("projects/ee-ahmadammargeo/assets/tp_regions_v2_range")
ROOT = "projects/ee-ahmadammargeo/assets/tp_snow_v1"
THRESHOLD = 40
PRODUCT = "MODIS/061/MOD10A1"
BAND = "NDSI_Snow_Cover"


def counts(start, end, prefix):
    """Snow days (NDSI at or above THRESHOLD) and valid days over a period, as two bands with a prefix."""
    col = ee.ImageCollection(PRODUCT).filterDate(start, end).select(BAND)
    snow = col.map(lambda i: i.gte(THRESHOLD).unmask(0).toUint16()).sum().rename(f"{prefix}snow")
    valid = col.map(lambda i: i.mask().gt(0).toUint16().rename(BAND)).sum().rename(f"{prefix}valid")
    return snow, valid


def year_image(y):
    """Snow and valid days of a hydrological year (October-September) and of late summer (1 August-15 September)."""
    s, v = counts(f"{y - 1}-10-01", f"{y}-10-01", "t_")
    ls, lv = counts(f"{y}-08-01", f"{y}-09-16", "l_")
    return ee.Image.cat(s.rename("snow_days"), v.rename("valid_days"), ls.rename("late_snow"), lv.rename("late_valid")).toUint16().set("year", y)


def main(mode, name, first=2001, last=2025):
    """Export the yearly snow metrics of a range, first to last year."""
    geom = REGIONS.filter(ee.Filter.eq("unit", name)).filter(ee.Filter.eq("level", "range")).geometry()
    print(name, "km2", round(geom.area(100).getInfo() / 1e6))
    proj = ee.ImageCollection(PRODUCT).first().select(BAND).projection()
    if mode == "check":
        img = year_image(int(first))
        print(img.reduceRegion(ee.Reducer.mean(), geom, 5000, bestEffort=True).getInfo())
        return
    try:
        ee.data.createAsset({"type": "FOLDER"}, ROOT)
    except ee.EEException:
        pass
    key = name.lower().replace(" ", "_")
    for y in range(int(first), int(last) + 1):
        t = ee.batch.Export.image.toAsset(image=year_image(y).clip(geom), description=f"tp_snow_{key}_{y}", assetId=f"{ROOT}/{key}_{y}",
                                          crs=proj.crs().getInfo(), crsTransform=proj.getInfo()["transform"], region=geom.bounds(), maxPixels=1e11,
                                          pyramidingPolicy={".default": "mean"})
        t.start()
        print("started", y, t.id, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
