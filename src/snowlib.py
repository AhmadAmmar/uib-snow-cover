"""Shared Earth Engine functions for snow cover area from MOD10A2. Areas in km2."""
import ee

COLLECTION = "MODIS/061/MOD10A2"
HYBAS = "WWF/HydroSHEDS/v1/Basins/hybas_{}"
TARBELA = (72.69, 34.09)
SCALE = 463.3127
CLASSES = {"snow": [200], "cloud": [50], "land": [25], "water": [37], "ice": [100], "night": [11], "nodecision": [1], "fill": [254, 255, 0]}


def upstream_ids(level, lon, lat):
    """HydroBASINS IDs draining through the basin containing a point (the point's basin included)."""
    basins = ee.FeatureCollection(HYBAS.format(level))
    start = basins.filterBounds(ee.Geometry.Point([lon, lat])).first().get("HYBAS_ID").getInfo()
    region = basins.filterBounds(ee.Geometry.BBox(60, 25, 100, 42))
    table = region.reduceColumns(ee.Reducer.toList(2), ["HYBAS_ID", "NEXT_DOWN"]).get("list").getInfo()
    children = {}
    for hid, down in table:
        children.setdefault(int(down), []).append(int(hid))
    found, stack = [], [int(start)]
    while stack:
        h = stack.pop()
        found.append(h)
        stack += children.get(h, [])
    return found


def upper_indus(level=7):
    """Upper Indus above Tarbela, as one geometry."""
    ids = upstream_ids(level, *TARBELA)
    return ee.FeatureCollection(HYBAS.format(level)).filter(ee.Filter.inList("HYBAS_ID", ids)).geometry().dissolve(maxError=100)


def class_areas(image, geometry):
    """Area of each snow-extent class in a geometry, as a dictionary of km2."""
    band = image.select("Maximum_Snow_Extent")
    area = ee.Image.pixelArea().divide(1e6)
    bands = [area.updateMask(band.remap(codes, [1] * len(codes), 0)).rename(name) for name, codes in CLASSES.items()]
    return ee.Image.cat(bands).reduceRegion(reducer=ee.Reducer.sum(), geometry=geometry, scale=SCALE, crs=image.select(0).projection(), maxPixels=1e10, tileScale=4)


def series(geometry, region, start="2000-02-01", end="2026-12-31"):
    """One feature per 8-day composite with the area of each class."""
    def one(image):
        """Class areas of one image over the geometry, as a dated feature."""
        return ee.Feature(None, class_areas(image, geometry)).set({"region": region, "date": image.date().format("YYYY-MM-dd")})
    return ee.ImageCollection(COLLECTION).filterDate(start, end).map(one)
