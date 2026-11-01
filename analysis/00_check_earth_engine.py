"""Which snow, climate and boundary datasets are available in Earth Engine, with their date ranges. Read-only."""
import ee

ee.Initialize(project="ee-ahmadammargeo")
for cid in ["MODIS/061/MOD10A1", "MODIS/061/MYD10A1", "MODIS/061/MOD10A2", "MODIS/061/MYD10A2", "MODIS/061/MOD10A1F", "MODIS/006/MOD10A1F", "NASA/SNOW/MOD10A1F", "ECMWF/ERA5_LAND/MONTHLY_AGGR", "ECMWF/ERA5/MONTHLY", "UCSB-CHG/CHIRPS/PENTAD", "NASA/GPM_L3/IMERG_MONTHLY_V07", "NASA/VIIRS/002/VNP10A1F", "NOAA/VIIRS/001/VNP10A1F", "NOAA/CDR/SST_PATHFINDER/V53", "NOAA/CDR/OISST/V2_1"]:
    try:
        c = ee.ImageCollection(cid)
        print(cid, c.size().getInfo(), ee.Date(c.aggregate_min("system:time_start")).format("YYYY-MM-dd").getInfo(), ee.Date(c.aggregate_max("system:time_start")).format("YYYY-MM-dd").getInfo(), c.first().bandNames().getInfo()[:6])
    except Exception as e:
        print(cid, "NOT AVAILABLE", str(e)[:70])
for lvl in (3, 4):
    b = ee.FeatureCollection(f"WWF/HydroSHEDS/v1/Basins/hybas_{lvl}").filterBounds(ee.Geometry.Point([74.5, 35.5]))
    print("basin level", lvl, b.first().get("HYBAS_ID").getInfo(), round(ee.Number(b.first().get("SUB_AREA")).getInfo()))
