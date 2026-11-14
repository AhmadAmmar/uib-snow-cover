"""Snow-covered area from the 8-day MOD10A2 product, as most published studies measure it, for the units they cover.

Per 8-day composite and unit (Maximum_Snow_Extent: 200 snow, 25 land without snow, 50 cloud):
  sca_clear   snow / (snow + land): cloud left out
  sca_all     snow / (snow + land + cloud): cloud counted as no snow, as when cloud is not filtered
  cloud       cloud / (snow + land + cloud)
Yearly values are means over composites starting in the hydrological year (1 Oct to 30 Sep), as in src/11.
Output: outputs/eightday_year.csv (unit, level, year, sca_clear, sca_all, cloud, composites). Resumes.
"""
import time
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
OUT = R / "outputs" / "eightday_year.csv"
REG = "projects/ee-ahmadammargeo/assets/tp_regions_v2_"
UNITS = [("Tibetan Plateau", "range"), ("Karakoram", "range"), ("Kashmir Himalayas", "division"),
         ("Central Himalayas", "division"), ("Kunlun Mountains", "division"), ("Himalaya", "range"),
         ("Hindu Kush", "range"), ("Pamir Mountains", "range"), ("Tian Shan", "range")]
COL = ee.ImageCollection("MODIS/061/MOD10A2").select("Maximum_Snow_Extent")


def per_composite(geom, tile_scale=4):
    def f(i):
        snow, land, cloud = i.eq(200), i.eq(25), i.eq(50)
        s = ee.Image.cat(snow.rename("s"), land.rename("l"), cloud.rename("c")).reduceRegion(
            ee.Reducer.sum(), geom, 2000, maxPixels=1e12, tileScale=tile_scale)
        return ee.Feature(None, s).set("date", i.date().format("YYYY-MM-dd"))
    return f


def main():
    done = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["unit", "year"])
    have = set(zip(done["unit"], done["year"]))
    for unit, level in UNITS:
        geom = ee.FeatureCollection(REG + level).filter(ee.Filter.eq("unit", unit)).geometry()
        for y in range(2001, 2026):
            if (unit, y) in have:
                continue
            for attempt, ts in enumerate((4, 8, 16, 16)):          # on a timeout, split the work into smaller tiles
                try:
                    fc = COL.filterDate(f"{y - 1}-10-01", f"{y}-10-01").map(per_composite(geom, ts)).getInfo()
                    break
                except ee.EEException as e:
                    if attempt == 3 or not any(s in str(e) for s in ("timed out", "Too Many Requests")):
                        raise
                    print("retry", unit, y, "tileScale", ts * 2 if ts < 16 else ts, str(e)[:40], flush=True)
                    time.sleep(60)
            d = pd.DataFrame([f["properties"] for f in fc["features"]])
            tot = d.s + d.l + d.c
            row = dict(unit=unit, level=level, year=y, sca_clear=(d.s / (d.s + d.l)).mean(), sca_all=(d.s / tot).mean(),
                       cloud=(d.c / tot).mean(), composites=len(d))
            pd.DataFrame([row]).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
            print(unit, y, round(row["sca_clear"], 3), flush=True)


if __name__ == "__main__":
    main()
