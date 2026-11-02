"""Gain and loss area in the Tibetan Plateau at 500 m, by batch, using the threshold from src/16. Writes outputs/tibet_change_500m.csv."""
import sys, time
from pathlib import Path
import ee, pandas as pd
ee.Initialize(project="ee-ahmadammargeo")
R = Path(__file__).resolve().parent.parent
A = "projects/ee-ahmadammargeo/assets/"; T = A + "tp_table_tibet_area"
thr = float(pd.read_csv(R / "outputs/tibet_fdr_500m.csv").threshold.iloc[0])
rng = ee.FeatureCollection(A + "tp_regions_v2_range").filter(ee.Filter.eq("unit", "Tibetan Plateau")).geometry()
m = ee.Image(A + "tp_snow_v1_mk/tibetan_plateau"); s, p, n = m.select("s"), m.select("p"), m.select("n")
sig = p.lte(thr).And(n.gte(20)); area = ee.Image.pixelArea().divide(1e6)
g = lambda mask: area.updateMask(mask).reduceRegion(ee.Reducer.sum(), rng, 463.3127, maxPixels=1e11, tileScale=16).get("area")
ee.batch.Export.table.toAsset(ee.FeatureCollection([ee.Feature(ee.Geometry.Point([85, 30]), {"gain": g(sig.And(s.gt(0))), "loss": g(sig.And(s.lt(0)))})]), "tp_table_tibet_area", T).start()
while True:
    x = [t for t in ee.data.getTaskList()[:6] if t["description"] == "tp_table_tibet_area"][0]
    if x["state"] in ("COMPLETED", "FAILED"): print(x["state"], x.get("error_message", "")); break
    time.sleep(60)
r = ee.Feature(ee.FeatureCollection(T).first()).toDictionary().getInfo(); print(r)
pd.DataFrame([r]).to_csv(R / "outputs/tibet_change_500m.csv", index=False)
