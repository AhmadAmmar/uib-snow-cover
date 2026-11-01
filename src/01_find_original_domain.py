"""Which outlet reproduces the original areas? Compares candidate basins with the original on four dates."""
import sys
from pathlib import Path

import ee
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import snowlib as S

ee.Initialize(project="ee-ahmadammargeo")
old = pd.read_csv(Path(__file__).parent.parent / "outputs" / "original_8day.csv", parse_dates=["date"]).set_index("date").area
dates = ["2001-02-26", "2006-08-13", "2010-04-07", "2015-12-27"]
dates = [d for d in dates if pd.Timestamp(d) in old.index] or [str(x.date()) for x in old.index[[40, 150, 400, 700]]]
print("original", {d: round(old[pd.Timestamp(d)]) for d in dates})
for name, pt in {"Tarbela": S.TARBELA, "Attock, with Kabul": (72.25, 33.90), "Kalabagh": (71.55, 32.96), "Panjnad, with Jhelum and Chenab": (70.98, 29.34), "Indus at Sukkur": (68.87, 27.70)}.items():
    ids = S.upstream_ids(7, *pt)
    g = ee.FeatureCollection(S.HYBAS.format(7)).filter(ee.Filter.inList("HYBAS_ID", ids)).geometry().dissolve(maxError=100)
    out = {}
    for d in dates:
        im = ee.ImageCollection(S.COLLECTION).filterDate(d, ee.Date(d).advance(1, "day")).first()
        out[d] = round(S.class_areas(im, g).get("snow").getInfo())
    print(name, "area", round(g.area(100).getInfo() / 1e6), out, flush=True)
