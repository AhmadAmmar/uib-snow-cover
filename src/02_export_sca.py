"""Export the 8-day snow area of the Upper Indus to Google Drive as CSV (folder uib_snow_cover).

    python src/02_export_sca.py check    one composite only, prints the result; exports nothing
    python src/02_export_sca.py start    starts the export task
"""
import sys
from pathlib import Path

import ee

sys.path.insert(0, str(Path(__file__).parent))
import snowlib as S

ee.Initialize(project="ee-ahmadammargeo")
geom = S.upper_indus(7)
print("Upper Indus above Tarbela, km2:", round(geom.area(100).getInfo() / 1e6))
if sys.argv[1] == "start-sukkur":
    # The Indus above Sukkur holds the snow of the original domain, whose exact boundary was lost.
    ids = S.upstream_ids(7, 68.87, 27.70)
    wide = ee.FeatureCollection(S.HYBAS.format(7)).filter(ee.Filter.inList("HYBAS_ID", ids)).geometry().dissolve(maxError=100)
    t = ee.batch.Export.table.toDrive(collection=S.series(wide, "indus_above_sukkur"), description="indus_sukkur_sca_8day", folder="uib_snow_cover", fileNamePrefix="indus_sukkur_sca_8day", fileFormat="CSV",
                                      selectors=["date", "region"] + list(S.CLASSES))
    t.start()
    print("started", t.id)
elif sys.argv[1] == "check":
    im = ee.ImageCollection(S.COLLECTION).filterDate("2000-02-26", "2000-02-27").first()
    print(im.date().format("YYYY-MM-dd").getInfo(), {k: round(v) for k, v in S.class_areas(im, geom).getInfo().items()})
else:
    t = ee.batch.Export.table.toDrive(collection=S.series(geom, "upper_indus"), description="uib_sca_8day", folder="uib_snow_cover", fileNamePrefix="uib_sca_8day", fileFormat="CSV",
                                      selectors=["date", "region"] + list(S.CLASSES))
    t.start()
    print("started", t.id)
