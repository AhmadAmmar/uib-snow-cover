"""Evidence map for the Third Pole work: search OpenAlex, save records, and count coverage by region, variable, method and period.

Counts rest on keyword matches in titles and abstracts, so they locate candidate gaps; each must be confirmed by reading.
Writes corpus.csv, corpus.ris, coverage_*.csv.
"""
import csv
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd

sys.stdout.reconfigure(encoding="utf8", errors="replace")
R = Path(__file__).parent
MAIL = "ahmadammargeo@gmail.com"
Q = [
    "snow cover trend High Mountain Asia MODIS", "snow cover Hindu Kush Himalaya trend", "snow cover Upper Indus Basin MODIS", "snow cover Karakoram", "snow cover Tibetan Plateau trend",
    "snow cover Tian Shan trend", "snow cover Pamir", "snow cover duration mountain trend satellite", "snow phenology onset melt High Mountain Asia", "snowline altitude Karakoram Himalaya",
    "cloud gap filled MODIS snow High Mountain Asia", "long term AVHRR snow cover Tibetan Plateau", "snow cover climate drivers Himalaya temperature precipitation", "snow water equivalent High Mountain Asia reanalysis",
    "western disturbances snowfall Karakoram Himalaya", "Indian monsoon snow cover Himalaya relationship", "Eurasian snow cover Indian summer monsoon", "ENSO snow cover Tibetan Plateau Himalaya", "Indian Ocean Dipole snow High Mountain Asia",
    "El Nino winter precipitation western Himalaya", "NAO western disturbance Karakoram snow", "Karakoram anomaly glacier mass balance", "High Mountain Asia glacier mass change", "glacier snow runoff Indus Basin contribution",
    "elevation dependent warming Himalaya snow", "snow cover sub-basin Indus Ganges Brahmaputra", "false discovery field significance snow cover trend", "snow cover dataset validation Himalaya Landsat Sentinel", "ICIMOD snow cover Hindu Kush Himalaya",
    "debris covered glacier snow cover Himalaya", "snow drought High Mountain Asia", "black carbon dust snow Himalaya albedo",
]
REGION = {"Hindu Kush": r"hindu.?kush(?!.{0,3}himalaya)", "Karakoram": "karakoram", "Himalaya": "himalaya", "Pamir": "pamir", "Tian Shan": r"tian.?shan|tien.?shan", "Tibetan Plateau": r"tibet|qinghai", "Hengduan": "hengduan", "Upper Indus": r"indus", "Ganges/Brahmaputra": r"ganges|brahmaputra", "Whole HMA": r"high mountain asia|third pole|hindu.?kush.?himalaya|\bhkh\b"}
VAR = {"Snow cover": r"snow.?cover|snow extent|\bsca\b", "Duration/phenology": r"snow.?(cover )?(duration|phenology|onset|melt.?out)|\bscd\b", "Snowline": "snow.?line", "SWE": r"snow water equivalent|\bswe\b", "Glacier mass": r"glacier.{0,30}mass|mass balance", "Runoff": r"runoff|streamflow|discharge"}
DRIVER = {"Temperature": "temperature|warming", "Precipitation": r"precipitation|snowfall", "Monsoon": "monsoon", "ENSO": r"enso|el ni.o", "IOD": r"indian ocean dipole|\biod\b", "Western disturbance": r"westerl|western disturbance", "NAO/AO": r"north atlantic oscillation|\bnao\b|arctic oscillation", "Aerosol": r"black carbon|aerosol|dust"}
RIGOR = {"Field significance or FDR": r"false discovery|field significance|benjamini", "Autocorrelation handled": r"pre.?whiten|autocorrelat", "Cloud sensitivity": r"cloud.{0,40}(sensitiv|gap.?fill|remov|filter)|gap.?fill", "Validation with stations or Landsat": r"validat|ground.?truth|landsat|station", "Open data or code": r"github|zenodo|open.?access data|publicly available|code is available"}
PERIOD = {"Before 1990 start": r"19[5-8]\d[-– ]+(to )?(19|20)\d\d", "Ends after 2020": r"(?:19|20)\d\d[-– ]+(to )?20(2[1-6])"}


def get(url):
    """JSON text of a URL, up to five attempts."""
    for i in range(5):
        try:
            with urlopen(Request(url, headers={"User-Agent": f"gap-map (mailto:{MAIL})"}), timeout=50) as r:
                return json.loads(r.read().decode("utf8"))
        except Exception:
            time.sleep(2 ** i)
    return {"results": []}


def abstract(inv):
    """Abstract text rebuilt from an OpenAlex inverted index."""
    if not inv:
        return ""
    w = sorted((p, k) for k, v in inv.items() for p in v)
    return " ".join(k for _, k in w)


def search(q):
    """OpenAlex works for a query, by relevance and by citations."""
    out = []
    for sort in ("relevance_score:desc", "cited_by_count:desc"):
        d = get(f"https://api.openalex.org/works?search={quote(q)}&per-page=100&sort={sort}&filter=from_publication_date:2000-01-01,type:article|review&mailto={MAIL}&select=id,doi,title,publication_year,cited_by_count,abstract_inverted_index,primary_location,authorships")
        out += [(q, w) for w in d.get("results", [])]
    return out


def main():
    """Run the gap searches in parallel and write the unique hits with abstracts."""
    with ThreadPoolExecutor(5) as ex:
        hits = [x for r in ex.map(search, Q) for x in r]
    keep = re.compile(r"snow|glacier|cryospher|monsoon|western disturbance", re.I)
    area = re.compile(r"himalaya|karakoram|hindu.?kush|tibet|pamir|tian.?shan|tien.?shan|indus|ganges|brahmaputra|high mountain asia|third pole|hengduan|qinghai", re.I)
    rec = {}
    for q, w in hits:
        t, a = w.get("title") or "", abstract(w.get("abstract_inverted_index"))
        if not (keep.search(t + a) and area.search(t + a)):
            continue
        src = ((w.get("primary_location") or {}).get("source") or {}).get("display_name", "")
        au = (w.get("authorships") or [{}])
        first = ((au[0].get("author") or {}).get("display_name", "") if au else "")
        rec[w["id"]] = dict(first_author=first, year=w.get("publication_year"), title=t, journal=src, doi=(w.get("doi") or "").replace("https://doi.org/", ""), cited_by=w.get("cited_by_count"), has_abstract=bool(a), abstract=a)
    d = pd.DataFrame(rec.values()).sort_values(["cited_by"], ascending=False)
    txt = (d.title + " " + d.abstract).str.lower()
    has = lambda p: txt.str.contains(p, regex=True, case=False)
    regions = {k: has(v) for k, v in REGION.items()}
    v_ = {k: has(v) for k, v in VAR.items()}
    dr = {k: has(v) for k, v in DRIVER.items()}
    ri = {k: has(v) for k, v in RIGOR.items()}
    print("records", len(d), "with abstract", int(d.has_abstract.sum()))
    a = pd.DataFrame({r: {k: int((m & x).sum()) for k, x in v_.items()} for r, m in regions.items()}); a.to_csv(R / "coverage_region_variable.csv"); print(a.to_string())
    sc = has(VAR["Snow cover"]) | has(VAR["Duration/phenology"])
    b = pd.DataFrame({r: {k: int((m & x & sc).sum()) for k, x in dr.items()} for r, m in regions.items()}); b.to_csv(R / "coverage_region_snow_driver.csv"); print("\nSNOW papers by driver\n", b.to_string())
    c = pd.DataFrame({r: {k: int((m & x & sc).sum()) for k, x in ri.items()} for r, m in regions.items()}); c.loc["All snow papers"] = [int((m & sc).sum()) for m in regions.values()]; c.to_csv(R / "coverage_region_snow_rigour.csv"); print("\nSNOW papers by rigour\n", c.to_string())
    for name, mask in (("snow + ENSO + region", sc & has(DRIVER["ENSO"])), ("snow + FDR", sc & has(RIGOR["Field significance or FDR"])), ("snow duration + subbasin", has(VAR["Duration/phenology"]) & has(r"sub.?basin")), ("snow + division/sub-range comparison across ranges", sc & has(r"(karakoram).*(himalaya|hindu)") & has(VAR["Duration/phenology"]))):
        print("\n", name, int(mask.sum()))
        print(d[mask].head(8)[["first_author", "year", "cited_by", "title"]].to_string(index=False, max_colwidth=85))
    d.drop(columns="abstract").to_csv(R / "corpus.csv", index=False, encoding="utf-8-sig")
    d.to_csv(R / "corpus_with_abstracts.csv", index=False, encoding="utf-8-sig")
    with open(R / "corpus.ris", "w", encoding="utf8") as f:
        for r in d.itertuples():
            f.write(f"TY  - JOUR\nTI  - {r.title}\nAU  - {r.first_author}\nPY  - {r.year}\nJO  - {r.journal}\nDO  - {r.doi}\nAB  - {r.abstract[:3000]}\nER  - \n\n")


if __name__ == "__main__":
    main()
