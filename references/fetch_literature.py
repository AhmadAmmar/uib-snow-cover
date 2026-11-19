"""Fetch verified records for the Third Pole reading list from Crossref by search, and write refs.bib, refs.ris and a check table.

A record is kept only if its title matches the search closely; everything else is listed for manual checking.
"""
import csv
import difflib
import json
import re
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

R = Path(__file__).parent
MAILTO = "ahmadammargeo@gmail.com"
LIST = [
    ("snow", "Immerzeel Droogers de Jong Bierkens 2009 Large-scale monitoring of snow cover and runoff simulation in Himalayan river basins using remote sensing"),
    ("snow", "Gurung Maharjan Shrestha Shrestha Bajracharya Murthy 2017 Climate and topographic controls on snow cover dynamics in the Hindu Kush Himalaya"),
    ("snow", "Maskey Uhlenbrook Ojha 2011 An analysis of snow cover changes in the Himalayan region using MODIS snow products and in-situ temperature data"),
    ("snow", "Tahir Chevallier Arnaud Neppel Ahmad 2011 Modeling snowmelt-runoff under climate scenarios in the Hunza River basin, Karakoram Range"),
    ("snow", "Hasson Lucarini Khan Petitta Bolch Gioli 2014 Early 21st century snow cover state over the western river basins of the Indus River system"),
    ("snow", "Notarnicola 2020 Hotspots of snow cover changes in global mountain regions over 2000-2018"),
    ("snow", "Pepin Arnone Gobiet et al 2022 Climate changes and their elevational patterns in the mountains of the world"),
    ("snow", "Immerzeel et al 2020 Importance and vulnerability of the world's water towers"),
    ("snow", "Immerzeel van Beek Bierkens 2010 Climate change will affect the Asian water towers"),
    ("snow", "Kraaijenbrink Stigter Yao Immerzeel 2021 Climate change decisive for Asia's snow meltwater supply"),
    ("snow", "Stigter Wanders Saloranta Shea Bierkens Immerzeel 2017 Assimilation of snow cover and snow depth into a snow model to estimate snow water equivalent and snowmelt runoff in a Himalayan catchment"),
    ("snow", "Liu Fang Margulis 2021 High Mountain Asia UCLA daily snow reanalysis"),
    ("snow", "Smith Bookhagen 2018 Changes in seasonal snow water equivalent distribution in High Mountain Asia (1987 to 2009)"),
    ("modis", "Hall Riggs Salomonson DiGirolamo Bayr 2002 MODIS snow-cover products"),
    ("modis", "Riggs Hall Roman 2017 Overview of NASA's MODIS and VIIRS snow-cover Earth System Data Records"),
    ("modis", "Hall Riggs 2007 Accuracy assessment of the MODIS snow products"),
    ("modis", "Gascoin Hagolle Huc Jarlan Dejoux Szczypta Marti Sanchez 2015 A snow cover climatology for the Pyrenees from MODIS snow products"),
    ("modis", "Hall Riggs DiGirolamo Roman 2019 Evaluation of MODIS and VIIRS cloud-gap-filled snow-cover products for production of an Earth science data record"),
    ("glacier", "Hugonnet McNabb Berthier et al 2021 Accelerated global glacier mass loss in the early twenty-first century"),
    ("glacier", "Farinotti Immerzeel de Kok Quincey Dehecq 2020 Manifestations and mechanisms of the Karakoram glacier anomaly"),
    ("glacier", "Hewitt 2005 The Karakoram anomaly? Glacier expansion and the elevation effect, Karakoram Himalaya"),
    ("glacier", "Gardelle Berthier Arnaud 2012 Slight mass gain of Karakoram glaciers in the early twenty-first century"),
    ("glacier", "Bolch et al 2012 The state and fate of Himalayan glaciers"),
    ("glacier", "Azam Wagnon Berthier Vincent Fujita Kargel 2018 Review of the status and mass changes of Himalayan-Karakoram glaciers"),
    ("glacier", "Brun Berthier Wagnon Kaab Treichler 2017 A spatially resolved estimate of High Mountain Asia glacier mass balances from 2000 to 2016"),
    ("glacier", "Forsythe Fowler Li Kilsby Pritchard 2017 Karakoram temperature and glacial melt driven by regional atmospheric circulation variability"),
    ("glacier", "Bookhagen Burbank 2010 Toward a complete Himalayan hydrological budget: spatiotemporal distribution of snowmelt and rainfall and their impact on river discharge"),
    ("geomorph", "Scherler Bookhagen Strecker 2011 Spatially variable response of Himalayan glaciers to climate change affected by debris cover"),
    ("geomorph", "Bookhagen Burbank 2006 Topography, relief, and TRMM-derived rainfall variations along the Himalaya"),
    ("geomorph", "Korner Jetz Paulsen Payne Rahbek Spehn 2017 A global inventory of mountains for bio-geographical applications"),
    ("geomorph", "Snethlage Geschke Spehn et al 2022 A hierarchical inventory of the world's mountains for global comparative mountain science"),
    ("geomorph", "Lehner Grill 2013 Global river hydrography and network routing: baseline data and new approaches to study the world's large river systems"),
    ("climate", "Yao Thompson Yang et al 2012 Different glacier status with atmospheric circulations in Tibetan Plateau and surroundings"),
    ("climate", "Palazzi von Hardenberg Provenzale 2013 Precipitation in the Hindu-Kush Karakoram Himalaya: observations and future scenarios"),
    ("climate", "Dimri Niyogi Barros et al 2015 Western disturbances: a review"),
    ("climate", "Cannon Carvalho Jones Bookhagen 2015 Multi-annual variations in winter westerly disturbance activity affecting the Himalaya"),
    ("climate", "Wester Mishra Mukherji Shrestha 2019 The Hindu Kush Himalaya Assessment"),
    ("climate", "Munoz Sabater et al 2021 ERA5-Land: a state-of-the-art global reanalysis dataset for land applications"),
    ("climate", "Mountain Research Initiative EDW Working Group 2015 Elevation-dependent warming in mountain regions of the world"),
    ("climate", "Gadgil 2003 The Indian monsoon and its variability"),
    ("climate", "Saji Goswami Vinayachandran Yamagata 1999 A dipole mode in the tropical Indian Ocean"),
    ("climate", "Kumar Rajagopalan Hoerling Bates Cane 2006 Unraveling the mystery of Indian monsoon failure during El Nino"),
    ("stats", "Yue Pilon Cavadias 2002 Power of the Mann-Kendall and Spearman's rho tests for detecting monotonic trends in hydrological series"),
    ("stats", "Wilks 2016 The stippling shows statistically significant grid points: how research results are routinely overstated and overinterpreted"),
    ("stats", "Benjamini Hochberg 1995 Controlling the false discovery rate: a practical and powerful approach to multiple testing"),
    ("stats", "Hamed Rao 1998 A modified Mann-Kendall trend test for autocorrelated data"),
    ("stats", "Sen 1968 Estimates of the regression coefficient based on Kendall's tau"),
    ("carto", "Crameri Shephard Heron 2020 The misuse of colour in science communication"),
    ("carto", "Brewer Hatchard Harrower 2003 ColorBrewer in print: a catalog of color schemes for maps"),
    ("carto", "Rougier Droettboom Bourne 2014 Ten simple rules for better figures"),
    ("carto", "Savric Patterson Jenny 2018 The Equal Earth map projection"),
]


def get(url):
    """Text of a URL (40 s timeout)."""
    with urlopen(Request(url, headers={"User-Agent": f"reading-list (mailto:{MAILTO})"}), timeout=40) as r:
        return r.read().decode("utf8")


def main():
    """Fetch the reading list (one OpenAlex search per topic) and write it as a table, BibTeX and RIS."""
    rows, bib, ris = [], [], []
    for topic, q in LIST:
        title = re.sub(r"^.*?\b(?:19|20)\d\d\b\s*", "", q)
        year = re.search(r"\b(19|20)\d\d\b", q)[0]
        try:
            items = json.loads(get(f"https://api.crossref.org/works?rows=8&query.bibliographic={quote(q)}"))["message"]["items"]
        except Exception as e:
            rows.append((topic, "", year, "", 0, f"error {e}", title))
            continue
        scored = sorted(((difflib.SequenceMatcher(None, title.lower(), (i.get("title") or [""])[0].lower()).ratio(), i) for i in items), key=lambda x: -x[0])
        s, i = scored[0]
        doi, found = i.get("DOI", ""), (i.get("title") or [""])[0]
        ok = s >= 0.85
        rows.append((topic, doi, year, found, round(s, 2), "ok" if ok else "CHECK", title))
        if ok:
            try:
                b = get(f"https://api.crossref.org/works/{quote(doi)}/transform/application/x-bibtex")
                bib.append(b.strip() + "\n")
                ris.append(get(f"https://api.crossref.org/works/{quote(doi)}/transform/application/x-research-info-systems").strip() + "\n")
            except Exception as e:
                rows[-1] = rows[-1][:5] + (f"export error {e}",) + rows[-1][6:]
        time.sleep(0.3)
    (R / "refs.bib").write_text("\n".join(bib), encoding="utf8")
    (R / "refs.ris").write_text("\n".join(ris), encoding="utf8")
    with open(R / "reading_list.csv", "w", newline="", encoding="utf8") as f:
        w = csv.writer(f)
        w.writerow(["topic", "doi", "year", "title_found", "match", "status", "title_wanted"])
        w.writerows(rows)
    print(len(bib), "verified;", [r[6][:70] for r in rows if r[5] != "ok"])


if __name__ == "__main__":
    main()
