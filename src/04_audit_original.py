"""Audit the original snow project, read-only: spreadsheet formulas, units, and the report's stated figures. Writes outputs/original_audit.txt."""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
O = R / "original"
out = []


def say(*a):
    """Print a line and keep it for the audit report."""
    s = " ".join(str(x) for x in a)
    out.append(s)
    print(s[:300])


wb = openpyxl.load_workbook(O / "UIBNDSI.xlsx")
for ws in wb.worksheets:
    say("SHEET", ws.title, ws.dimensions)
    for row in ws.iter_rows(min_row=1, max_row=4):
        say("  ", [str(c.value)[:28] for c in row if c.value is not None][:14])
    seen = {}
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str) and c.value.startswith("="):
                k = re.sub(r"\d+", "#", c.value)
                seen.setdefault(k, [c.coordinate, c.value, 0])[2] += 1
    for k, (co, v, n) in list(seen.items())[:12]:
        say("   F", co, "x", n, v[:110])
    say("   charts", len(ws._charts))
old = pd.read_csv(R / "outputs/original_8day.csv", parse_dates=["date"])
say("series", len(old), old.date.min().date(), old.date.max().date(), list(old.columns))
say("area min, median, max", old.area.min(), old.area.median(), old.area.max())
gap = old.date.diff().dt.days
say("steps of other than 8 days", gap[gap != 8].value_counts().head(6).to_dict(), "duplicates", int(old.date.duplicated().sum()))
old["y"] = old.date.dt.year
a = old.groupby("y").area.mean()
full = a[old.groupby("y").size() >= 44]
say("annual mean trend", stats.theilslopes(full.values, full.index.values)[0], "tau, p", stats.kendalltau(full.index, full.values))
say("linear r2 on raw 8-day", stats.linregress(old.date.map(pd.Timestamp.toordinal), old.area).rvalue ** 2)
for f in sorted(O.glob("*.docx")):
    t = subprocess.run(["pandoc", str(f), "-t", "plain", "--wrap=none"], capture_output=True).stdout.decode("utf8", "replace")
    say("\nDOC", f.name, len(t.split()), "words")
    for line in t.splitlines():
        if re.search(r"\d", line) and re.search(r"km|trend|increase|decrease|percent|%|R2|R²|correlat|cloud|threshold|NDSI|area of|mean", line, re.I):
            say("  >", line.strip()[:280])
(R / "outputs/original_audit.txt").write_text("\n".join(out), encoding="utf8")
