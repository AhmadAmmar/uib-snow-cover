"""Summarise the surviving original 8-day snow areas: coverage, seasonal cycle and trend. Writes outputs/original_8day.csv."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
d = pd.read_excel(R / "original" / "UIBNDSI.xlsx")
d["year"] = d.Day.str[:4].astype(int)
d["doy"] = d.Day.str[4:7].astype(int)
d["date"] = pd.to_datetime(d.year.astype(str) + d.doy.astype(str).str.zfill(3), format="%Y%j")
d = d.rename(columns={"SUM_AREA_G": "area"})[["date", "year", "doy", "area", "FREQUENCY"]].sort_values("date")
(R / "outputs").mkdir(exist_ok=True)
d.to_csv(R / "outputs" / "original_8day.csv", index=False)
print("rows", len(d), d.date.min().date(), d.date.max().date(), "duplicates", d.date.duplicated().sum(), "frequency", d.FREQUENCY.unique()[:5])
print("area min/median/max", d.area.min(), d.area.median(), d.area.max())
print("composites per year", d.groupby("year").size().to_dict())
print("step between composites, days", d.date.diff().dt.days.value_counts().head(5).to_dict())
month = d.groupby(d.date.dt.month).area.agg(["mean", "min", "max"]).round(0)
print(month)
full = d[d.groupby("year").year.transform("size") >= 44]
for name, s in [("annual mean", full.groupby("year").area.mean()), ("annual max", full.groupby("year").area.max()), ("annual min", full.groupby("year").area.min())]:
    t = stats.theilslopes(s.values, s.index.values)
    k = stats.kendalltau(s.index.values, s.values)
    print(f"{name}: years {s.index.min()}-{s.index.max()} n={len(s)} Sen slope {t.slope:.0f} [{t.low_slope:.0f}, {t.high_slope:.0f}] per year; Kendall tau {k.statistic:.2f} p={k.pvalue:.3f}")
d["jump"] = d.area.diff().abs()
print("largest 8-day jumps", d.nlargest(5, "jump")[["date", "area", "jump"]].round(0).to_string(index=False))
