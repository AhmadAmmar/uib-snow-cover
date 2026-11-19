"""Association of yearly snow duration with large-scale climate indices, by season and with a one-year lag.

Indices (downloaded to data/indices, 2026-11-14):
  ONI   ENSO, NOAA CPC  https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt (3-month seasons; month = centre)
  DMI   Indian Ocean Dipole, HadISST, NOAA PSL  https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data
  NAO   North Atlantic Oscillation, NOAA CPC  norm.nao.monthly.b5001.current.ascii; a proxy often linked to
        western disturbances, not a measure of them.
Seasons of hydrological year Y: autumn = Sep-Nov of Y-1, winter = Dec of Y-1 to Mar of Y. Lag 1 uses year Y-1's season.
Spearman correlation for each unit, index, season and lag; Benjamini-Hochberg false discovery control over all tests.
Associations, not causes. Output: outputs/indices_attribution.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
D = R / "data" / "indices"


def oni():
    """Monthly ONI (ENSO) values, seasons centred on their middle month."""
    t = pd.read_csv(D / "oni.txt", sep=r"\s+")
    centre = {"DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12}
    return pd.DataFrame({"year": t.YR, "month": t.SEAS.map(centre), "value": t.ANOM})


def dmi():
    """Monthly Dipole Mode Index (IOD)."""
    rows = []
    for line in (D / "dmi.txt").read_text().splitlines()[1:]:
        p = line.split()
        if len(p) == 13 and p[0].isdigit():
            rows += [(int(p[0]), m + 1, float(v)) for m, v in enumerate(p[1:]) if float(v) > -999]
    return pd.DataFrame(rows, columns=["year", "month", "value"])


def nao():
    """Monthly NAO index."""
    return pd.read_csv(D / "nao.txt", sep=r"\s+", header=None, names=["year", "month", "value"])


def season(ix, y, months):
    """Mean over (year offset, month) pairs relative to hydrological year y."""
    v = [ix.loc[(ix.year == y + dy) & (ix.month == m), "value"] for dy, m in months]
    v = [x.iloc[0] for x in v if len(x)]
    return np.mean(v) if len(v) == len(months) else np.nan


SEASONS = {"autumn": [(-1, 9), (-1, 10), (-1, 11)], "winter": [(-1, 12), (0, 1), (0, 2), (0, 3)]}


def main():
    """Correlation of each unit's snow duration with ENSO, IOD and NAO."""
    u = pd.read_csv(R / "outputs" / "unit_year.csv")
    idx = {"ENSO (ONI)": oni(), "IOD (DMI)": dmi(), "NAO": nao()}
    rows = []
    for name, ix in idx.items():
        for s, months in SEASONS.items():
            for lag in (0, 1):
                v = {y: season(ix, y - lag, months) for y in range(2000, 2026)}
                for (unit, level), g in u.groupby(["unit", "level"]):
                    x = g.year.map(v)
                    ok = x.notna()
                    r, p = stats.spearmanr(x[ok], g.duration[ok])
                    rows.append(dict(unit=unit, level=level, index=name, season=s, lag=lag, rho=r, p=p, n=int(ok.sum())))
    t = pd.DataFrame(rows)
    t["q"] = stats.false_discovery_control(t.p)
    t.to_csv(R / "outputs" / "indices_attribution.csv", index=False)
    print(len(t), "tests;", (t.p < 0.05).sum(), "with p < 0.05;", (t.q < 0.05).sum(), "with q < 0.05")
    print(t.sort_values("p").head(12).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
