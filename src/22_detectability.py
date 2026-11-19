"""Smallest trend in snow duration that each region's record can detect.

For every range, division and sub-basin: remove the Sen slope from yearly duration, measure the residual
standard deviation and lag-1 autocorrelation, then simulate AR(1) noise plus a linear trend and find the trend
a two-sided Mann-Kendall test (alpha 0.05) detects in 80% of runs, for records of 10, 15, 20, 25 and 40 years.

Inputs: outputs/unit_year.csv, outputs/subbasin_year.csv.
Outputs: outputs/detectability.csv, outputs/detectability_checks.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent.parent
OUT = R / "outputs"
LENGTHS = [10, 15, 20, 25, 40]
SIMS = 4000
POWER = 0.8
ALPHA = 0.05
rng = np.random.default_rng(20261114)


def mk_z(y):
    """Mann-Kendall Z for each row of y (no tie correction; series are continuous)."""
    n = y.shape[1]
    i, j = np.triu_indices(n, 1)
    s = np.sign(y[:, j] - y[:, i]).sum(axis=1)
    var = n * (n - 1) * (2 * n + 5) / 18
    return (s - np.sign(s)) / np.sqrt(var)


def sen(x, y):
    """Theil-Sen slope: median of all pairwise slopes."""
    i, j = np.triu_indices(len(y), 1)
    return np.median((y[j] - y[i]) / (x[j] - x[i]))


def ar1_noise(n, sd, phi, sims):
    """sims series of AR(1) noise of length n, standard deviation sd and lag-1 correlation phi (clipped to 0-0.9)."""
    phi = float(np.clip(phi, 0, 0.9))
    e = rng.standard_normal((sims, n)) * sd * np.sqrt(1 - phi**2)
    out = np.empty_like(e)
    out[:, 0] = rng.standard_normal(sims) * sd
    for t in range(1, n):
        out[:, t] = phi * out[:, t - 1] + e[:, t]
    return out


def power(n, sd, phi, slope, noise=None):
    """Power: share of AR(1) series with the given slope in which the Mann-Kendall test is significant
    at ALPHA (two-sided)."""
    noise = ar1_noise(n, sd, phi, SIMS) if noise is None else noise
    y = noise + slope * np.arange(n)
    return np.mean(np.abs(mk_z(y)) > stats.norm.ppf(1 - ALPHA / 2))


def mdt(n, sd, phi):
    """Smallest slope (days per year) with power >= 0.8, by bisection on common noise."""
    noise = ar1_noise(n, sd, phi, SIMS)
    lo, hi = 0.0, 4 * sd / np.sqrt(n)
    while power(n, sd, phi, hi, noise) < POWER:
        hi *= 2
    for _ in range(25):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if power(n, sd, phi, mid, noise) < POWER else (lo, mid)
    return hi


PHIS = np.round(np.arange(0, 0.91, 0.05), 2)


def mdt_table():
    """MDT for unit noise (sd 1) on a grid of record length and lag-1 autocorrelation. MDT scales with sd."""
    return {n: np.array([mdt(n, 1.0, p) for p in PHIS]) for n in LENGTHS}


def mdt_lookup(table, n, sd, phi):
    """Minimum detectable trend for series length n, noise sd and phi, interpolated from the table."""
    return sd * np.interp(np.clip(phi, 0, 0.9), PHIS, table[n])


def checks():
    """Known-answer checks: false positive rate near 5% for n=25, Mann-Kendall p against scipy's Kendall tau,
    Sen slope of 2x+1, and the white-noise MDT against the regression value 0.078."""
    rows = []
    z = mk_z(rng.standard_normal((20000, 25)))
    rows.append(("False positive rate, no trend, n=25", np.mean(np.abs(z) > 1.96), "close to 0.05"))
    y = rng.standard_normal(25)
    rows.append(("MK Z matches scipy Kendall tau p", 2 * stats.norm.sf(abs(mk_z(y[None])[0])),
                 f"scipy tau p {stats.kendalltau(np.arange(25), y).pvalue:.4f}"))
    rows.append(("Sen slope of 2x+1", sen(np.arange(10.0), 2 * np.arange(10.0) + 1), "2"))
    # Known case: a t-test on a regression slope needs (1.96+0.84)*sd/sqrt(sum((t-tbar)^2)); MK needs slightly more
    rows.append(("MDT, white noise sd 1, n 25", mdt(25, 1.0, 0.0), "slightly above the regression-slope value 2.80/sqrt(1300) = 0.078"))
    pd.DataFrame(rows, columns=["check", "value", "expected"]).to_csv(OUT / "detectability_checks.csv", index=False)
    print(pd.DataFrame(rows, columns=["check", "value", "expected"]).to_string(index=False))


def main():
    """Minimum detectable trends for each unit's length, noise and autocorrelation."""
    checks()
    table = mdt_table()
    pd.DataFrame({f"n{n}": v for n, v in table.items()}, index=PHIS).rename_axis("lag1").to_csv(
        OUT / "detectability_unit_noise.csv")
    u =pd.read_csv(OUT / "unit_year.csv")
    s = pd.read_csv(OUT / "subbasin_year.csv").assign(level="sub-basin")
    rows = []
    for (unit, level, parent, km2), g in pd.concat([u, s]).groupby(["unit", "level", "parent", "km2"]):
        g = g.sort_values("year")
        x, y = g.year.to_numpy(float), g.duration.to_numpy(float)
        if len(y) < 20 or np.nanmean(y) < 5:
            continue
        b = sen(x, y)
        res = y - b * x
        res -= res.mean()
        sd = res.std(ddof=1)
        phi = np.corrcoef(res[:-1], res[1:])[0, 1]
        z = mk_z(y[None])[0]
        row = dict(unit=unit, level=level, parent=parent, km2=km2, years=len(y), mean_days=y.mean(),
                   sen_days_per_decade=10 * b, mk_p=2 * stats.norm.sf(abs(z)), resid_sd=sd, lag1=phi)
        for n in LENGTHS:
            row[f"mdt_{n}y_days_per_decade"] = 10 * mdt_lookup(table, n, sd, phi)
        row["mdt25_pct_of_mean_per_decade"] = 100 * row["mdt_25y_days_per_decade"] / y.mean()
        row["observed_over_mdt25"] = abs(row["sen_days_per_decade"]) / row["mdt_25y_days_per_decade"]
        rows.append(row)
    d = pd.DataFrame(rows).sort_values(["level", "parent", "unit"])
    d.to_csv(OUT / "detectability.csv", index=False)
    print(d[["unit", "level", "mean_days", "sen_days_per_decade", "mk_p", "resid_sd", "lag1",
             "mdt_25y_days_per_decade", "observed_over_mdt25"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
