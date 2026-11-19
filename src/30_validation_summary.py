"""Agreement of MOD10A1 with same-day Landsat 8/9, from outputs/landsat_validation_points.csv (src/29).

Landsat is taken as reference: a MODIS cell is snow in Landsat when at least half its 30 m pixels are snow. MODIS snow
is NDSI_Snow_Cover >= 40. Overall accuracy, snow omission (Landsat snow, MODIS not) and commission (MODIS snow, Landsat
not), by range, season and 1000 m elevation band, with the number of cells and scenes.
Also the mean absolute difference between MODIS NDSI snow cover/100 and Landsat snow fraction.
Outputs: outputs/landsat_validation_summary.csv, figures/landsat_validation.png/.pdf
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent.parent


def scores(g):
    """Agreement of MODIS (NDSI 40 or more) with Landsat (snow fraction 0.5 or more): cells, scenes, accuracy, omission, commission, Landsat snow
    share and mean absolute error of the fraction."""
    ls, md = g.ls_frac >= 0.5, g.modis >= 40
    tp, tn = (ls & md).sum(), (~ls & ~md).sum()
    return pd.Series(dict(cells=len(g), scenes=g.scene.nunique(), accuracy=(tp + tn) / len(g),
                          omission=((ls & ~md).sum() / ls.sum()) if ls.sum() else np.nan,
                          commission=((md & ~ls).sum() / md.sum()) if md.sum() else np.nan,
                          landsat_snow_share=ls.mean(), mae_fraction=(g.modis.clip(0, 100) / 100 - g.ls_frac).abs().mean()))


def main():
    """Validation scores overall and by range and season."""
    d = pd.read_csv(R / "outputs" / "landsat_validation_points.csv")
    d = d[d.modis <= 100]
    d["band"] = (d.elev // 1000 * 1000).astype(int)
    out = [scores(d).to_frame().T.assign(group="all", value="all")]
    for col in ("range", "season", "band"):
        s = d.groupby(col).apply(scores, include_groups=False).reset_index().rename(columns={col: "value"}).assign(group=col)
        out.append(s)
    t = pd.concat(out)[["group", "value", "cells", "scenes", "accuracy", "omission", "commission", "landsat_snow_share", "mae_fraction"]]
    t.to_csv(R / "outputs" / "landsat_validation_summary.csv", index=False)
    print(t.round(3).to_string(index=False))

    b = t[t.group == "band"].astype({"value": int}).sort_values("value")
    b = b[b.cells >= 50]
    fig, ax = plt.subplots(figsize=(5, 3.6), dpi=200)
    mid = b.value + 500
    ax.plot(b.accuracy, mid, "-o", color="#0072B2", label="overall accuracy", ms=4)
    ax.plot(b.omission, mid, "-s", color="#D55E00", label="snow missed by MODIS (omission)", ms=4)
    ax.plot(b.commission, mid, "-^", color="#009E73", label="false snow in MODIS (commission)", ms=4)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share")
    ax.set_ylabel("Elevation (m, band midpoint)")
    ax.set_title(f"MOD10A1 against same-day Landsat 8/9: {int(t.iloc[0].cells)} cells, {int(t.iloc[0].scenes)} scenes", fontsize=8)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(R / "figures" / f"landsat_validation.{ext}")


if __name__ == "__main__":
    main()
