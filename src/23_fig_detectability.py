"""Figure: smallest detectable trend in snow duration against record length, per range, with observed trends.

Input: outputs/detectability.csv (src/22). Output: figures/detectability.png and .pdf.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

R = Path(__file__).resolve().parent.parent
d = pd.read_csv(R / "outputs" / "detectability.csv")
r = d[d.level == "range"].sort_values("unit")
L = [10, 15, 20, 25, 40]
colours = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
markers = ["o", "s", "^", "D", "v", "P"]

fig, ax = plt.subplots(figsize=(7, 4.5), dpi=200)
for (_, row), c, m in zip(r.iterrows(), colours, markers):
    ax.plot(L, [row[f"mdt_{n}y_days_per_decade"] for n in L], color=c, marker=m, lw=1.5, label=row.unit)
    ax.plot(25, abs(row.sen_days_per_decade), marker=m, color=c, mfc="white", ms=8, ls="none")
ax.set_yscale("log")
ax.set_yticks([2, 3, 5, 10, 20, 30, 50])
ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
ax.minorticks_off()
ax.set_xticks(L)
ax.set_xlabel("Length of record (years)")
ax.set_ylabel("Trend in snow-cover duration (days per decade)")
ax.set_title("Smallest trend detected 80% of the time (lines) and observed trend 2001-2025 (open)", fontsize=9)
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8, frameon=False, ncol=2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(R / "figures" / f"detectability.{ext}")
print("saved figures/detectability.png")
