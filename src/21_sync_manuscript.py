"""Copy generated numbers, tables and figures into the private manuscript folder, rebuild it and its two-column preview.

Run after src/20_figures.py. The manuscript text in paper/main.tex is never written here.
"""
import shutil
from pathlib import Path

from research_helpers.paperkit import build, preview

R = Path(__file__).resolve().parent.parent
MS = Path(r"D:\OneDrive - Ulster University\PhD\Manuscripts\11_ThirdPole_Snow")
shutil.copytree(R / "docs/generated", MS / "paper/generated", dirs_exist_ok=True)
for f in (R / "figures").glob("fig_*.pdf"):
    shutil.copy2(f, MS / "figures" / f.name)
shutil.copy2(R / "references/refs.bib", MS / "paper/references.bib")
print("main", build(MS / "paper/main.tex")[:8])
print("preview", preview(MS / "paper/main.tex"))
