"""Reusable AHA-17 cohort figures."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LEVEL_ORDER = ["basal", "mid", "apical"]


def grouped_boxplots(levels: pd.DataFrame, output: str | Path) -> Path:
    data = levels[(levels.timepoint != "rt") & levels.delta_hu_mean.notna()].copy()
    arm = "modality" if "modality" in data else "group"
    data[arm] = data[arm].astype(str).str.lower()
    tps = [x for x in ("fu1", "fu2", "fu3") if x in set(data.timepoint)]
    fig, axes = plt.subplots(len(LEVEL_ORDER), len(tps), figsize=(6 * len(tps), 4.4 * 3), sharey=True, squeeze=False)
    for i, level in enumerate(LEVEL_ORDER):
        for j, tp in enumerate(tps):
            ax = axes[i, j]; part = data[(data.level == level) & (data.timepoint == tp)]
            vals = [part[part[arm].str.contains(x)].delta_hu_mean.dropna().to_numpy() for x in ("proton", "photon")]
            ax.boxplot(vals, labels=[f"Proton\n(n={len(vals[0])})", f"Photon\n(n={len(vals[1])})"], patch_artist=True, boxprops={"facecolor":"#d9eaf4"})
            for k, values in enumerate(vals, 1): ax.scatter(np.full(len(values), k) + np.linspace(-.07, .07, len(values)), values, s=18)
            ax.axhline(0, color="#3182bd", ls="--", lw=.8); ax.set_title(f"{tp.upper()}: {level.title()} LV")
            if j == 0: ax.set_ylabel("Mean attenuation change, ΔHU")
    fig.tight_layout(); target = Path(output); target.parent.mkdir(parents=True, exist_ok=True); fig.savefig(target, dpi=300); plt.close(fig); return target


def dose_response_panels(levels: pd.DataFrame, output: str | Path) -> Path:
    data = levels[(levels.timepoint != "rt") & levels.delta_hu_mean.notna() & levels.mean_dose_eqd2.notna()].copy()
    arm = "modality" if "modality" in data else "group"; data[arm] = data[arm].astype(str).str.lower()
    tps = [x for x in ("fu1", "fu2", "fu3") if x in set(data.timepoint)]
    fig, axes = plt.subplots(len(tps), 3, figsize=(15, 4.2 * len(tps)), sharex=True, sharey=True, squeeze=False)
    for i, tp in enumerate(tps):
        for j, level in enumerate(LEVEL_ORDER):
            ax = axes[i,j]; part = data[(data.timepoint == tp) & (data.level == level)]
            for name, marker, color in (("proton","o","#3182bd"),("photon","^","#e6550d")):
                p = part[part[arm].str.contains(name)]; ax.scatter(p.mean_dose_eqd2, p.delta_hu_mean, marker=marker, color=color, label=name.title())
            if len(part) >= 3:
                fit = np.polyfit(part.mean_dose_eqd2, part.delta_hu_mean, 1); x = np.linspace(part.mean_dose_eqd2.min(), part.mean_dose_eqd2.max(), 100); ax.plot(x, np.polyval(fit,x), color="black")
            ax.axhline(0,color="grey",ls="--",lw=.7); ax.set_title(f"{tp.upper()}: {level.title()}"); ax.set_xlabel("Mean dose (Gy EQD2)")
            if j == 0: ax.set_ylabel("Mean attenuation change, ΔHU")
    axes[-1,1].legend(loc="upper center", bbox_to_anchor=(.5,-.22), ncol=2); fig.tight_layout(); target=Path(output); target.parent.mkdir(parents=True,exist_ok=True); fig.savefig(target,dpi=300,bbox_inches="tight"); plt.close(fig); return target


def bullseye(values: dict[int, float], output: str | Path, title: str = "AHA-17 profile") -> Path:
    fig, ax = plt.subplots(figsize=(7,7), subplot_kw={"projection":"polar"}); cmap=plt.cm.coolwarm; finite=np.array([v for v in values.values() if np.isfinite(v)]); lim=max(abs(finite.min()),abs(finite.max())) if len(finite) else 1
    for start, count, ids in ((2/3,6,range(1,7)),(1/3,6,range(7,13)),(.1,4,range(13,17))):
        width=2*np.pi/count
        for k, seg in enumerate(ids):
            val=values.get(seg,np.nan); ax.bar(k*width,width=width,bottom=start,height=(1/3 if start>0.2 else .233),align="edge",color=cmap((val+lim)/(2*lim)) if np.isfinite(val) else "lightgrey",edgecolor="white"); ax.text((k+.5)*width,start+(.165 if start>0.2 else .116),str(seg),ha="center",va="center",weight="bold")
    val=values.get(17,np.nan); ax.bar(0,width=2*np.pi,bottom=0,height=.1,color=cmap((val+lim)/(2*lim)) if np.isfinite(val) else "lightgrey",edgecolor="white"); ax.text(0,.05,"17",ha="center",va="center",weight="bold"); ax.set_axis_off(); ax.set_title(title); target=Path(output); target.parent.mkdir(parents=True,exist_ok=True); fig.savefig(target,dpi=300,bbox_inches="tight"); plt.close(fig); return target
