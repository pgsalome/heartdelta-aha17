"""Cohort summaries and dose-response models for AHA levels."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def group_tests(levels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    data = levels[(levels.timepoint != "rt") & levels.delta_hu_mean.notna()].copy()
    arm_col = "modality" if "modality" in data else "group"
    data[arm_col] = data[arm_col].astype(str).str.lower()
    for (tp, level), group in data.groupby(["timepoint", "level"]):
        proton = group[group[arm_col].str.contains("proton")].delta_hu_mean
        photon = group[group[arm_col].str.contains("photon")].delta_hu_mean
        test = stats.mannwhitneyu(proton, photon, alternative="two-sided") if len(proton) and len(photon) else None
        rows.append({"timepoint": tp, "level": level, "n_proton": len(proton), "n_photon": len(photon), "proton_median_delta_hu": proton.median(), "photon_median_delta_hu": photon.median(), "mannwhitney_u": test.statistic if test else np.nan, "p_value": test.pvalue if test else np.nan})
    return pd.DataFrame(rows)


def dose_response(levels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    data = levels[(levels.timepoint != "rt") & levels.delta_hu_mean.notna() & levels.mean_dose_eqd2.notna()]
    for (tp, level), group in data.groupby(["timepoint", "level"]):
        if len(group) < 3: continue
        slope, intercept, r, p, se = stats.linregress(group.mean_dose_eqd2, group.delta_hu_mean)
        rows.append({"timepoint": tp, "level": level, "n": len(group), "slope_hu_per_gy_eqd2": slope, "intercept": intercept, "r": r, "p_value": p, "slope_se": se})
    return pd.DataFrame(rows)
