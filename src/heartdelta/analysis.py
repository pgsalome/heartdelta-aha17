"""Cohort summaries and dose-response models for AHA levels."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


def canonical_arm(value, fallback=""):
    text=f"{value} {fallback}".lower()
    if "proton" in text or "pencil" in text or "pbs" in text: return "proton"
    if any(x in text for x in ("photon","vmat","imrt","3d","conformal")): return "photon"
    return "unknown"


def complete_level_cohort(levels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply manuscript wall-band QC and retain identical levels per follow-up."""
    data=levels[(levels.timepoint!="rt") & levels.segment.isna()].copy()
    current=(data.hu_median.between(0,100)&(data.fraction_hu_0_100>=.8))
    baseline=(data.baseline_hu_median.between(0,100)&(data.baseline_fraction_hu_0_100>=.8))
    data["level_qc_pass"]=current&baseline
    audit=[]; keep=[]
    for (case,tp),group in data.groupby(["case_id","timepoint"]):
        levels_present=set(group.loc[group.level_qc_pass,"level"]); passed=levels_present=={"basal","mid","apical"}
        audit.append({"case_id":case,"timepoint":tp,"levels_present":len(levels_present),"complete_three_level_qc":passed})
        if passed: keep.append(group[group.level_qc_pass])
    return (pd.concat(keep,ignore_index=True) if keep else data.iloc[0:0]),pd.DataFrame(audit)


def mixed_effects_models(levels: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Fit the manuscript random-intercept attenuation models separately by FU."""
    data=levels.rename(columns={"delta_hu_mean":"delta_hu","mean_dose_eqd2":"dose"}).copy()
    data["level"]=pd.Categorical(data.level,categories=["basal","mid","apical"])
    data["arm"]=[canonical_arm(a,b) for a,b in zip(data.get("modality",pd.Series("",index=data.index)),data.get("group",pd.Series("",index=data.index)))]
    formulas={"common_slope":"delta_hu ~ dose + C(level)","dose_by_level":"delta_hu ~ dose * C(level)","dose_by_modality":"delta_hu ~ dose * C(arm) + C(level)"}
    terms=[]; diagnostics=[]
    for tp,part in data.groupby("timepoint"):
        for name,formula in formulas.items():
            try:
                result=smf.mixedlm(formula,part,groups=part.case_id,re_formula="1").fit(reml=True,method="powell",maxiter=4000,disp=False)
                diagnostics.append({"timepoint":tp,"model":name,"formula":formula+" + (1 | case_id)","n_patients":part.case_id.nunique(),"n_observations":len(part),"converged":bool(result.converged),"random_intercept_variance":float(result.cov_re.iloc[0,0]),"residual_variance":float(result.scale)})
                ci=result.conf_int()
                for term,value in result.fe_params.items(): terms.append({"timepoint":tp,"model":name,"term":term,"estimate":value,"ci95_low":ci.loc[term,0],"ci95_high":ci.loc[term,1],"p_value":result.pvalues[term],"n_patients":part.case_id.nunique(),"n_observations":len(part)})
                if name=="dose_by_level":
                    names=list(result.fe_params.index); covariance=result.cov_params().loc[names,names].to_numpy()
                    for label,weights in (("basal dose slope",{"dose":1}),("mid dose slope",{"dose":1,"dose:C(level)[T.mid]":1}),("apical dose slope",{"dose":1,"dose:C(level)[T.apical]":1})):
                        vector=np.array([weights.get(x,0) for x in names]); estimate=float(vector@result.fe_params.to_numpy()); se=float(np.sqrt(vector@covariance@vector)); terms.append({"timepoint":tp,"model":name,"term":label,"estimate":estimate,"ci95_low":estimate-1.96*se,"ci95_high":estimate+1.96*se,"p_value":float(2*stats.norm.sf(abs(estimate/se))),"n_patients":part.case_id.nunique(),"n_observations":len(part)})
            except (ValueError,np.linalg.LinAlgError) as exc:
                diagnostics.append({"timepoint":tp,"model":name,"formula":formula+" + (1 | case_id)","n_patients":part.case_id.nunique(),"n_observations":len(part),"converged":False,"error":str(exc)})
        if tp=="fu1":
            photon=part[part.arm=="photon"]
            if photon.case_id.nunique()>=3:
                try:
                    result=smf.mixedlm(formulas["common_slope"],photon,groups=photon.case_id,re_formula="1").fit(reml=True,method="powell",maxiter=4000,disp=False); ci=result.conf_int()
                    for term,value in result.fe_params.items(): terms.append({"timepoint":tp,"model":"photon_only_common_slope","term":term,"estimate":value,"ci95_low":ci.loc[term,0],"ci95_high":ci.loc[term,1],"p_value":result.pvalues[term],"n_patients":photon.case_id.nunique(),"n_observations":len(photon)})
                except (ValueError,np.linalg.LinAlgError): pass
    return pd.DataFrame(terms),pd.DataFrame(diagnostics)


def group_tests(levels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    data = levels[(levels.timepoint != "rt") & levels.delta_hu_mean.notna()].copy()
    data["_arm"]=[canonical_arm(a,b) for a,b in zip(data.get("modality",pd.Series("",index=data.index)),data.get("group",pd.Series("",index=data.index)))]
    for (tp, level), group in data.groupby(["timepoint", "level"]):
        proton = group[group._arm.eq("proton")].delta_hu_mean
        photon = group[group._arm.eq("photon")].delta_hu_mean
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
