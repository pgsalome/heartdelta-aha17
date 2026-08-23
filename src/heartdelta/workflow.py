"""Restartable raw-to-results heart workflow."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
from .analysis import dose_response, group_tests
from .figures import dose_response_panels, grouped_boxplots
from .native import aggregate_levels, extract_case
from .registry import load
from .segmentation import qc_labels, segment_scan
from .slab import extract_case as extract_slab_case


def run(config_path: str | Path) -> dict[str, int]:
    config_path = Path(config_path); cfg = json.loads(config_path.read_text())
    resolve=lambda p: str((config_path.parent / p).resolve()) if p and not Path(p).is_absolute() else p
    registry = resolve(cfg["registry"]); processed = Path(resolve(cfg.get("processed_dir", "processed"))); output = Path(resolve(cfg.get("output_dir", "outputs")))
    processed.mkdir(parents=True,exist_ok=True); output.mkdir(parents=True,exist_ok=True)
    cases=load(registry); do_segment=bool(cfg.get("segment", True)); overwrite=bool(cfg.get("overwrite",False)); rows=[]; qc=[]; status=[]; slab_rows=[]; slab_qc=[]
    for case in cases:
        for tp in ("rt","fu1","fu2","fu3"):
            ct = case.get(f"ct_{tp}") if tp != "rt" else case.get("ct_rt") or case.get("baseline_ct")
            aha_key=f"aha17_{tp}"; aha=case.get(aha_key) if tp != "rt" else case.get("aha17_rt") or case.get("aha17")
            if ct and (not aha or not Path(aha).is_file()) and do_segment:
                aha=segment_scan(ct, processed / str(case["case_id"]) / tp, overwrite)["aha17"]; case[aha_key]=aha
            if aha and Path(aha).is_file():
                case[aha_key]=aha; qc.append({"case_id":case["case_id"],"timepoint":tp,"aha17":aha,**qc_labels(aha)})
            elif aha: case[aha_key]=None
        try:
            frame=extract_case(case,float(cfg.get("alpha_beta",3.0)))
            if not frame.empty: rows.append(frame)
            slab,slab_status=extract_slab_case(case,float(cfg.get("slab_mm",8.0)),tuple(cfg.get("myocardial_band_mm",[2.0,5.0])))
            if not slab.empty: slab_rows.append(slab)
            if not slab_status.empty: slab_qc.append(slab_status)
            status.append({"case_id":case["case_id"],"status":"complete" if not frame.empty else "no_evaluable_timepoints","rows":len(frame),"error":""})
        except (OSError,RuntimeError,ValueError,KeyError) as exc:
            status.append({"case_id":case["case_id"],"status":"failed","rows":0,"error":str(exc)})
    segments=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); segments.to_csv(output/"segment_metrics.csv",index=False)
    slab=pd.concat(slab_rows,ignore_index=True) if slab_rows else pd.DataFrame(); slab.to_csv(output/"native_slab_metrics.csv",index=False)
    (pd.concat(slab_qc,ignore_index=True) if slab_qc else pd.DataFrame()).to_csv(output/"native_slab_qc.csv",index=False)
    levels=aggregate_levels(segments) if not segments.empty else pd.DataFrame(); levels.to_csv(output/"level_metrics.csv",index=False); pd.DataFrame(qc).to_csv(output/"aha_qc.csv",index=False); pd.DataFrame(status).to_csv(output/"processing_status.csv",index=False)
    analysis_levels=levels
    if not slab.empty:
        analysis_levels=slab[slab.segment.isna()].copy()
        if not levels.empty:
            dose=levels[["case_id","timepoint","level","mean_dose_eqd2"]].drop_duplicates(["case_id","timepoint","level"])
            analysis_levels=analysis_levels.merge(dose,on=["case_id","timepoint","level"],how="left")
        analysis_levels.to_csv(output/"native_slab_level_metrics.csv",index=False)
    if not analysis_levels.empty:
        group_tests(analysis_levels).to_csv(output/"group_tests.csv",index=False); dose_response(analysis_levels).to_csv(output/"dose_response.csv",index=False)
        grouped_boxplots(analysis_levels,output/"level_group_boxplots.png"); dose_response_panels(analysis_levels,output/"level_dose_response.png")
    manifest={"cases":len(cases),"cases_complete":sum(x["status"]=="complete" for x in status),"cases_failed":sum(x["status"]=="failed" for x in status),"segment_rows":len(segments),"level_rows":len(levels),"native_slab_rows":len(slab),"qc_rows":len(qc)}; (output/"run_summary.json").write_text(json.dumps(manifest,indent=2)+"\n"); return manifest


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("config"); args=parser.parse_args(argv); print(json.dumps(run(args.config),indent=2)); return 0
