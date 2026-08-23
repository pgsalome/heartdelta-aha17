from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from .aha import generate
from .metrics import aggregate_levels,extract
from .registry import load,validate
from .workflow import run as run_workflow

def run_validate(a):
    errors=validate(a.registry,not a.no_check_files)
    if errors: print('\n'.join('ERROR: '+x for x in errors)); return 1
    print(f'Registry is valid: {len(load(a.registry))} cases'); return 0
def run_generate(a):
    cases=load(a.registry); root=Path(a.output); root.mkdir(parents=True,exist_ok=True); manifest={}
    for c in cases:
        target=root/c['case_id']/ 'aha17.nii.gz'; target.parent.mkdir(parents=True,exist_ok=True)
        generate(c['lv_mask'],c['rv_mask'],target,a.shell_mm); manifest[c['case_id']]=str(target)
    (root/'aha_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); print(f'Generated {len(manifest)} AHA maps'); return 0
def run_extract(a):
    frames=[]
    for c in load(a.registry): frames.append(extract(c,a.alpha_beta))
    segments=pd.concat(frames,ignore_index=True); root=Path(a.output); root.mkdir(parents=True,exist_ok=True)
    segments.to_csv(root/'segment_metrics.csv',index=False); aggregate_levels(segments).to_csv(root/'level_metrics.csv',index=False)
    print(f'Extracted {len(segments)} segment rows from {segments.case_id.nunique()} cases'); return 0
def run_all(a):
    print(json.dumps(run_workflow(a.config),indent=2)); return 0
def parser():
    p=argparse.ArgumentParser(prog='heartdelta'); s=p.add_subparsers(dest='command',required=True)
    x=s.add_parser('registry-validate'); x.add_argument('registry'); x.add_argument('--no-check-files',action='store_true'); x.set_defaults(func=run_validate)
    x=s.add_parser('aha-generate'); x.add_argument('--registry',required=True); x.add_argument('--output',required=True); x.add_argument('--shell-mm',type=float,default=8); x.set_defaults(func=run_generate)
    x=s.add_parser('extract'); x.add_argument('--registry',required=True); x.add_argument('--output',required=True); x.add_argument('--alpha-beta',type=float,default=2); x.set_defaults(func=run_extract)
    x=s.add_parser('run'); x.add_argument('config'); x.set_defaults(func=run_all)
    return p
def main(argv=None):
    a=parser().parse_args(argv)
    try: return int(a.func(a))
    except (ValueError,FileNotFoundError,RuntimeError) as e: print(f'ERROR: {e}'); return 1
if __name__=='__main__': raise SystemExit(main())
