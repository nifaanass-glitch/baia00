#!/usr/bin/env python3
from __future__ import annotations
import json, os, time, traceback
from pathlib import Path
from khoca import InteractiveCalculator

DATA=json.loads(Path('candidates.json').read_text())
name=os.environ['CANDIDATE']
controls=DATA['controls']
pd=DATA['candidates'][name]['pd']
out={'schema':'isolated-cork-2r-khoca-v1','candidate':name,'metadata':DATA['candidates'][name],'controls':{}}

def extract_s(result):
    active=[]
    for t,q,torsion,coefficient in result[1]:
        c=int(coefficient)
        if c:
            active.extend([int(q)]*abs(c))
    active.sort()
    if len(active)!=2 or active[1]-active[0]!=2:
        return {'ok':False,'active_q':active,'s':None}
    return {'ok':True,'active_q':active,'s':(active[0]+active[1])//2}

def compute(link_pd):
    start=time.time()
    calc=InteractiveCalculator(1,(0,-1),0)
    result=calc(link_pd,verbose=True,progress=True)
    return {'elapsed_seconds':time.time()-start,'extracted':extract_s(result),'result':result}

try:
    expected={'positive_trefoil':2,'negative_trefoil':-2,'figure_eight':0}
    for label,c_pd in controls.items():
        print('CONTROL_START',label,flush=True)
        out['controls'][label]=compute(c_pd)
        print('CONTROL_DONE',label,out['controls'][label]['extracted'],flush=True)
        Path('result.json').write_text(json.dumps(out,indent=2,sort_keys=True))
    observed={k:v['extracted']['s'] for k,v in out['controls'].items()}
    out['calibration_pass']=observed==expected
    if not out['calibration_pass']:
        raise RuntimeError(f'control calibration failed: {observed}')
    print('TARGET_START',name,'crossings',len(pd),flush=True)
    out['target_Q']=compute(pd)
    out['s_Q']=out['target_Q']['extracted']['s']
    out['terminal_nonzero_s']=out['s_Q'] not in (None,0)
    print('TARGET_DONE',name,out['target_Q']['extracted'],flush=True)
except BaseException as exc:
    out['error']={'type':type(exc).__name__,'repr':repr(exc),'traceback':traceback.format_exc()}
    print('FATAL',repr(exc),flush=True)
finally:
    Path('result.json').write_text(json.dumps(out,indent=2,sort_keys=True))
