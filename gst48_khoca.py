#!/usr/bin/env python3
from __future__ import annotations
import json, time, traceback
from pathlib import Path
from khoca import InteractiveCalculator

CONTROLS = {
    'positive_trefoil': [[1,5,2,4],[5,3,0,2],[3,1,4,0]],
    'negative_trefoil': [[1,4,2,5],[5,2,0,3],[3,0,4,1]],
    'figure_eight': [[3,7,4,6],[1,4,2,5],[7,3,0,2],[5,0,6,1]],
}
EXPECTED = {'positive_trefoil': 2, 'negative_trefoil': -2, 'figure_eight': 0}
GST_PD = [[8,35,9,36],[9,24,10,25],[10,81,11,82],[11,55,12,54],[53,13,54,12],[82,13,83,14],[25,14,26,15],[36,15,37,16],[16,27,17,28],[17,38,18,39],[69,19,70,18],[19,84,20,85],[20,52,21,51],[21,67,22,66],[55,23,56,22],[80,23,81,24],[37,26,38,27],[28,39,29,40],[70,30,71,29],[30,85,31,86],[31,51,32,50],[32,66,33,65],[56,34,57,33],[79,34,80,35],[71,41,72,40],[41,86,42,87],[42,50,43,49],[43,65,44,64],[57,45,58,44],[78,45,79,46],[46,77,47,78],[47,59,48,58],[48,63,49,64],[67,52,68,53],[2,60,3,59],[91,61,92,60],[61,93,62,92],[62,4,63,3],[68,84,69,83],[87,73,88,72],[73,4,74,5],[74,93,75,94],[90,75,91,76],[1,76,2,77],[5,89,6,88],[94,90,95,89],[95,7,0,6],[7,1,8,0]]

out = {
    'schema': 'gst48-calibrated-rational-lee-v1',
    'provenance': 'Regina ExampleLink::gst(), converted exactly from Link::fromData signs/component data',
    'crossings': len(GST_PD),
    'controls': {},
}

def extract_s(result):
    active = []
    for t, q, torsion, coefficient in result[1]:
        c = int(coefficient)
        if c:
            active.extend([int(q)] * abs(c))
    active.sort()
    if len(active) != 2 or active[1] - active[0] != 2:
        return {'ok': False, 'active_q': active, 's': None, 'reason': 'unexpected Lee E_infinity shape'}
    return {'ok': True, 'active_q': active, 's': (active[0] + active[1]) // 2}

def compute(pd):
    start = time.time()
    calc = InteractiveCalculator(1, (0, -1), 0)
    result = calc(pd, verbose=True, progress=True)
    return {'elapsed_seconds': time.time() - start, 'extracted': extract_s(result), 'result': result}

try:
    for name, pd in CONTROLS.items():
        print('CONTROL_START', name, flush=True)
        out['controls'][name] = compute(pd)
        print('CONTROL_DONE', name, out['controls'][name]['extracted'], flush=True)
        Path('gst48_result.json').write_text(json.dumps(out, indent=2, sort_keys=True))
    observed = {k: v['extracted']['s'] for k, v in out['controls'].items()}
    out['calibration_pass'] = observed == EXPECTED
    if not out['calibration_pass']:
        raise RuntimeError(f'control calibration failed: {observed}')
    print('GST_START', len(GST_PD), flush=True)
    out['GST_Q'] = compute(GST_PD)
    out['s_Q'] = out['GST_Q']['extracted']['s']
    out['terminal_nonzero_s'] = out['s_Q'] not in (None, 0)
    print('GST_DONE', out['GST_Q']['extracted'], flush=True)
except BaseException as exc:
    out['error'] = {'type': type(exc).__name__, 'repr': repr(exc), 'traceback': traceback.format_exc()}
    print('FATAL', repr(exc), flush=True)
finally:
    Path('gst48_result.json').write_text(json.dumps(out, indent=2, sort_keys=True))
