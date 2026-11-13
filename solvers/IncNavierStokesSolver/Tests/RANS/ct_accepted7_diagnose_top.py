#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only analytical weak localization of the retained RN preparations."""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted4_analysis import weak
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    labels=['A95-p3','A93-p3','E94E95-p3','A95-p4']
    m=RequestManifest(a.output,dict(schema='CT7ReadOnlyTopWeakLocalizationV1',source=sha(Path(__file__)),
        packets={name:sha(a.input/name/'native.json') for name in labels},constructions={name:sha(a.input/name/'construction.json') for name in labels}))
    result=dict(schema='CT7ReadOnlyTopWeakLocalizationV1',execution='complete',cases={},admission=False)
    for name in labels:
        folder=a.input/name;packet=json.loads((folder/'native.json').read_text());c=json.loads((folder/'construction.json').read_text())
        require(packet['status']=='pass','failed native preparation cannot be interpreted as a weak residual')
        q=weak(packet,c);result['cases'][name]=dict(continuous_weak=q,consumed_configuration_sha256=c['configuration_sha256'],
            packet_sha256=sha(folder/'native.json'),native_initialization=packet['initialization_receipt'],
            actual_stored_slots=len(packet['observation']['state']),
            actual_orders={k:v for k,v in packet.get('resolution',{}).items() if 'order' in k.lower()})
        print(name,{k:v['scaled_free_maximum'] for k,v in q['equations'].items()},flush=True)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0

if __name__=='__main__':raise SystemExit(main())
