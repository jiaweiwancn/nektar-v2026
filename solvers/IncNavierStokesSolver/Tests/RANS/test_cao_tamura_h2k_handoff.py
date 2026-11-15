#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Synthetic routing tests; not a physical hill certificate or native import."""
import argparse, tempfile, copy
from pathlib import Path
from unittest.mock import patch
import qualify_cao_tamura_budget_h2k as driver
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();tests={}
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp);context=root/'context.json';atomic(context,dict(synthetic=True));ledger=root/'aggregate.json';atomic(ledger,dict(synthetic=True));cert=root/'certificate.json';receipt=root/'hill-import.json'
        c=dict(state_sha256='a'*64,H2k=dict(version='H2kSingleFreshEpochV1',context_sha256=sha(context),aggregate_ledger_file=str(ledger),aggregate_ledger_sha256=sha(ledger)))
        argv=['qualification','pilot','--root',str(root),'--context',str(context),'--certificate',str(cert)]
        def invoke(document,change=None):
            atomic(cert,document);r=dict(passed=True,certificate_sha256=sha(cert),state_sha256=document.get('state_sha256'))
            if change:r.update(change)
            atomic(receipt,r)
            with patch('sys.argv',argv),patch.object(driver,'pilot') as native:
                try:driver.main()
                except ValueError:return native.call_count
                return native.call_count
        assert invoke(c)==1;tests['complete_synthetic_routing_prerequisites_delegate_once']=True
        for label,change in [('failed_import',dict(passed=False)),('changed_certificate_hash',dict(certificate_sha256='0'*64)),('changed_import_state',dict(state_sha256='b'*64))]:assert invoke(c,change)==0;tests[label+'_stops_before_native']=True
        bad=copy.deepcopy(c);bad['H2k']['version']='unknown';assert invoke(bad)==0;tests['unknown_epoch_stops_before_native']=True
        bad=copy.deepcopy(c);bad['H2k']['context_sha256']='0'*64;assert invoke(bad)==0;tests['changed_context_stops_before_native']=True
        atomic(ledger,dict(changed=True));assert invoke(c)==0;tests['changed_epoch_ledger_stops_before_native']=True
    atomic(a.output,dict(version='H2kHandoffRoutingTestsV1',passed=True,scope='synthetic precondition routing only; real producer/import evidence remains separately required',tests=tests))


if __name__=='__main__':main()
