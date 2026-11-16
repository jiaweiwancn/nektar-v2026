#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind CT15 qualification requests without changing scientific thresholds."""
import argparse,json
from pathlib import Path
from ct_accepted_feasibility import write_json
from ct_accepted_open import require
from rans_request_manifest import sha

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--ct14',type=Path,required=True);p.add_argument('--phase',choices=['timestep','pair'],required=True);p.add_argument('--runtime',type=Path,required=True);a=p.parse_args()
    template=a.ct14/'source/docs/plans/cao-tamura-accepted-14-evidence/short-pair-contract-v2.json';q=json.loads(template.read_text())
    q={k:q[k] for k in ('dt','T','physical_manifest_sha256','actual_initialized_state_sha256','slots','elements','limits','classifier','admissibility','observations','quadrature','forbidden','alpha','steps')}
    q.update(schema='CT15AdmittedShortPairContractV2',static_sha256=sha(a.root/'static-v1/results.json'),component_sha256=sha(a.root/'small-v2/results.json'),capacity_sha256=sha(a.ct14/'capacity-v4/results.json'),component_identity='CT15-ImmutableTensorV1',governing_sha256=sha(a.runtime/'IncNavierStokesSolver'),probe_sha256=sha(a.runtime/'RANSPreservationProbe'),oracle_complete_pair_sha256=sha(a.ct14/'pair-v2/results.json'),scope='fresh CT15 component only; unchanged CT11 scientific classifier; old capacity reused only for unchanged sparse algebra with new exact canonical matrix/mask equality')
    if a.phase=='timestep':
        q['previous_pair_sha256']=sha(a.ct14/'pair-v1/results.json');q.pop('timestep_admission_sha256',None)
    else:q['timestep_admission_sha256']=sha(a.root/'timestep-v1/results.json')
    dest=a.root/(a.phase+'-contract-v1.json');require(not dest.exists(),'immutable request contract already exists');write_json(dest,q)
if __name__=='__main__':main()
