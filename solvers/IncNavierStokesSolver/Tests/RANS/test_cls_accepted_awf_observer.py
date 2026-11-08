#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused observer identity, crash-window and signed-budget checks."""
from pathlib import Path
import tempfile
from observe_cls_accepted_awf import MetricLedger,budget,window_statistics
from run_cls_epawf_picard import require


def main():
    with tempfile.TemporaryDirectory(prefix='accepted-awf-observer-') as folder:
        path=Path(folder)/'metrics.jsonl';ledger=MetricLedger(path)
        def record(step,positive,negative):
            return dict(global_step=step,source_identity=dict(kind='native_checkpoint',sha256='checkpoint-'+str(step)),
                integral_signed=positive+negative,integral_positive=positive,integral_negative=negative,
                integral_absolute=positive-negative,bound_activations=1,bound_mass_norm=.01,
                correction_raw_ratio=.5,increment_norm_identity_defect=0.,bound_seconds=1.,physical_bands=[])
        first=record(1,.1,-.05);second=record(2,.2,-.2);ledger.append(first);ledger.append(second)
        # Recover after observer metric append but before cursor publication.
        for i in range(2):
            ledger=MetricLedger(path);ledger.append(second)
            require(len(ledger.records)==2 and len(path.read_text().splitlines())==2,'idempotent observer recovery duplicated a step')
        totals=budget(ledger.records,0,2,.000625)
        require(totals['complete_coverage'] and abs(totals['signed']-.05)<1e-15 and
            abs(totals['positive']-.3)<1e-15 and abs(totals['negative']+.25)<1e-15 and
            abs(totals['absolute']-.55)<1e-15,'signed integral concealed positive/negative redistribution')
        require(abs(totals['mean_correction_rate']['signed']-40)<1e-12,'correction amount/rate units mixed')
        incomplete=budget(ledger.records,0,3,.000625);require(not incomplete['complete_coverage'],'missing prefix budget was declared complete')
        with path.open('a') as f:f.write('{"global_step":3')
        repaired=MetricLedger(path);require(len(repaired.records)==2 and len(path.read_text().splitlines())==2,'partial final observer record was not repaired')
        conflict=dict(second,source_identity=dict(kind='native_checkpoint',sha256='different'))
        try:repaired.append(conflict)
        except ValueError as error:require('conflicting' in str(error),'unexpected conflict rejection')
        else:raise AssertionError('accepted checkpoint identity conflict was hidden')
    observations={n:dict(diagnostic=dict(mean_wall_stress=float(n),fixed_bands=dict(regions=[dict(
        name='fixed',x_region='all',k_integral=2.+n,production=3.+n,destruction=1.)]))) for n in (0,1,2)}
    window=window_statistics(observations,0,2,.25)
    require(window['sample_count']==3 and window['maximum_sample_gap']==.25 and
        window['wall_stress']['mean']==1. and abs(window['wall_stress']['rms_fluctuation']-2**-.5)<1e-15,
        'physical window moments or time weights differ')
    band=window['fixed_bands']['fixed:all']
    require(band['k_integral']['mean']==3. and band['production_time_integral_estimate']==2. and
        band['destruction_time_integral_estimate']==.5,'sampled source rates confused with integrated amounts')
    require(window_statistics(observations,0,3,.25) is None,'missing physical endpoint hidden')
    print('accepted-AWF observer idempotence, partial write, identity and signed/rate accounting passed')


if __name__=='__main__':main()
