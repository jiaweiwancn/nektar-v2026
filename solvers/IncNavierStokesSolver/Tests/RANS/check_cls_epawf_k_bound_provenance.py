#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check effective bound defaults and intentional checkpoint incompatibility."""
import argparse
from pathlib import Path
import platform

from run_cls_epawf_pimple import prepare, identity
from run_cls_epawf_picard import load, require, run, sha, write_json
from compare_cls_epawf_pimple import same_history


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', type=Path, required=True)
    p.add_argument('--checkpoint', type=Path, required=True)
    p.add_argument('--checkpoint-report', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--old-off-report', type=Path)
    p.add_argument('--current-off-report', type=Path)
    a = p.parse_args(); require(platform.system() == 'Linux', 'authoritative Linux check')
    a.output = a.output.resolve(); a.output.mkdir(parents=True, exist_ok=False)
    a.solver = a.solver.resolve(strict=True); a.checkpoint = a.checkpoint.resolve()
    reference = load(a.checkpoint_report)
    require(reference['k_bound_mode'] == 'MassProjection', 'bound checkpoint required')
    result = dict(status='running', solver=identity(a.solver), cases={})
    for name, mode, floor in (('mode-mismatch', 'Off', None),
                              ('floor-mismatch', 'MassProjection', 1e-6),
                              ('explicit-default', 'MassProjection', 1e-8)):
        xml = prepare(a.output/name, reference['dt'], 0, True, restart=a.checkpoint,
                      k_bound=mode, floor_relative=floor)
        receipt = run([str(a.solver), 'case.xml'], xml.parent, 'solver.log', 900)
        r = load(xml.with_suffix('.evolution.json'))
        if name == 'explicit-default':
            require(receipt['exit_code'] == 0 and same_history(r, reference), 'explicit effective default changed restart')
        else:
            require(receipt['exit_code'] != 0 and 'signature' in r['failure'], 'changed method accepted by checkpoint')
        result['cases'][name] = dict(receipt=receipt, failure=r['failure'], mode=mode, floor=floor,
                                    xml_sha256=sha(xml), report_sha256=sha(xml.with_suffix('.evolution.json')))
    if a.old_off_report or a.current_off_report:
        require(a.old_off_report and a.current_off_report, 'both Off reports required')
        require(same_history(load(a.old_off_report), load(a.current_off_report)), 'old Off checkpoint compatibility changed')
        result['old_off_bitwise_state_time_all_native_history'] = dict(
            old_report=str(a.old_off_report.resolve()), old_sha256=sha(a.old_off_report),
            current_report=str(a.current_off_report.resolve()), current_sha256=sha(a.current_off_report))
    result['status'] = 'passed'; write_json(a.output/'summary.json', result)
    print('k-bound checkpoint provenance checks passed')


if __name__ == '__main__':
    main()
