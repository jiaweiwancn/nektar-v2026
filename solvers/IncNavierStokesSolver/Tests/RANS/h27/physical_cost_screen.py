#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare warm source19 startup work and a bounded same-mesh stock control."""

import argparse
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path('/home/jay/build/wj-h27-coupling-r1-20261117')
NAME = 'cp3-short-screen-v19'
OUT = ROOT/'science'/NAME
REFERENCE = ROOT/'science/cutoff-panels-v19/physical.xml'
TABLE = ROOT/'science/reference-tables-v1.json'
CASE_CAPS = {
    'paired': (.005, 12, 3, 8, 3e-12, 130),
    'dt02': (.02, 12, 4, 8, 3e-12, 180),
    'dt04': (.04, 6, 2, 8, 3e-12, 120),
    'dt08': (.08, 3, 1, 8, 3e-12, 90),
    'dt04-age24': (.04, 6, 2, 24, 3e-12, 120),
    'dt04-loose': (.04, 6, 2, 8, 1e-10, 120),
}
STOCK = tuple(f'stock-{i}' for i in range(1, 4))
COMPARISONS = (('dt04', 'dt02'), ('dt08', 'dt04'), ('dt04-age24', 'dt04'), ('dt04-loose', 'dt04'))
STEP = re.compile(r'^\s*Steps:\s+(\d+)\s+Time:\s+(\S+)\s+CPU Time:\s+([^\s]+)s', re.M)
WALL = re.compile(r'^H27 child elapsed seconds: ([.0-9]+)$', re.M)


def fixtures():
    from physical_oracle import Reference
    from physical_reference_refinement import read, require
    from resources import atomic, identity

    OUT.mkdir()
    ref = Reference(read(TABLE))
    r = f'(y-{ref.d0:.17g})'
    # Reconstruct the same continuous historical velocity trace in the stock
    # expression parser. The irrelevant outer log below H stays finite.
    side = f'sign({r}-{ref.H:.17g})'
    inner = f'{ref.tau/ref.a:.17g}*log(1+{ref.a/ref.D0:.17g}*{r})'
    outer = f'{ref.utau/ref.kappa:.17g}*log(max({r},{ref.H:.17g})/{ref.z0:.17g})'
    startup = f'0.5*((1-{side})*({inner})+(1+{side})*({outer}))'
    reference_geometry = ET.tostring(ET.parse(REFERENCE).find('GEOMETRY'))
    cases = {}
    for case in (*CASE_CAPS, *STOCK):
        root = ET.parse(REFERENCE).getroot()
        require(ET.tostring(root.find('GEOMETRY')) == reference_geometry and
                len(root.findall('GEOMETRY/ELEMENT/Q')) == 700, 'retained geometry differs')
        c = root.find('CONDITIONS')
        params = {a.text.split('=')[0].strip(): a for a in c.findall('PARAMETERS/P')}
        opts = {a.get('PROPERTY'): a for a in c.findall('SOLVERINFO/I')}
        path = OUT/(case+'.xml')
        if case in STOCK:
            for key, item in list(opts.items()):
                if key.startswith('WJ'):
                    c.find('SOLVERINFO').remove(item)
            for key, item in list(params.items()):
                if key.startswith('WJ'):
                    c.find('PARAMETERS').remove(item)
            ET.SubElement(c.find('PARAMETERS'), 'P').text = 'CheckAbortSteps = 1'
            opts['SolverType'].set('VALUE', 'VCSWeakPressure')
            for element in root.findall('EXPANSIONS/E'):
                element.set('NUMPOINTS', '8,12')
                element.set('POINTSTYPE', 'GaussGaussLegendre,GaussGaussLegendre')
                if element.get('FIELDS') != 'p':
                    element.set('FIELDS', 'u,v')
            variables = c.find('VARIABLES')
            variables.clear()
            for index, var in enumerate(('u', 'v', 'p')):
                ET.SubElement(variables, 'V', ID=str(index)).text = var
            for region in c.findall('BOUNDARYCONDITIONS/REGION'):
                for element in list(region):
                    if element.get('VAR') in ('k', 'epsilon'):
                        region.remove(element)
                side = int(region.get('REF'))
                pressure = region.find("N[@VAR='p']")
                if side == 1:
                    pressure.tag = 'D'
                else:
                    pressure.set('USERDEFINEDTYPE', 'H')
                if side == 0:
                    region.find("D[@VAR='u']").set('VALUE', startup)
            initial = c.find("FUNCTION[@NAME='InitialConditions']")
            for element in list(initial):
                if element.get('VAR') in ('k', 'epsilon'):
                    initial.remove(element)
                if element.get('VAR') == 'u':
                    element.set('VALUE', startup)
        else:
            dt, steps, cadence, age, tolerance, _ = CASE_CAPS[case]
            options = dict(TimeStep=dt, NumSteps=steps, IO_CheckSteps=cadence, WJScientificSteps=cadence,
                           WJStepPCMaximumUses=age, WJCoupledTolerance=tolerance)
            for key, value in options.items():
                params[key].text = f'{key} = {value}'
            opts['WJOutputPrefix'].set('VALUE', str(OUT/case))
        ET.indent(root)
        path.write_text(ET.tostring(root, encoding='unicode')+'\n')
        option_path = path.with_suffix('.xml.opt')
        shutil.copy2(REFERENCE.with_suffix('.xml.opt'), option_path)
        cases[case] = dict(session=identity(path), options=identity(option_path))
    atomic(OUT/'inputs.json', dict(schema='WJH27Source19ShortCostInputsV1', reference=identity(REFERENCE),
                                  session_generator=identity(__file__), runtime=identity(ROOT/'science/runtime-v19.json'),
                                  source_recipe=CASE_CAPS, stock_velocity_expression=startup, cases=cases))
    print('SOURCE19_COST_SCREEN_FIXTURES_READY', flush=True)


def assess(output):
    from physical_reference_refinement import SOURCE, read, require
    from resources import atomic, identity

    known = read(OUT/'inputs.json')['cases']
    rows = {}
    for case in known:
        operation = ROOT/'operations'/(NAME+'-'+case)
        receipt = read(operation/'receipt.json')
        result = dict(receipt=identity(operation/'receipt.json'), status=receipt['status'],
                      whole_seconds=receipt['seconds'], peak_tree_rss_bytes=receipt.get('peak_tree_rss_bytes'))
        rows[case] = result
        if receipt['status'] != 'PASS':
            continue
        if case in STOCK:
            log = (operation/'command.log').read_text()
            steps = [(int(n), float(t), float(s)) for n, t, s in STEP.findall(log)]
            elapsed = WALL.findall(log)
            fields = [OUT/(case+'.fld'), *sorted(OUT.glob(case+'_*.chk'))]
            result.update(stable_complete=len(steps) == 12 and [s[0] for s in steps] == list(range(1, 13)) and
                          abs(steps[-1][1]-.06) <= 1e-12 and all(math.isfinite(s[2]) and s[2] > 0 for s in steps) and
                          len(fields) >= 5 and all(p.is_file() and p.stat().st_size for p in fields) and len(elapsed) == 1,
                          stock_child_wall_seconds=float(elapsed[0]) if len(elapsed) == 1 else None,
                          step_internal_seconds=steps, files=[identity(p) for p in fields if p.is_file()])
            if result['stable_complete']:
                result['warm_internal_seconds_per_step'] = [statistics.mean(s[2] for s in steps[i:i+3])
                                                            for i in (3, 6, 9)]
            continue
        path = OUT/(case+'.result.json')
        run = read(path)
        dt, steps, cadence, _, _, _ = CASE_CAPS[case]
        state = run['state']
        require(run['normal_main'] and state['source_identity'] == SOURCE and state['step'] == steps and
                state['time'] == steps*dt, f'unexpected coupled output {path}')
        internal = [json.loads(line) for line in (OUT/(case+'.steps.jsonl')).read_text().splitlines()]
        complete = run['complete_step_seconds']
        require(len(internal) == len(complete) == steps, 'missing complete accepted-step timing')
        disjoint = {k: sum(s['disjoint_seconds'][k] for s in internal) for k in internal[-1]['disjoint_seconds']}
        require(abs(sum(disjoint.values())-sum(s['seconds'] for s in internal)) < 1e-7, 'disjoint row sum differs')
        outside_engine = sum(complete)-sum(disjoint.values())
        scientific = sum(run['scientific_frame_seconds']) if isinstance(run['scientific_frame_seconds'], list) else run['scientific_frame_seconds']
        recovery = sum(run['recovery_seconds'])
        require(outside_engine >= scientific+recovery-1e-7, 'native recurring timing overlaps engine categories')
        process_overhead = receipt['seconds']-sum(complete)-sum(run[key] for key in
                          ('construction_seconds', 'initialisation_seconds', 'final_export_seconds'))
        require(process_overhead >= 0, 'complete process accounting has a negative remainder')
        newton = [s['coupled_newton'] for s in internal]
        require(all(row['counts_complete'] for row in newton), 'Newton callback counts are incomplete')
        other_advancing = outside_engine-scientific-recovery
        unattributed = disjoint['other']+other_advancing+process_overhead
        result.update(native=identity(path), state_sha256=state['state_sha256'], dt=dt, time=state['time'],
                      steps=steps, cadence=cadence, complete_step_seconds=complete, disjoint_engine=disjoint,
                      initialization=run['initialisation_seconds'], construction=run['construction_seconds'],
                      science=scientific, recovery=recovery, other_advancing=other_advancing,
                      export=run['final_export_seconds'], process_remainder=process_overhead,
                      unattributed_complete_seconds=unattributed,
                      attributed_complete_fraction=1-unattributed/receipt['seconds'],
                      nonlinear_iterations=[s['nonlinear_iterations'] for s in internal],
                      coupled_krylov_iterations=[row['krylov_iterations'] for row in newton],
                      coupled_Jacobian_actions=[row['jacobian_actions'] for row in newton],
                      coupled_residual_evaluations=[row['residual_evaluations'] for row in newton],
                      inverse_factorizations=[s['rebuilds'] for s in internal],
                      inverse_refreshes=[row['refreshes'] for row in newton],
                      inverse_reuses=[row['reuses'] for row in newton],
                      inverse_refresh_reasons=[row['refresh_reasons'] for row in newton],
                      rejected_trials=[row['rejected_trials'] for row in newton],
                      retries=[s['retries'] for s in internal], residual_max=max(max(s['residuals']) for s in internal),
                      files=[identity(p) for p in sorted(OUT.glob(case+'*')) if p.is_file() and
                             (p.name.startswith(case+'.') or p.name.startswith(case+'-') and p.suffix == '.fld')])
        block = max(1, steps//4)
        means = [statistics.mean(complete[i:i+block]) for i in range(block, steps-block+1, block)]
        largest = max(means)
        advancement = math.ceil(505.02466521278404/dt)*largest
        result.update(warm_block_size=block, warm_complete_seconds_per_step=means,
                      warm_seconds_per_unit_physical_time=largest/dt,
                      conditional_advancement_seconds=advancement, conditional_margin_seconds=1.25*advancement,
                      equation_or_solve_share_of_whole=sum(disjoint[k] for k in
                            ('constitutive', 'wall', 'assembly', 'flow', 'scalar', 'forcing', 'residual'))/receipt['seconds'])
        if case == 'paired':
            baseline = REFERENCE.with_suffix('.result.json')
            previous = read(baseline)['state']
            roles = ('state_sha256', 'previous_sha256', 'advection_sha256', 'previous_advection_sha256')
            result['reviewed_startup_binding'] = dict(previous=identity(baseline),
                history_roles={k: state[k] == previous[k] for k in roles},
                independent=identity(REFERENCE.with_name('physical-current.json')))
    comparisons = {f'{a}-to-{b}': identity(p) for a,b in COMPARISONS
                   if (p := OUT/f'{a}-to-{b}.json').is_file()}
    atomic(output, dict(schema='WJH27ShortCostAndStartupScreenV2', source_identity=SOURCE,
                                      inputs=identity(OUT/'inputs.json'), measurements=rows, comparisons=comparisons,
                                      assessor=identity(__file__),
                                      scope='Bounded retained-inlet startup and conventional cost comparison',
                                      selected_recipe=None, complete_trajectory_forecast=None, TASK3_READY=False))


def run():
    from physical_reference_refinement import ELF, H27, TABLE as TABLE_SHA, read, require
    from resources import identity, size

    account = read(ROOT/'accounting.json')
    settled = [v['seconds'] for v in account['entries'].values() if v['status'] != 'RESERVED']
    require(len(settled) == len(account['entries']) and not account['closed'] and
            sum(settled)+1705 <= account['ceiling']-account['evidence_protected'], 'cost screen cap needs an open settled account')
    require(size(ROOT/'science')+.3*1024**3 < 8*1024**3 and shutil.disk_usage(ROOT).free > 15.3*1024**3,
            'cost screen storage headroom differs')
    binary = ROOT/'runtime-v19/IncNavierStokesSolver'
    require(identity(binary.with_suffix('.elf'))['sha256'] == ELF and identity(TABLE)['sha256'] == TABLE_SHA,
            'wrong selected binary or table identity')
    py = [sys.executable, '-B']
    runner = [*py, ROOT/'tools/run.py', '--root', ROOT]

    def operation(label, phase, cap, *argv):
        return subprocess.run(list(map(str, [*runner, '--label', NAME+'-'+label, '--phase',
                                             phase, '--seconds', cap, '--', *argv]))).returncode == 0

    require(operation('fixtures', 'analysis', 30, *py, __file__, 'fixtures'), 'screen fixtures failed')
    passed = set()
    order = (STOCK[0], 'paired', *STOCK[1:], *(c for c in CASE_CAPS if c != 'paired'))
    for case in order:
        xml = OUT/(case+'.xml')
        timer = ['/usr/bin/time', '-f', 'H27 child elapsed seconds: %e'] if case in STOCK else []
        if operation(case, 'science', CASE_CAPS[case][-1] if case in CASE_CAPS else 45,
                     *timer, binary, '--use-opt-file', xml.with_suffix('.xml.opt'), xml):
            passed.add(case)
    for first, second in COMPARISONS:
        if {first, second} <= passed:
            operation(first+'-to-'+second, 'analysis', 160, *py, H27/'physical_transient_oracle.py',
                      '--native', OUT/(first+'.result.json'), '--compare-finer', OUT/(second+'.result.json'),
                      '--table', TABLE, '--out', OUT/f'{first}-to-{second}.json', '--seconds', 145)
    require(operation('assessment', 'analysis', 50, *py, __file__, 'assess'), 'cost assessment failed')
    print('SOURCE19_PHYSICAL_COST_SCREEN_FINISHED', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'fixtures', 'assess'))
    parser.add_argument('--assessment-out', type=Path, default=OUT/'assessment.json')
    args = parser.parse_args()
    if not sys.platform.startswith('linux') or not ROOT.is_dir():
        parser.error('requires the authoritative Linux H27 run root')
    sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27')]
    if args.action == 'assess':
        assess(args.assessment_out)
    else:
        {'run': run, 'fixtures': fixtures}[args.action]()


if __name__ == '__main__':
    main()
