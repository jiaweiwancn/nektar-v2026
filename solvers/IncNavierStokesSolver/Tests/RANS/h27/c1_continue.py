#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Continue a physical C1 hill using independently checked implicit states."""

import argparse
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from resources import identity


ROOT = Path(__file__).resolve().parent
REGIONS = ('inner', 'outer', 'matching_inner', 'matching_outer')


def read(path):
    return json.loads(path.read_text())


def accepted_state(path):
    packet = read(path)
    state = packet['state']
    if (not packet['normal_main'] or state['manufactured'] or
            max(state['residuals']) > 1e-10):
        raise ValueError(f'Physical implicit residual acceptance missing in {path}')
    if packet['steady_mode']:
        if not packet['steady_solver']['converged'] or state['step'] != 0:
            raise ValueError(f'Converged steady state missing in {path}')
    elif (state['step'] < 2 or state['step_advection'] != 'BDF2ImplicitAdvectionV1' or
          state['accepted_by_fixed_budget']):
        raise ValueError(f'Residual-controlled implicit BDF2 state missing in {path}')
    return state


def checked(result, residual_file, margin_file, table):
    audit, margins = read(residual_file), read(margin_file)
    native_sha = identity(result)['sha256']
    for check in (audit, margins):
        if check['native']['sha256'] != native_sha or check['table']['sha256'] != table['sha256']:
            raise ValueError(f'Independent evidence binding differs for {result}')
    return (audit['residual_pass'] and audit['epsilon_trace']['interior_projection_defect'] <= 1e-10
            and all(margins['regions'][key]['sample_pass'] for key in REGIONS))


def run(args, cwd, name):
    with (cwd/name).open('x') as log:
        return subprocess.run([str(x) for x in args], cwd=cwd, stdin=subprocess.DEVNULL,
                              stdout=log, stderr=subprocess.STDOUT, check=False).returncode


def fail_report(case):
    failure = case/'state.failure.json'
    if not failure.is_file():
        return {}
    doc = read(failure)
    return dict(failure=identity(failure), reason=doc['reason'],
                last_accepted_step=doc['state']['step'], last_accepted_time=doc['state']['time'],
                failed_residuals=doc['failed_attempt']['residuals'])


def advance(args, directory, height, parent_path, candidate, nx, nr):
    setup = [sys.executable, args.fixture, 'hill', directory, '--height', candidate,
             '--from-height', height, '--initial-fields', parent_path, '--quadrature', nx, nr,
             '--dt', args.dt, '--steps', args.steps, '--step-seconds', args.step_seconds,
             '--science', args.steps, '--recovery', 1]
    if args.seed_taper is not None:
        setup.extend(('--seed-taper', args.seed_taper))
    if run(setup, args.root, f'{directory.name}-fixture.log'):
        raise RuntimeError(f'Fixture generation failed for {directory}')
    began = time.monotonic()
    status = run(['/usr/bin/time', '-v', 'mpirun', '--stdin', 'none', '--bind-to', 'none', '-np', 1,
                  args.solver, '--use-opt-file', 'run.xml.opt', 'run.xml'], directory, 'run.log')
    return dict(status=status, wall_seconds=time.monotonic()-began,
                native_messages=[line for line in (directory/'run.log').read_text().splitlines()
                                 if 'WJ-H27' in line])


def audit_stage(args, directory, native, q):
    independent = directory/'independent.json'
    margins = directory/'independent-margins.json'
    for program, output, log in (('c1_oracle.py', independent, 'audit.log'),
                                 ('c1_margins.py', margins, 'margins.log')):
        commands = [sys.executable, ROOT/program, '--table', args.table, '--native', native,
                    '--order', q, '--seconds', args.audit_seconds, '--out', output]
        if run(commands, directory, log):
            raise RuntimeError(f'{program} failed for {directory}')
    return independent, margins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('solver', 'table', 'start', 'start-audit', 'start-margins', 'root'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--fixture', type=Path, default=ROOT/'c1_fixtures.py')
    parser.add_argument('--increment', type=Decimal, default=Decimal('.00315'))
    parser.add_argument('--minimum', type=Decimal, default=Decimal('.0000625'))
    parser.add_argument('--maximum', type=Decimal, default=Decimal('.075'))
    parser.add_argument('--target', type=Decimal, default=Decimal('1'))
    parser.add_argument('--attempts', type=int, default=50)
    parser.add_argument('--dt', type=float, default=.02)
    parser.add_argument('--steps', type=int, default=2)
    parser.add_argument('--step-seconds', type=float, default=180.)
    parser.add_argument('--audit-seconds', type=float, default=300.)
    parser.add_argument('--quadrature', nargs=2, type=int, default=[10, 16])
    parser.add_argument('--seed-taper', type=float)
    args = parser.parse_args()
    for option in ('solver', 'table', 'start', 'start_audit', 'start_margins', 'root', 'fixture'):
        setattr(args, option, getattr(args, option).resolve())
    if args.steps < 2 or args.attempts < 1 or not (0 < args.minimum <= args.increment <= args.maximum):
        parser.error('Require BE plus BDF2, a positive attempt cap and ordered positive height increments')
    parent = args.start
    preceding = accepted_state(parent)
    parent_audit, parent_margins = args.start_audit, args.start_margins
    table = identity(args.table)
    if not checked(parent, parent_audit, parent_margins, table):
        raise ValueError('Starting state lacks a passing independent check')
    height = Decimal(str(preceding['geometry']['hill_height']))
    if not 0 < height < args.target <= 1:
        parser.error('Require 0 < starting height < target <= 1')
    args.root.mkdir()
    increment = args.increment
    nx, nr = args.quadrature
    with (args.root/'stages.jsonl').open('x') as evidence:
        for attempt in range(1, args.attempts+1):
            candidate = min(args.target, height+increment)
            case = args.root/f'{attempt:02d}-height-{candidate}-q{nx}x{nr}'
            result = advance(args, case, height, parent, candidate, nx, nr)
            row = dict(attempt=attempt, height=float(candidate), source_height=float(height),
                       quad=[nx, nr], directory=str(case), seed=identity(parent),
                       source_time=preceding['time'], **result)
            if result['status']:
                row.update(decision='failed-native', **fail_report(case))
                increment /= 2
            else:
                native = case/'state.result.json'
                state = accepted_state(native)
                if (state['geometry_seed']['source_step'] != preceding['step'] or
                        state['geometry_seed']['source_time'] != preceding['time'] or
                        state['step'] != args.steps):
                    raise RuntimeError('Geometry transfer must build fresh time histories')
                independent, margins = audit_stage(args, case, native, max(20, nx+4, nr+4))
                outer = read(margins)['regions']['outer']
                row.update(result=identity(native), residuals=state['residuals'],
                           covariance=state['minimum_covariance_over_k'],
                           principal=state['minimum_principal_over_mu'],
                           whole_K_over_Kin=state['interface_K_over_Kin'],
                           K_margin=state['interface_K_table_margin_over_Kin'],
                           cutoff=state['interface_cutoff_range'], minimum_k=state['minimum_k'],
                           minimum_epsilon=state['minimum_epsilon'], audit=identity(independent),
                           margins=identity(margins), independent_outer_principal=outer['principal_over_mu'])
                if checked(native, independent, margins, table):
                    row['decision'] = 'accepted'
                    height, parent, preceding = candidate, native, state
                    parent_audit, parent_margins = independent, margins
                    if outer['principal_over_mu'] > 5e-4:
                        increment = min(args.maximum, increment*Decimal('1.6'))
                else:
                    row['decision'] = 'failed-independent-check'
                    if nx < 18:
                        nx, nr = nx+2, nr+4
                    else:
                        increment /= 2
            evidence.write(json.dumps(row, allow_nan=False)+'\n')
            evidence.flush()
            os.fsync(evidence.fileno())
            print(json.dumps(row, allow_nan=False), flush=True)
            if height == args.target or increment < args.minimum:
                break
    summary = dict(height=float(height), target=float(args.target), attempts=attempt,
                   complete_height_ladder=height == args.target,
                   last_accepted_state=identity(parent), audit=identity(parent_audit),
                   margins=identity(parent_margins), positive_time_at_last_height=preceding['time'])
    (args.root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
