#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific figures for the bounded, unchanged WJ-6 finite-strip problem."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compact', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    load = lambda name: list(csv.DictReader((args.compact / (name + '.csv')).open()))
    colors = dict(coarse='#315c8c', fine='#b16a28', baseline='#43834d')
    labels = dict(coarse='L2/Nx2, dt=32 us', fine='L2/Nx2, dt=16 us', baseline='L2/Nx1, dt=32 us')

    def save(fig, name):
        fig.tight_layout()
        for suffix in ('png', 'pdf', 'svg'):
            fig.savefig(out / (name + '.' + suffix), dpi=160, bbox_inches='tight')
        plt.close(fig)

    divergence = load('divergence')
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, policy in zip(axes, ('EP', 'Budget')):
        for tag, color in colors.items():
            for scope, style in (('full', '-'), ('central', '--')):
                q = sorted((r for r in divergence if r['policy'] == policy and r['tag'] == tag and
                            r['scope'] == scope and r['field'] == 'divergence'), key=lambda r: float(r['time']))
                ax.plot([float(r['time']) for r in q], [2*float(r['sampled_Linf']) for r in q],
                        style, marker='.', color=color, label=labels[tag] + ' ' + scope)
        ax.axhline(.001, color='#a33232', ls=':', label='unchanged budget')
        ax.set_title(policy)
        ax.set_xlabel('physical time [s]')
        ax.set_ylabel(r'sampled $|\mathrm{div}\,U|H/u_*$')
        ax.grid(alpha=.2)
        ax.legend(fontsize=6)
    fig.suptitle('WJ-6 uniform refinement; dimensionless time = 0.5 t [s]')
    save(fig, 'strong-divergence')

    fields = load('fields')
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
    for ax, field in zip(axes.flat, ('velocity', 'k', 'outer_epsilon', 'mechanical_pressure')):
        for policy, marker in (('EP', 'o'), ('Budget', 'x')):
            for tag, color in colors.items():
                for kind, style in (('continuum_error', '-'), ('initial_drift', '--')):
                    q = sorted((r for r in fields if r['policy'] == policy and r['tag'] == tag and
                                r['scope'] == 'full' and r['field'] == field and r['kind'] == kind),
                               key=lambda r: float(r['time']))
                    ax.semilogy([float(r['time']) for r in q],
                                [max(float(r['normalized_L2']), 1e-18) for r in q], style,
                                color=color, marker=marker, ms=3, label=policy+' '+tag+' '+kind)
        ax.axhline(1e-4, color='#a33232', ls=':')
        ax.set_title(field)
        ax.set_ylabel('normalized physical L2')
        ax.set_xlabel('physical time [s]')
        ax.grid(alpha=.2)
        ax.legend(fontsize=5, ncol=2)
    fig.suptitle('Total piecewise-reference error and represented-initial drift')
    save(fig, 'reference-error-and-drift')

    sensitivity = load('sensitivity')
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, group in zip(axes, ('horizontal', 'timestep')):
        for policy, color in (('EP', '#315c8c'), ('Budget', '#b16a28')):
            q = [r for r in sensitivity if r['policy'] == policy and r['phase'] == 'primary' and
                 r['quantity'].startswith(group+'_sensitivity/')]
            times = sorted({int(r['quantity'].split('/')[1]) for r in q})
            ax.semilogy([s*3.2e-5 for s in times],
                        [max(float(r['fraction']) for r in q if int(r['quantity'].split('/')[1]) == s) for s in times],
                        'o-', color=color, label=policy)
        ax.axhline(1, color='#a33232', ls='--')
        ax.set_title(group+' sensitivity: worst required field/scope')
        ax.set_ylabel('(difference + both uncertainty contributions) / limit')
        ax.set_xlabel('common physical time [s]')
        ax.grid(alpha=.2)
        ax.legend()
    save(fig, 'matched-sensitivity')

    projection = load('pressure-projection')
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, policy in zip(axes, ('EP', 'Budget')):
        for tag, color in colors.items():
            for key, style in (('Q1_projection_L2', '--'), ('complement_L2', '-')):
                q = sorted((r for r in projection if r['policy'] == policy and r['tag'] == tag and
                            float(r['region']) == 0), key=lambda r: float(r['time']))
                ax.semilogy([float(r['time']) for r in q], [max(float(r[key]), 1e-25) for r in q],
                            style, color=color, label=tag+' '+key)
        ax.set_title(policy+' inner-region divergence')
        ax.set_ylabel('physical L2')
        ax.set_xlabel('physical time [s]')
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    save(fig, 'pressure-test-projection')

    interfaces = load('interfaces')
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    # These table fields already use their declared physical amplitudes.
    for ax, (name, budget) in zip(axes, [('traction', 1e-3), ('k_flux', 1e-4), ('pressure_offset', 1e-4)]):
        for policy, marker in (('EP', 'o'), ('Budget', 'x')):
            for tag, color in colors.items():
                q = sorted((r for r in interfaces if r['policy'] == policy and r['tag'] == tag), key=lambda r: float(r['time']))
                ax.semilogy([float(r['time']) for r in q], [max(abs(float(r[name])), 1e-18) for r in q],
                            marker+'-', ms=3, color=color, label=policy+' '+tag)
        ax.axhline(budget, color='#a33232', ls=':')
        ax.set_title(name.replace('_', ' '))
        ax.set_ylabel('normalized sampled maximum')
        ax.set_xlabel('physical time [s]')
        ax.grid(alpha=.2)
    axes[-1].legend(fontsize=6)
    save(fig, 'interface-diagnostics')
    manifest = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in out.iterdir() if f.is_file()}
    (out/'figures.json').write_text(json.dumps(dict(schema='WJ6ScientificFiguresV1', files=manifest,
        scope='serial flat2D finite reference-Dirichlet strip; residuals and actual-stage defects remain separate tables'), indent=2)+'\n')


if __name__ == '__main__':
    main()
