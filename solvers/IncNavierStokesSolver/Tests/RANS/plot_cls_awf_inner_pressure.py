#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot already verified matched readbacks; never evaluate a constitutive law."""
import argparse
from pathlib import Path

from analyse_cls_awf_inner_pressure import verified_physical
from analyse_cls_accepted_awf import analytical_column, dense_column
from run_cls_epawf_picard import load, require, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    d = load(a.summary); old = {}
    for receipt in d['baseline_summaries']:
        file = Path(receipt['path']); require(sha(file) == receipt['sha256'], 'baseline summary differs')
        for label, row in load(file)['cases'].items():
            if label.startswith('refined-'): old[int(row['selection']['native_step'])] = row
    steps = sorted(map(int, d['observations'])); require(len(steps) > 1, 'evolved readback required')
    a.output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    for name, rows, color in [('P1 inner', old, 'tab:blue'), ('P2 inner', {int(k): v for k, v in d['observations'].items()}, 'tab:orange')]:
        times = [rows[n]['physical']['time'] for n in steps]
        axes[0, 0].plot(times, [100*(rows[n]['physical']['mean_wall_stress']/.16-1) for n in steps], 'o-', label=name, color=color)
        for region, style in [(0, '-'), (1, '--')]:
            axes[0, 1].plot(times, [rows[n]['complete_audit']['strong_continuity'][region] for n in steps],
                style, marker='o', color=color, label=name+(' inner' if region == 0 else ' outer'))
        axes[1, 0].plot(times, [max(rows[n]['physical']['common_normalized_spatial_first4']) for n in steps], 'o-', label=name, color=color)
        axes[1, 1].plot(times, [rows[n]['physical']['matching_trace_defect'] for n in steps], 'o-', label=name, color=color)
    axes[0, 0].set(ylabel='Mean traction reference error (%)')
    axes[0, 1].set(ylabel='Strong divergence physical L2', yscale='log')
    axes[1, 0].set(ylabel='Largest common-normalized spatial block', yscale='log')
    axes[1, 1].set(ylabel='Refreshed epsilon trace defect', yscale='log')
    for axis in axes.flat: axis.set_xlabel('Physical time'); axis.grid(alpha=.2); axis.legend(fontsize=8)
    fig.savefig(a.output/'pressure-space-evolution.png', dpi=150); plt.close(fig)
    n = steps[-1]; last = d['observations'][str(n)]; q = verified_physical(last); q1 = verified_physical(old[n]); q0 = verified_physical(d['observations']['0'])
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout='constrained')
    for i, station in enumerate([-3.5, 1.6]):
        for name, field, color, style in [('native h0', q0, 'gray', ':'), ('P1 inner', q1, 'tab:blue', '-'), ('P2 inner', q, 'tab:orange', '-')]:
            x, column = dense_column(field, station); z = np.array(column['z']); ref = analytical_column(field, z)
            values = [np.array(column['u'])-ref['u'], np.array(column['k'])-ref['k'],
                      100*(np.array(column['du_dz'])/np.array(ref['du_dz'])-1)]
            for j, value in enumerate(values): axes[i, j].plot(z, value, style, color=color, label=name)
        for j, label in enumerate(['U minus analytical reference', 'k minus analytical reference', 'dU/dz reference error (%)']):
            axes[i, j].set(xlabel='Physical z', ylabel=label, xscale='symlog', xlim=(0, .5), title=f'Native column x={x:.6f}, T={q["time"]:.6g}')
            axes[i, j].set_xscale('symlog', linthresh=.0025); axes[i, j].grid(alpha=.2); axes[i, j].axhline(0, color='k', lw=.5)
            axes[i, j].legend(fontsize=8)
    fig.savefig(a.output/'pressure-space-profiles.png', dpi=150); plt.close(fig)


if __name__ == '__main__': main()
