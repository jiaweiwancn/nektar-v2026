#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot immutable common-panel fields and native one-sided ground traction."""
import argparse
from pathlib import Path
import os

from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import THREADS, load, require, sha


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--diagnosis', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--labels', nargs='+', required=True)
    p.add_argument('--station', type=float, default=1.6); a = p.parse_args(); os.environ.update(THREADS)
    import numpy as np
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    require(not a.output.exists(), 'preserve previous plot generations')
    source = load(a.diagnosis/'summary.json'); fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout='constrained')
    receipts = []; reference = False
    for label in a.labels:
        row = source['cases'][label]; result = Path(row['receipt']['output'])
        require(sha(result) == row['receipt']['output_sha256'], 'plot input differs from retained native result')
        physical = load(result)['diagnostics']['physical']; s = physical['common_samples']
        x = np.array(s['x']); z = np.array(s['z']); station = min(x, key=lambda v: abs(v-a.station)); mask = (x == station) & (z<=physical['H'])
        for axis, key, ref in [(axes[0, 0], 'u', 'reference_u'), (axes[0, 1], 'du_dz', 'reference_du_dz'), (axes[0, 2], 'k', 'reference_k')]:
            axis.plot(np.array(s[key])[mask], z[mask], label=label)
            if not reference: axis.plot(np.array(s[ref])[mask], z[mask], 'k:', label='analytical')
        ez = z[z >= physical['H']]; ex = x[z >= physical['H']]; emask = ex == station
        axes[1, 0].plot(np.array(s['epsilon'])[emask], ez[emask], label=label)
        if not reference: axes[1, 0].plot(np.array(s['reference_epsilon'])[emask], ez[emask], 'k:', label='analytical')
        reference = True
        wall = row.get('native_ground_traction')
        if wall: axes[1, 1].plot([q['coordinate'][0] for q in wall], [q['traction'] for q in wall], label=label)
        axes[1, 2].semilogy(range(4), physical['common_normalized_spatial_first4'], 'o-', label=label)
        receipts.append(dict(label=label, result=str(result), sha256=sha(result), actual_station=float(station)))
    for ax, xlabel in zip(axes[0], ('U', 'dU/dz', 'k')):
        ax.set(xlabel=xlabel, ylabel='physical z', ylim=(0, .5)); ax.set_yscale('symlog', linthresh=.0025)
        sec = ax.twinx(); sec.set_yscale('symlog',linthresh=.005);sec.set_ylim(0,1);sec.set_yticks([0,.01,.1,1]);sec.set_ylabel('z/H')
    axes[0, 1].set_xscale('log')
    axes[1, 0].set(xlabel='transported outer epsilon', ylabel='physical z', ylim=(.5, 4))
    axes[1, 1].axhline(.16, color='black', linestyle=':'); axes[1, 1].set(xlabel='native ground x', ylabel='wall traction')
    axes[1, 2].axhline(1e-3, color='black', linestyle=':'); axes[1, 2].set_xticks(range(4), ['Ux', 'W', 'k', 'epsilon'])
    axes[1, 2].set(ylabel='common-normalized selected spatial residual')
    axes[0, 0].legend(fontsize=7); fig.suptitle(f'H=.5 calibrated budget AWF, x={station:.6g}: reference error and evolution')
    a.output.parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.output, dpi=170); plt.close(fig)
    atomic(a.output.with_suffix('.json'), dict(script_sha256=sha(Path(__file__)), inputs=receipts, plot_sha256=sha(a.output),
        governing_solves=0, state_or_history_changes=False))


if __name__ == '__main__': main()
