#!/usr/bin/env python3
"""Plot the measured standard-AWF preservation and refinement CSVs on Linux."""
import argparse
import csv
import json
from pathlib import Path
import sys


def rows(path):
    return list(csv.DictReader(path.open()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not sys.platform.startswith('linux'):
        parser.error('Generate scientific results in the authoritative Linux environment')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    report = json.loads((args.study/'validation-study.json').read_text())
    if report['status'] != 'passed':
        parser.error('Use a completed study report')
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .25, 'savefig.dpi': 160})
    fields = [('u', 'Velocity'), ('k', 'Turbulent kinetic energy'), ('eps', 'Dissipation'), ('rans_nut', 'Eddy viscosity')]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharex=True)
    for dim in (2, 3):
        for mode, color, style in [('fixed', '#d17c1f', '-'), ('switched', '#127569', '--')]:
            name = f'five_d{dim}_{mode}'
            data = rows(args.study/name/(name+'.rans-errors.csv'))
            spec = report['cases'][name]['spec']
            for column, (field, title) in enumerate(fields):
                group = [r for r in data if r['field'] == field]
                axes[dim-2, column].semilogy([float(r['time'])*spec['uref']/spec['length'] for r in group],
                    [max(1e-10, 100*float(r['normalized_l2'])) for r in group],
                    color=color, ls=style, label=mode.capitalize())
        for column, (_, title) in enumerate(fields):
            ax = axes[dim-2, column]
            ax.set_title(f'{dim}D {title}')
            ax.axhline(.1, color='black', lw=.8, ls=':', label='0.1% target')
            ax.set_ylim(1e-8, .3)
            if column == 0:
                ax.set_ylabel('Relative L2 error (%)')
            if dim == 3:
                ax.set_xlabel('Flow-through times')
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output/'preservation.png')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for field, title in fields:
        values = report['spatial_refinement'][field]
        axes[0].loglog([6, 12, 24], values['h_errors'], 'o-', label=title)
        axes[1].semilogy([3, 4, 5], values['p_errors'], 'o-', label=title)
    axes[0].set_xlabel('Vertical layers (degree 4)')
    axes[1].set_xlabel('Polynomial degree (12 layers)')
    axes[0].set_xticks([6, 12, 24], ['6', '12', '24'])
    axes[1].set_xticks([3, 4, 5])
    for ax in axes:
        ax.set_ylabel('Final relative L2 error')
        ax.set_title('2D; fixed H=2 m, P=1 m, N=4 m')
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output/'refinement.png')
    plt.close(fig)
    fig, axes = plt.subplots(1, 4, figsize=(13, 4), sharey=True)
    for dim, style, color in [(2, '-', '#127569'), (3, '--', '#7850a2')]:
        name = f'five_d{dim}_switched'
        data = rows(args.study/name/'profiles.csv')
        for ax, (field, title) in zip(axes, fields):
            group = [r for r in data if r['field'] == field and float(r['x']) == 285]
            ax.plot([100*(float(r['value'])-float(r['reference']))/float(r['reference']) for r in group],
                    [float(r['height']) for r in group], ls=style, color=color, label=f'{dim}D')
            ax.set_title(title)
            ax.set_xlabel('Profile error (%)')
            ax.xaxis.set_major_locator(MaxNLocator(4))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((0, 0))
            ax.xaxis.set_major_formatter(formatter)
    axes[0].set_ylabel('Physical height (m)')
    axes[0].legend()
    fig.suptitle('Switched standard AWF after five flows; x/L=0.95')
    fig.tight_layout()
    fig.savefig(args.output/'profiles.png')
    plt.close(fig)
    print('Created preservation.png, refinement.png and profiles.png')


if __name__ == '__main__':
    main()
