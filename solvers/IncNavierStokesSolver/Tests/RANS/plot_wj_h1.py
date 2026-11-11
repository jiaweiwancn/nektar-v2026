#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Static scientific H1 figures from immutable compact physical observations."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(); p.add_argument('data', type=Path); p.add_argument('output', type=Path)
    args = p.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'figure.dpi': 140, 'savefig.dpi': 170})
    s = np.linspace(-5, 5, 1001); f = np.where(abs(s) < 2.5, np.cos(np.pi * s / 5) ** 2, 0)
    fp = np.where(abs(s) < 2.5, -np.pi / 5 * np.sin(2 * np.pi * s / 5), 0); g = np.hypot(1, fp)
    def curve(r): return s - (.0875 + r) * fp / g, f + (.0875 + r) / g
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for a in ax:
        a.plot(s, f, 'k--', label='physical hill')
        for r, label in [(0, 'effective ground'), (.375, 'matching surface')]: a.plot(*curve(r), label=label)
        for foot in [-5, -2.5, 0, 2.5, 5]:
            j = np.argmin(abs(s - foot)); x, z = curve(.375)
            a.plot([s[j], x[j], foot], [f[j], z[j], 20], color='.65', lw=.7)
        a.set(xlabel=r'$x/H_h$', ylabel=r'$z/H_h$'); a.grid(alpha=.2)
    ax[0].set_ylim(0, 20.4); ax[0].set_title('Certified composite map: full tunnel')
    ax[1].set_ylim(0, 1.8); ax[1].set_title('Normal inner region; independent outer continuation')
    ax[1].legend(loc='upper left', fontsize=8)
    fig.tight_layout(); fig.savefig(args.output / 'composite-geometry.png'); plt.close(fig)
    rows = list(csv.DictReader((args.data / 'hill-field-samples.csv').open()))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
    for a, field, label in zip(ax, ['u', 'k'], [r'$U_x/U_\infty$', r'$k/U_\infty^2$']):
        selected = [r for r in rows if r['phase'] == 'final' and float(r['z']) < 2]
        dots = a.scatter([float(r['x']) for r in selected], [float(r['z']) for r in selected],
                         c=[float(r[field]) for r in selected], s=5, cmap='viridis', rasterized=True)
        a.plot(*curve(0), 'k', lw=.7); a.plot(*curve(.375), 'r--', lw=.7)
        a.set(xlabel=r'$x/H_h$', ylabel=r'$z/H_h$', ylim=(0, 2), title=label)
        fig.colorbar(dots, ax=a, shrink=.85)
    fig.suptitle('Controlled native Budget prototype, 3 accepted steps, '
                 r'$t U_\infty/H_h=3\times10^{-6}$; red: common interface', fontsize=10)
    fig.tight_layout(); fig.savefig(args.output / 'native-coupled-fields.png'); plt.close(fig)
    rows = list(csv.DictReader((args.data / 'interfaces.csv').open()))
    rows = sorted([r for r in rows if r['case'] == 'hill' and int(r['step']) == 3], key=lambda r: float(r['s']))
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.8))
    for a, names, title in zip(axes, [('K',), ('cutoff',), ('traction_x', 'traction_z', 'k_flux_sum')],
                                  [r'Endpoint $K-K_0$', 'Metric cutoff / matching height', 'One-sided pointwise defects']):
        for name in names:
            values = [float(r[name]) / (.375 if name == 'cutoff' else 1) for r in rows]
            if name == 'K': values = np.array(values) - 81 / (2 * np.sqrt(143)) * .051 ** 2
            a.plot([float(r['x']) for r in rows], values,
                   label=name, lw=1)
        a.set(xlabel=r'$x/H_h$', title=title); a.grid(alpha=.25); a.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))
        if len(names) > 1: a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(args.output / 'interface-diagnostics.png'); plt.close(fig)


if __name__ == '__main__': main()
