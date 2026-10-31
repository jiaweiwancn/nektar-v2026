#!/usr/bin/env python3
"""Plot recorded terrain verification and short transition-width pilots."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wg2',type=Path,required=True)
    parser.add_argument('--wg4',type=Path,help='Also plot completed transition-width pilots')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    wg2 = json.loads(args.wg2.read_text())['suites']
    args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,2,figsize=(8,3.1),constrained_layout=True)
    for model,label,colour in [('standard','Standard log wall','#1769aa'),
                               ('standard_fixed','Standard AWF','#bd3d32')]:
        for suite,style,geometry in [('curved','--','compact C1'),('periodic-curved','-','smooth periodic')]:
            cases = wg2[suite]['cases']
            for ax,metric,title in zip(axes,('normal_pointwise_max','normal_l2'),
                                       ('Global maximum','Global strong L2')):
                ax.loglog([1,.5,.25],[cases[model+'_h'+str(n)][metric] for n in (1,2,4)],
                          style+'o',color=colour,label=label+', '+geometry)
                ax.set(xlabel='Relative mesh size',ylabel='Normal-condition error (m/s)',title=title)
                ax.set_xticks([.25,.5,1],['1/4','1/2','1'])
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.grid(True,which='both',alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle('Normal-condition refinement with integrated strip continuity\n'
                 'One step per grid: log law dt=2e-4 s; AWF dt=2e-5 s',fontsize=10)
    for suffix in ('png','pdf'): fig.savefig(args.output/('normal-refinement.'+suffix),dpi=180)
    plt.close(fig)
    if not args.wg4:
        return
    suites = json.loads(args.wg4.read_text())['suites']
    pilot = suites.get('pilot',suites.get('pilot-2d'))
    if not pilot: raise RuntimeError('No completed hill/width pilot evidence')
    fig,axes = plt.subplots(1,2,figsize=(8,3.1),constrained_layout=True)
    for width in (1,2,4,8):
        start,end=-.2-width*.04,-.2
        xs=[-.56+i*.5/400 for i in range(401)]
        rs=[max(0.,min(1.,(x-start)/(end-start))) for x in xs]
        axes[0].plot(xs,[1-10*r**3+15*r**4-6*r**5 for r in rs],label=str(width))
    axes[0].axvline(-.1,color='.4',linestyle=':',label='Hill foot')
    axes[0].set(xlabel='Physical wall foot x (m)',ylabel='Prescribed correction mask',
                title='Upstream C2 correction mask\nTransition length / hill height')
    axes[0].legend(title='Length / height',fontsize=7)
    for closure,colour in [('standard','#1769aa'),('cls','#bd3d32')]:
        study=pilot['comparisons'][closure+'_width_study']
        for terrain,style in [('hill','-'),('floor','--')]:
            baseline = study[0][terrain]['incident_tau']
            # All supplied pilots use u_tau=0.3 m/s. A common friction-velocity
            # scale shows the small width response in both hill/floor controls.
            axes[1].plot([r['width_over_hill_height'] for r in study],
                [1e6*(r[terrain]['incident_tau']-baseline)/.3**2 for r in study],style+'o',color=colour,
                label=('Standard' if closure == 'standard' else 'CLS')+', '+terrain)
    axes[1].set(xlabel='Transition length / hill height',ylabel=r'$10^6\Delta(\tau/\rho)/u_\tau^2$ from width 1',
                title='Two-step startup (t=4e-5 s)\nNo optimal width inferred')
    axes[1].set_xticks([1,2,4,8]);axes[1].legend(fontsize=7)
    for ax in axes:ax.grid(True,alpha=.2)
    for suffix in ('png','pdf'):fig.savefig(args.output/('transition-widths.'+suffix),dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
