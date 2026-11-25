#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot the one-sided support profiles and compare the declared fixed strip."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

from c1_support_witness import FIELD_NAMES, choose
from resources import identity


ORDER = ("M0", "Mr", "Mx", "Mxr", "Mshift", "Mref", "Mshift-ref")
STYLES = (("#24384c","-"), ("#288e84",":"), ("#a05584","--"),
          ("#9b3330","-"), ("#bc7322","--"), ("#366bc5","-"), ("#536b22","-."))
PANELS = (("margin", r"Principal response / $\mu$"),
          ("Gxz", r"$G_{xz}$"), ("divergence", r"$\mathrm{tr}(G)$"),
          ("pressure", "Sampled pressure"), ("k_load", "Sampled k"),
          ("convection", "Sampled convection"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("profiles",type=Path)
    p.add_argument("image",type=Path)
    p.add_argument("--summary",type=Path,required=True)
    p.add_argument("--half-length",type=float,required=True)
    p.add_argument("--pair",nargs=2,action="append",required=True)
    args=p.parse_args()
    curves=defaultdict(lambda:defaultdict(list))
    strip=defaultdict(lambda:defaultdict(list))
    L=args.half_length
    for k in args.pair:
        if len(set(k))<2:
            p.error("comparison requires distinct field labels")
    compare_fields={name for pair in args.pair for name in pair}
    with args.profiles.open() as f:
        for row in csv.DictReader(f):
            if row["fixed"]!="True":
                continue
            name,a,r,offset=row["mesh"],float(row["a"]),float(row["r"]),float(row["r_minus_h"])
            values={k:float(row[k]) for k in FIELD_NAMES}
            if name in ORDER and offset==0.:
                curves[name][row["owner"]].append((a,values))
            if name in compare_fields and .35*L<=a<=L and abs(offset)<=.0015+1e-15:
                strip[name][a,r].append(dict(axial_side=int(row["axial_side"]),
                    radial_side=int(row["radial_side"]), **values))
    result={}
    for a,b in args.pair:
        if not strip[a] or strip[a].keys()!=strip[b].keys():
            raise ValueError("missing common coordinates in declared support strip")
        changes={k:0. for k in FIELD_NAMES}
        signs=both_adverse=0
        minima=[float("inf"),float("inf")]
        for xy,rows in strip[a].items():
            for sx,sr in ((-1,-1),(-1,1),(1,-1),(1,1)):
                u,v=choose(rows,sx,sr),choose(strip[b][xy],sx,sr)
                for k in FIELD_NAMES:
                    changes[k]=max(changes[k],abs(u[k]-v[k]))
                minima=[min(minima[0],u["margin"]),min(minima[1],v["margin"])]
                signs+=int((u["margin"]<0.)!=(v["margin"]<0.))
                both_adverse+=int(max(u["margin"],v["margin"])<0.)
        result[f"{a}/{b}"]=dict(maximum_absolute_difference=changes, smallest_margins=minima,
            sign_changes=signs,common_adverse_approaches=both_adverse,approaches=4*len(strip[a]))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,sharex=True,figsize=(10,8))
    handles=[]
    for name,(color,style) in zip(ORDER,STYLES):
        for axis,(key,unit) in zip(axes.ravel(),PANELS):
            for j,points in enumerate(curves[name].values()):
                points.sort(key=lambda row:row[0])
                line,=axis.plot([a for a,_ in points],[v[key] for _,v in points],
                                color=color,linestyle=style,lw=1.4,
                                marker="." if len(points)==1 else None,label=name if j==0 else None)
                if j==0 and key=="margin":
                    handles.append(line)
            axis.set_ylabel(unit)
            axis.grid(alpha=.23)
            axis.ticklabel_format(axis="y",style="sci",scilimits=(-2,2))
            axis.spines[["top","right"]].set_visible(False)
    axes[0,0].axhline(0.,color=".3",lw=.6)
    for ax in axes[-1,:]:
        ax.set_xlabel("Chart abscissa a (code units)")
    fig.legend(handles=handles,loc="upper center",ncol=7,frameon=False,fontsize=9)
    fig.tight_layout(rect=(0,0,1,.94))
    fig.savefig(args.image,dpi=175)
    plt.close(fig)
    with args.summary.open("x") as f:
        json.dump(dict(profile=identity(args.profiles),image=identity(args.image),
                  a_interval=[.35*L,L],absolute_radial_offset=.0015,results=result),
                  f,indent=2,allow_nan=False)
        f.write("\n")


if __name__=="__main__":
    main()
