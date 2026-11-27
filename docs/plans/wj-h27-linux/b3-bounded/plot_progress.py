#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot the matched finite trials and subsequently accepted B3 lineages."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(evidence, output):
    seeds = [t for t in evidence["frozen_trials"]["trials"] if t["metrics"]["fraction"] == 1]
    proposed = [t["metrics"]["physical_model"]*1.e6 for t in seeds]
    actual = [t["metrics"]["physical_actual"]*1.e6 for t in seeds]
    figure, (left, right) = plt.subplots(1,2,figsize=(10.,4.3),constrained_layout=True)
    frozen = max(evidence["frozen_trials"]["base_residuals"])*1.e6
    x = np.arange(len(proposed))
    left.bar(x-.17,proposed,.34,color="#c4dce6",label="Current linear model")
    left.bar(x+.17,actual,.34,color="#177088",label="Fresh original residual")
    left.axhline(frozen,color="#80726d",ls="--",label="Frozen base")
    left.set_xticks(x,[".005",".025",".1","Old ray\n.025"])
    left.set_xlabel("Full trial direction and scaled radius")
    left.set_ylabel(r"Original physical maximum ($10^{-6}$)")
    left.set_title("Matched finite trials",loc="left")
    left.legend(frameon=False,fontsize=8)

    styles = (("bounded-main","Current Krylov",0,"#264e79","o"),
              ("bounded-recycled","Recycled span",4,"#92793f","s"),
              ("bounded-enriched","Scalar-enriched span",4,"#298272","o"),
              ("bounded-active","Active-gradient seed",6,"#b45039","D"))
    for case,label,offset,color,mark in styles:
        if case not in evidence["lineages"]:
            continue
        d = evidence["lineages"][case]["run"]
        values = [d["initial_max"]]+[max(a["residuals"]) for a in d["accepted"]]
        right.plot(np.arange(len(values))+offset,np.asarray(values)*1.e6,marker=mark,
                   label=label,color=color,lw=1.3,ms=4)
    right.set_title("Accepted original progress",loc="left")
    right.set_xlabel("Accepted update along each lineage")
    right.legend(frameon=False,fontsize=8)
    for axes in (left,right):
        axes.grid(axis="y",alpha=.2)
    figure.savefig(output,dpi=175)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    plot(json.loads(args.evidence.read_text()),args.output)
