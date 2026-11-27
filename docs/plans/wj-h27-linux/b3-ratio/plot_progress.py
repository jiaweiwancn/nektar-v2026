#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot fixed-base scalar-span measurements and durable original residual progress."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(data, output):
    figure, (left,right) = plt.subplots(1,2,figsize=(10.5,4.1),layout="constrained")
    trials = data["lineages"].get("scalars",{}).get("trials",[])
    panel = [t for t in trials if t["base"]==0 and t["span"] is not None and
             t["proposal"]<2 and "physical_actual" in t["measurements"]]
    if panel:
        model = [100*(m["physical_base"]-m["physical_model"])/m["physical_base"]
                 for m in (t["measurements"] for t in panel)]
        actual = [100*(m["physical_base"]-m["physical_actual"])/m["physical_base"]
                  for m in (t["measurements"] for t in panel)]
        x = np.arange(len(panel))
        left.bar(x-.17,model,.34,color="#93b9cb",label="Current model")
        left.bar(x+.17,actual,.34,color="#196d87",label="Actual trial")
        if len(panel)==8:
            left.set_xticks([.5,2.5,4.5,6.5],["Retained","k","Outer ε","Joint"])
        left.axhline(0,color="#666666",lw=.9)
        left.set_xlabel("Scalar residual span, radii .00625 and .025")
        left.set_ylabel("Decrease from matched base (%)")
        left.set_title("Matched original-residual trials",loc="left")
        left.legend(frameon=False,fontsize=8)
    else:
        trials = data["frozen"]["refined"]["trials"]
        x = np.arange(len(trials))
        rhos = [t["metrics"]["physical_decrease_ratio"] for t in trials]
        left.bar(x,rhos,color="#196d87")
        left.set_xticks(x,["Capped .025","Refined .025",".00625",".003125"],rotation=15)
        left.axhline(.1,color="#905046",lw=1,ls="--")
        left.set_ylabel("Actual / predicted decrease")
        left.set_title("Refined actual trials",loc="left")
    retained = data["lineages"]["retained"]
    baseline = retained["initial_max"]
    packet_index = {}
    for name,label,color in (("retained","Retained span","#38577b"),
                             ("scalars","Coupled scalars","#a54e39"),
                             ("momentum","Momentum comparison","#287f70")):
        if name not in data["lineages"]:
            continue
        arm = data["lineages"][name]
        values = [arm["initial_max"]] + [max(p["residuals"]) for p in arm["accepted"]]
        start = packet_index.get(arm["seed"]["sha256"],0)
        indices = np.arange(len(values))+start
        packet_index.update({p["packet"]["sha256"]:i for p,i in zip(arm["accepted"],indices[1:])})
        right.plot(indices,100*(1.-np.asarray(values)/baseline),"o-",label=label,color=color,ms=4)
    right.set_title("Accepted original progress",loc="left")
    right.set_xlabel("R11 accepted update along lineage")
    right.set_ylabel("Cumulative decrease from R11 base (%)")
    right.legend(frameon=False,fontsize=8)
    for axis in (left,right):
        axis.grid(axis="y",alpha=.2)
    figure.savefig(output,dpi=160)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest",type=Path)
    parser.add_argument("figure",type=Path)
    config = parser.parse_args()
    plot(json.loads(config.manifest.read_text()),config.figure)
