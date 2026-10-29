#!/usr/bin/env python3
"""Render the recorded CLS convergence and original-AWF baseline evidence."""
import argparse
import json
import os
from pathlib import Path
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS"):
    os.environ[key]="1"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("validation",type=Path);p.add_argument("output",type=Path)
    args=p.parse_args();r=json.loads(args.validation.read_text());args.output.mkdir(parents=True,exist_ok=True)
    fig,axis=plt.subplots(figsize=(6,3.8),constrained_layout=True)
    steps=[.025,.0125,.00625]
    for field in ("u","v","k","eps_tilde","p"):
        errors=[r["cases"][f"cls_dt{dt:g}"]["field_audit"]["fields"][field]["l2_error"] for dt in steps]
        axis.loglog(steps,[e/errors[0] for e in errors],"o-",label=field)
    axis.loglog(steps,[(dt/steps[0])**2 for dt in steps],"k--",label="order 2")
    axis.set(xlabel="Timestep (s)",ylabel="L2 error / coarsest-step error",title="Full CLS transport: unsteady manufactured solution")
    axis.grid(True,which="both",alpha=.2);axis.legend(ncol=2,fontsize=8)
    for ext in ("png","pdf"):fig.savefig(args.output/("cls-temporal."+ext),dpi=180)
    plt.close(fig)
    rows=r["cases"]["awf"]["wall_stages"]
    fig,axes=plt.subplots(1,2,figsize=(9,3.6),constrained_layout=True)
    times=[v["time"] for v in rows]
    axes[0].plot(times,[v["mean_kP"] for v in rows],"o-",label="coupled original AWF")
    axes[0].axhline(3.190675081624288*.4**2,color="black",linestyle="--",label="CLS reference")
    axes[0].set(xlabel="Time (s)",ylabel="Strip representative k (m²/s²)",title="Uncorrected wall-cell evolution")
    axes[0].legend(fontsize=8)
    axes[1].plot(times,[v["mean_production"] for v in rows],"o-",label="original mean production")
    axes[1].plot(times,[v["mean_dissipation"] for v in rows],"s-",label="original mean dissipation")
    axes[1].set(xlabel="Time (s)",ylabel="Source density (m²/s³)",title="Original AWF retains its source mismatch")
    axes[1].legend(fontsize=8)
    for axis in axes:axis.grid(True,alpha=.2)
    for ext in ("png","pdf"):fig.savefig(args.output/("original-awf-baseline."+ext),dpi=180)


if __name__=="__main__":main()
