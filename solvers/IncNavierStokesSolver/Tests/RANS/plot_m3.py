#!/usr/bin/env python3
"""Plot the recorded standard-RANS temporal and ABL acceptance results."""
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
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("validation",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args();report=json.loads(args.validation.read_text())
    args.output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(9,3.7),constrained_layout=True)
    timesteps=[.05,.025,.0125]
    for axis,order in zip(axes,(1,2)):
        for field in ("u","v","k","eps","p"):
            errors=[report["cases"][f"mms_o{order}_dt{dt:g}"]["field_audit"]["fields"][field]["l2_error"] for dt in timesteps]
            axis.loglog(timesteps,[e/errors[0] for e in errors],"o-",label=field)
        axis.loglog(timesteps,[(dt/timesteps[0])**order for dt in timesteps],"k--",label=f"order {order}")
        axis.set(title=f"IMEX{order}: manufactured URANS",xlabel="Timestep (s)",ylabel="L2 error / coarsest-step error")
        axis.grid(True,which="both",alpha=.2);axis.legend(ncol=2,fontsize=8)
    fig.savefig(args.output/"temporal-convergence.png",dpi=180)
    fig.savefig(args.output/"temporal-convergence.pdf")
    plt.close(fig)
    degrees=[4,6,8]
    fig,axis=plt.subplots(figsize=(5.5,3.7),constrained_layout=True)
    for field in ("u","k","eps"):
        errors=[report["cases"][f"abl_p{degree}"]["field_audit"]["fields"][field]["l2_error"] for degree in degrees]
        axis.semilogy(degrees,[e/errors[0] for e in errors],"o-",label=field)
    axis.set(xlabel="Velocity/scalar polynomial degree",ylabel="L2 error / degree-4 error",title="Standard neutral ABL: spatial refinement, t=1 s",xticks=degrees)
    axis.grid(True,which="both",alpha=.2);axis.legend()
    fig.savefig(args.output/"abl-refinement.png",dpi=180)
    fig.savefig(args.output/"abl-refinement.pdf")


if __name__=="__main__":main()
