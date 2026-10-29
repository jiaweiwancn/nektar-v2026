#!/usr/bin/env python3
"""Corrected-AWF transport fixtures and an independent local reference oracle."""
import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from generate_m4 import awf_reference, generate_cls, info, parameter


def reference(k, height=2., point_p=1., roughness=.1, nu=1.5e-5,
              q=3.190675081624288, kappa=.440553905186659):
    original = awf_reference(k, height, roughness, nu)
    star = math.sqrt(k/q); numerator = star**3/kappa
    p0 = (k/q)**2/height*original["production_weight"]
    lam = height*p0/numerator
    cutoff = height*math.exp(1-lam)
    ze = 5.1*nu/math.sqrt(k)
    return dict(original, ustar=star, numerator=numerator, p0=p0, lam=lam,
                minimum_lambda=1+math.log(height/point_p), cutoff=cutoff,
                matching_velocity=k/q*original["j0"],
                log_z0=math.log(height)-kappa*(k/q*original["j0"])/star,
                epsilon_p=k**1.5/(2.55*max(point_p,ze)))


def generate_equilibrium(path, mode="Fixed", weight=1., on=.02, off=.1,
                         failure="RequirePreservation", point_p=1., point_n=4.,
                         audit=True, point_diagnostics=True, original=False, **options):
    generate_cls(path, kind="awf", **options)
    tree=ET.parse(path);root=tree.getroot()
    policy="AWFOriginal" if original else "AWFEquilibrium"
    info(root,"RANSWallModel",policy)
    for node in root.findall("CONDITIONS/BOUNDARYCONDITIONS/REGION/*"):
        if node.attrib.get("USERDEFINEDTYPE")=="AWFOriginal":node.attrib["USERDEFINEDTYPE"]=policy
    info(root,"RANSAWFWeight",mode);info(root,"RANSAWFFailurePolicy",failure)
    for name,value in {"RANSAWFFixedWeight":weight,"RANSAWFOn":on,"RANSAWFOff":off,
                       "RANSStripPointP":point_p,"RANSStripPointN":point_n,
                       "RANSStageAudit":int(audit),"RANSAWFPointDiagnostics":int(point_diagnostics),
                       "RANSAWFMinimumRe":100}.items():parameter(root,name,value)
    angle=options.get("span_angle",0.)
    parameter(root,"RANSAWFDirectionX",math.cos(angle));parameter(root,"RANSAWFDirectionY",math.sin(angle))
    ET.indent(root,space="  ");tree.write(path,encoding="utf-8",xml_declaration=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("path",type=Path)
    parser.add_argument("--mode",choices=("Fixed","Switched"),default="Switched")
    parser.add_argument("--weight",type=float,default=1.)
    parser.add_argument("--dt",type=float,default=.025);parser.add_argument("--steps",type=int,default=4)
    parser.add_argument("--modes",type=int,default=7);parser.add_argument("--points",type=int,default=12)
    args=vars(parser.parse_args());path=args.pop("path");generate_equilibrium(path,**args)
