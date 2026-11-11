#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused regional-pressure case generation and default compatibility checks."""
import argparse
from copy import deepcopy
import importlib.util
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

from cls_epawf_m8_r1_cases import parameters,solver_info
from run_cls_accepted_awf import case_xml
from run_cls_epawf_picard import load,require


def canonical(xml):
    tree=ET.parse(xml).getroot()
    for row in tree.findall('./CONDITIONS/SOLVERINFO/I'):
        if row.get('PROPERTY') in ('RANSWallCheckpoint','RANSWallRestart'):row.set('VALUE','TRANSIENT_PATH')
    return ET.tostring(tree)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--baseline-script',type=Path);a=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='inner-pressure-cases-') as tmp:
        root=Path(tmp);common=dict(inner_refinement='cutoff-v1',companion='regularized',companion_reason='matched pressure-space qualification')
        base=case_xml('Accepted-AWF-05G',root/'p1',0,**common)
        quadratic=case_xml('Accepted-AWF-05G',root/'p2',0,pressure_space='inner-p2-v1',**common)
        left,right=[ET.parse(p).getroot() for p in (base,quadratic)];dl,dr=[load(p.parent/'definition.json') for p in (base,quadratic)]
        require(dr['tensor_field_size_estimates']+[1]==[784,784,784,400,187,54,1] and dr['tensor_unknown_estimate']==2994,'incorrect inner-P2 native-size prediction')
        require(dr['pressure_degree']==dict(inner=2,outer=1) and dr['case_identity'].endswith(':inner-p2-v1'),'undeclared pressure-space identity')
        require([(p.get('COMPOSITE'),p.get('NUMMODES')) for p in right.findall('./EXPANSIONS/E[@FIELDS="p"]')]==[('C[1]','3,3'),('C[2]','2,2')],'wrong regional pressure expansion')
        for tag in ('VERTEX','EDGE','ELEMENT'):
            require(ET.tostring(left.find('./GEOMETRY/'+tag))==ET.tostring(right.find('./GEOMETRY/'+tag)),'changed physical mesh: '+tag)
        require(parameters(left)==parameters(right),'pressure selection changed physical/temporal parameters')
        sl,sr=solver_info(left),solver_info(right)
        for k in ('RANSWALLCHECKPOINT','RANSWALLRESTART'):sl.pop(k,None);sr.pop(k,None)
        require(sl==sr,'pressure selection changed a physical solver policy')
        for name in ('radial_knots','tangential_knots','model_choices','H','z_K','z_m','h','Z','native_points','source_normal_order','source_tangential_order','work'):
            require(dl[name]==dr[name],'changed matched-case setting: '+name)
        require(right.find('./GEOMETRY/DOMAIN').text=='C[1,2]','regional expansions missing from native DOMAIN')
        for case,kw in [('Accepted-AWF-2G',common),('Accepted-AWF-05G',dict(common,inner_refinement='none'))]:
            try:case_xml(case,root/('invalid-'+case+kw['inner_refinement']),0,pressure_space='inner-p2-v1',**kw)
            except ValueError as error:require('matched cutoff-v1 H=.5' in str(error),'unexpected scope rejection')
            else:raise AssertionError('unsupported pressure control accepted')
        if a.baseline_script:
            spec=importlib.util.spec_from_file_location('archived_accepted_case_generator',a.baseline_script);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
            for case,kw in [('Accepted-AWF-05G',common),('Accepted-AWF-2G',dict(inner_refinement='cutoff-v1',outer_refinement='pressure-p2-v1'))]:
                prior=old.case_xml(case,root/('old-'+case),0,**kw)
                current=case_xml(case,root/('current-'+case),0,**kw)
                require(canonical(prior)==canonical(current),'default or existing outer-P2 XML changed')
    print('inner-P2 regional policy, unchanged physical mesh/parameters, pilot size, scope guards and default XML compatibility passed')


if __name__=='__main__':main()
