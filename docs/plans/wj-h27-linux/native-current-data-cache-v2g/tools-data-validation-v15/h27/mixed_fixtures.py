#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate tiny native mixed-component/Q3-scalar fixtures on Linux."""
import argparse,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'h26'))
from fixtures import mesh

def generate(out):
    out.mkdir(parents=True,exist_ok=False)
    for name in('closed','slip','physicaloutlet'):
        p=out/(name+'.xml');mesh(p,[0.,.5,1.],[0.,.15,.5,1.],quadrature=12)
        tree=ET.parse(p);root=tree.getroot();exp=root.find('EXPANSIONS')
        attrs=dict(exp[0].attrib);attrs['FIELDS']='k';ET.SubElement(exp,'E',attrs)
        cond=root.find('CONDITIONS');ET.SubElement(cond.find('VARIABLES'),'V',ID='3').text='k'
        ET.SubElement(cond.find("FUNCTION[@NAME='InitialConditions']"),'E',VAR='k',VALUE='0.008809015145515192')
        for reg in cond.find('BOUNDARYCONDITIONS'):
            boundary=int(reg.attrib['REF'])
            for bc in reg:
                var=bc.attrib['VAR']
                if var in('u','v'):
                    natural=(name!='closed'and boundary==3 and var=='u')or(name=='physicaloutlet'and boundary==1)
                    bc.tag='N'if natural else'D'
            ET.SubElement(reg,'D'if boundary==0 else'N',VAR='k',VALUE='0')
        ET.indent(root);tree.write(p,encoding='unicode');p.write_text(p.read_text()+'\n')
        with(out/(name+'.json')).open('x')as f:json.dump(dict(case=name,mode='mixed_foundation_only',H=.15,h=.10854503542535461),f,indent=2)
    tree=ET.parse(out/'physicaloutlet.xml');root=tree.getroot()
    for vertex in root.find('GEOMETRY/VERTEX'):
        x,z,d=map(float,vertex.text.split());vertex.text=f'{x:.17g} {z+.0875:.17g} {d:.17g}'
    ET.indent(root);tree.write(out/'scalar-shifted.xml',encoding='unicode')
    print(out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);generate(p.parse_args().out)
