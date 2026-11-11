#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Continuous cubic input minima; read the exact sealed emitted expressions."""
import argparse,ast,json,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import write_json
import numpy as np
from numpy.polynomial import Polynomial as P
from rans_request_manifest import sha


def poly(n):
    if isinstance(n,ast.Constant):return P([float(n.value)])
    if isinstance(n,ast.Name) and n.id=='y':return P([0.,1.])
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -poly(n.operand)
    if isinstance(n,ast.BinOp):
        a,b=poly(n.left),poly(n.right)
        if isinstance(n.op,ast.Add):return a+b
        if isinstance(n.op,ast.Sub):return a-b
        if isinstance(n.op,ast.Mult):return a*b
        if isinstance(n.op,ast.Pow) and b.degree()==0:return a**int(b.coef[0])
    raise ValueError('unsupported polynomial input syntax')


def terms(n):
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Add):return terms(n.left)+terms(n.right)
    return [n]


def interval(n,lo,hi):
    if isinstance(n,ast.Constant):assert n.value==1;return lo,hi
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
        lo,hi=interval(n.left,lo,hi);return interval(n.right,lo,hi)
    assert isinstance(n,ast.Compare) and isinstance(n.left,ast.Name) and n.left.id=='y' and len(n.ops)==1
    value=poly(n.comparators[0]).coef[0]
    if isinstance(n.ops[0],(ast.Gt,ast.GtE)):lo=max(lo,value)
    elif isinstance(n.ops[0],(ast.Lt,ast.LtE)):hi=min(hi,value)
    else:raise ValueError('unsupported interval comparison')
    return lo,hi


def audit(expression,lower,upper):
    body=ast.parse(expression.replace('^','**'),mode='eval').body;rows=[]
    for term in terms(body):
        assert isinstance(term,ast.BinOp) and isinstance(term.op,ast.Mult)
        lo,hi=interval(term.left,lower,upper);p=poly(term.right);assert hi>lo
        candidates=[lo,hi]+[float(z.real) for z in p.deriv().roots() if abs(z.imag)<=64*np.finfo(float).eps*max(1.,abs(z.real)) and lo<z.real<hi]
        values=[float(p(x)) for x in candidates];i=int(np.argmin(values))
        error=64*np.finfo(float).eps*sum(abs(c)*max(1.,abs(lo),abs(hi))**j for j,c in enumerate(p.coef))
        rows.append(dict(lower=lo,upper=hi,minimum=values[i],location=candidates[i],arithmetic_bound=float(error),stationary_points=candidates[2:]))
    rows.sort(key=lambda x:x['lower']);assert rows[0]['lower']==lower and rows[-1]['upper']==upper
    assert all(a['upper']==b['lower'] for a,b in zip(rows[:-1],rows[1:]))
    return dict(minimum=min(x['minimum'] for x in rows),arithmetic_adjusted_lower_estimate=min(x['minimum']-x['arithmetic_bound'] for x in rows),pieces=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=dict(status='pass',cases={})
    if a.output.exists():raise FileExistsError('immutable audit destination already exists')
    for folder in sorted(a.pilot.glob('*-zero')):
        c=json.loads((folder/'construction.json').read_text());root=ET.parse(folder/'case.xml').getroot();rows={}
        for region,lo,hi,variables in [('Inner',0.,c['H'],('u','k')),('Outer',c['H'],c['Z'],('u','k','eps_tilde'))]:
            for variable in variables:
                expression=root.find("./CONDITIONS/FUNCTION[@NAME='InitialConditions"+region+"']/E[@VAR='"+variable+"']").get('VALUE')
                row=audit(expression,lo,hi);rows[region+'/'+variable]=row
                if variable!='u' and row['arithmetic_adjusted_lower_estimate']<=0:result['status']='fail'
        result['cases'][folder.name]=dict(xml_sha256=sha(folder/'case.xml'),fields=rows,
            scope='emitted cubic input and derivative stationary roots with conservative evaluation-roundoff estimate; physical endpoints included; floating numerical audit, not an interval-arithmetic theorem or evolved positivity claim')
    assert len(result['cases'])==4
    write_json(a.output,result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
