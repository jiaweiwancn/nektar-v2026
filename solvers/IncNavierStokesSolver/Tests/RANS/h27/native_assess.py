#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Cross-check the native adapter/action against independent planar algebra."""
import argparse,json,math,os,subprocess,time
from pathlib import Path
import numpy as np
from resources import THREADS,atomic,identity
import inlet_principal as planar
from reference_oracle import Column,Parameters

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--version',type=int,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--tag',default='v1')
    a=p.parse_args();r=a.root;science=r/'science';column=Column(Parameters())
    assert a.tag.isalnum(), 'simple immutable result tag required'
    ref=json.loads((science/'reference-kin-v1.json').read_text())['cases'][0]['roots'][0]['audit']
    inlet=json.loads((science/'fixed-inlet-v1.json').read_text())
    states=[]
    def add(label,row,reference=False):
        radius=float(row['y'])*column.p.H if reference else float(row['r'])
        states.append(dict(label=label,k=row['k'],e=row['e']if reference else row['epsilon'],
                           shear=row['g']if reference else row['G'],nu=column.p.nu,
                           nuA=column.auxiliary(radius/column.p.H),r=radius,H=column.p.H))
    for i,row in enumerate(ref['states']):
        if i<3 or row['y']in(.4,.5,.6,.75,1.):add(f'reference_{i}',row,True)
    for i,row in enumerate(inlet['reference_selected_cutoffs']):
        add(f'inlet_selected_{i}',row['min_witness'])
        add(f'inlet_selected_interior_{i}',row['positive_r_witness'])
    calibration=Path(planar.__file__).resolve().parents[5]/'docs/plans/wj-h5-linux/calibration.json'
    historical=next(v for v in json.loads(calibration.read_text())['candidates']if v['H_w_mm']==6.)
    old=planar.NeutralReference(historical['H_w'],historical['z0_e'])
    for i,radius in enumerate((.01*old.H,.02*old.H,.04*old.H)):
        add(f'inlet_positive_width_{i}',dict(r=radius,k=old.k_inner(radius),
            epsilon=column.A/max(radius,ref['c']),G=old.u_tau**2/old.D(radius)))
    scans=inlet['scans']
    for i in range(0,len(scans),max(1,len(scans)//6)):
        add(f'inlet_cutoff_{i}',scans[i]['positive_r_witness'])
    tails=inlet['cutoff_independent_negative_tail']
    if tails:
        add('inlet_tail_largest_r',max(tails,key=lambda x:x['r']))
        add('inlet_tail_worst',min(tails,key=lambda x:x['normalized_principal']))
    add('outer_inlet',inlet['outer_inlet']['witness'])
    control=inlet.get('infinite_inner_diffusion_reference_control',{})
    for i,row in enumerate(control.get('roots',[])):add(f'infinite_diffusion_{i}',row['witness'])
    for i,row in enumerate(inlet['ground_small_cutoff_control'][-2:]):add(f'small_cutoff_ground_{i}',row)
    inputs=science/f'native-states-{a.tag}.tsv'
    with inputs.open('x')as f:
        for s in states:f.write(s['label']+' '+' '.join(format(s[k],'.17g')for k in('k','e','shear','nu','nuA','r','H'))+'\n')
    runtime=r/f'runtime-v{a.version}';raw=science/f'native-principal-raw-{a.tag}.json';libraries=set()
    with raw.open('x')as output:
        child=subprocess.Popen([str(runtime/'WJH27PrincipalProbe'),str(inputs)],
            stdout=output,stderr=subprocess.PIPE,env=dict(os.environ,**THREADS))
        while child.poll()is None:
            try:
                for line in(Path('/proc')/str(child.pid)/'maps').read_text().splitlines():
                    pieces=line.split(None,5)
                    if len(pieces)==6 and pieces[-1].startswith('/'):
                        path=Path(pieces[-1]);
                        if path.is_file():libraries.add(path.resolve())
            except FileNotFoundError:pass
            time.sleep(.1)
        error=child.stderr.read().decode();assert child.returncode==0,error
    native=json.loads(raw.read_text());assert len(native['states'])==len(states)
    comparisons=[]
    for expected,got in zip(states,native['states']):
        assert expected['label']==got['label']
        k,e,g=expected['k'],expected['e'],expected['shear'];radius=expected['r'];y=radius/expected['H']
        gradient=np.zeros((3,3));gradient[0,2]=g
        rw=column.tensor(y,k,e,gradient)
        theta=got['theta'];n=np.array([math.cos(theta),math.sin(theta)])
        l=np.array([-math.sin(theta),math.cos(theta)]);h=np.outer(l,n)
        matrix=np.array([[0.,g],[0.,0.]])
        dr=planar.tensor(k,e,matrix,expected['nuA'],column.chi(y),h)
        action=expected['nu']*(h+h.T)-dr
        quality=float(np.sum(h*action))
        mu=expected['nu']+(1-column.chi(y))*expected['nuA']
        if column.chi(y):mu+=column.chi(y)*(-column.tensor(1,k,e,gradient)[0,2]/g)
        production=-float(np.sum(rw*gradient))
        errors=dict(tensor_over_k=float(np.max(abs(np.asarray(got['R'])-rw))/k),
            action_over_mu=abs(quality-got['planar_min'])/mu,
            split_mu_relative=abs(mu-got['mu'])/mu,
            production_over_E=abs(production-got['production'])/e)
        assert max(errors.values())<=1e-10,(got['label'],errors)
        fd_error=min(abs(v['difference_over_mu'])for v in got['finite_difference'][-2:])
        assert fd_error<=2e-8,(got['label'],fd_error)
        covariance=float(np.linalg.eigvalsh(np.asarray(got['R'])).min()/k)
        comparisons.append(dict(label=got['label'],native_principal_over_mu=got['planar_min_over_mu'],
            native_principal=got['planar_min'],theta=theta,r=radius,
            covariance_min_over_k=covariance,independent_errors=errors,
            finest_FD_error_over_mu=fd_error))
    selected=[x for x in comparisons if x['label'].startswith(('reference_0','inlet_selected'))]
    assert selected and all(x['native_principal_over_mu']<0 for x in selected)
    interior=[x for x in comparisons if x['label'].startswith('inlet_positive_width_')]
    assert any(x['native_principal_over_mu']<0 for x in interior)
    result=dict(schema='WJH27NativeIndependentComparisonV1',
        diagnostic_verification='PASS',model_gate='NOT_MET',scope='Constitutive necessary-condition rejection; no coupled solver/pilot',
        source_manifest=identity(r/f'source-v{a.version}.json'),
        assessment_script=identity(__file__),
        executable=identity(runtime/'WJH27PrincipalProbe.elf'),inputs=identity(inputs),
        raw=identity(raw),comparisons=comparisons,
        loaded_files=[identity(f)for f in sorted(libraries)],
        maximum_independent_error=max(max(x['independent_errors'].values())for x in comparisons))
    atomic(a.output,result);print(json.dumps(dict(diagnostic_verification='PASS',model_gate='NOT_MET',
        states=len(comparisons),maximum_independent_error=result['maximum_independent_error'])))

if __name__=='__main__':main()
