#!/usr/bin/env python3
"""Linux CV2 verification. Short functionality checks, not experimental fits."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from cao_tamura_experiment import (ROOT, DATA, MODELS, THREADS, allocation, config,
    environment, execute, frozen_data, generate, inlet_values, observe, prepare,
    probes, records, reference, sha, statistics, compare as score, write_csv, write_json)
from cao_tamura_experiment import stage_metrics
from generate_terrain import plane, hill
from generate_m3 import element
from generate_m4 import parameter, info
from generate_m6 import original_j0
from validate_terrain import compare, corner_compatibility
from terrain_initial import prepare as terrain_prepare
from terrain_continuity import audit as continuity


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Validation:
    def __init__(self, build, output, cpus):
        self.build = environment(build); self.output = output.resolve(); self.cpus = cpus
        allocation(4, cpus)
        require(not self.output.is_relative_to(ROOT), 'Store numerical output outside the source checkout')
        self.output.mkdir(parents=True, exist_ok=False)
        self.results = {}
        changed=subprocess.check_output(['git','-C',ROOT,'ls-files','--modified','--others','--exclude-standard'],text=True).splitlines()
        self.source_hashes={p:sha(ROOT/p) for p in changed if (ROOT/p).is_file() and Path(p).suffix in ('.cpp','.h','.py','.txt')}

    def save(self, name, result):
        self.results[name] = result
        write_json(self.output/'progress.json', self.results)
        print(name, 'PASS', flush=True)

    def native(self, binary, args, work, label, ranks=1, timeout=900):
        return execute(allocation(ranks, self.cpus)+[self.build/binary]+list(args), work, label, timeout)

    def field(self, session, label='field'):
        result = session.parent/(label+'.csv')
        self.native('RANSWallFieldAudit', [session, session.with_suffix('.fld'), result], session.parent, label)
        return result

    def solve(self, name, model='standard_fixed', diffusion='HighRe', smooth=False,
              steps=2, modify=None, curved=False, degree=5, layers=8, ranks=1, prepared=None):
        directory = self.output/name; directory.mkdir(); path = directory/'neutral.xml'
        if prepared:
            prepared(path)
        elif curved:
            hill(path, model=model, degree=degree, layers=layers, steps=steps, dt=2e-5)
        else:
            plane(path, model=model, degree=degree, layers=layers, nx=2,
                  steps=steps, dt=.001, diffusion=diffusion)
        tree = ET.parse(path); root = tree.getroot()
        parameter(root, 'RANSDebug', 1); parameter(root, 'RANSStageAudit', 1)
        if smooth:
            q = 10/3 if model.startswith('standard') else 3.190675081624288
            parameter(root, 'RANSRoughnessHeight', 0.)
            parameter(root, 'MatchingVelocity', .09*original_j0(q*.09, .008, 0., 1.5e-5))
        if modify:
            modify(root)
        tree.write(path, encoding='utf-8', xml_declaration=True)
        preparation = terrain_prepare(path, self.build, int(self.cpus.split(',')[0])) if curved else None
        self.native('IncNavierStokesSolver', [path.name, '--use-scotch'], directory, 'solver', ranks, 1800)
        field = self.field(path)
        stages = records(directory/'neutral.rans.csv')
        div = max(float(r['weak_divergence_l2']) for r in stages)
        require(div <= 1e-7 and min(float(r['min_k']) for r in stages) > 0 and min(float(r['min_eps']) for r in stages) > 0,
                name+' stage/positivity gate')
        wall = records(directory/'neutral.terrain-wall.csv')
        projected = max(float(r['normal_projected_max']) for r in wall)
        require(projected <= 1e-7, name+' projected-normal gate')
        free = max(float(r['free_coefficient_l2']) for r in records(directory/'neutral.rans-residual.csv') if float(r['phase']) == 1)
        require(free <= 1e-7, name+' free-equation residual gate')
        result = dict(session=str(path), field=str(field), ranks=ranks, case_sha256=sha(path),
                      maximum_weak_divergence=div, maximum_projected_normal=projected,
                      maximum_free_residual=free, preparation=preparation, continuity=continuity(directory))
        self.save(name, result)
        return path, field

    def probes_suite(self):
        from scipy.integrate import quad
        # Polynomial physical fields on affine meshes are represented exactly.
        # Off-grid and shared-face queries make nearest-GLL extraction fail.
        for dim, angle in ((2, 0.), (2, .35), (3, 0.), (3, .35)):
            directory = self.output/f'polynomial-d{dim}-a{angle}'; directory.mkdir()
            path = directory/'sample.xml'
            plane(path, model='standard_original', degree=4, nx=3, layers=3, dimension=dim, angle=angle)
            tree = ET.parse(path); root = tree.getroot(); coord = 'y' if dim == 2 else 'z'
            group = root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
            formulas = dict(u=f'1+.2*x+.3*{coord}+.1*x*{coord}+.02*x^2',
                            v='.07+.03*y' if dim == 3 else f'.1+.02*x-.04*{coord}',
                            w=f'.1+.02*x-.04*{coord}', k='.3', eps='.2', p='.2')
            for e in group:
                e.attrib['VALUE'] = formulas[e.attrib['VAR']]
            tensor = element(root.find('CONDITIONS'), 'FUNCTION', NAME='RANSProbeTestTensor')
            for i, name in enumerate(('xx','yy','zz','xy','xz','yz')):
                element(tensor, 'E', VAR='rans_R_'+name, VALUE=f'{i+1}+.2*x+.3*{coord}'+('+.4*y' if dim == 3 else ''))
            tree.write(path, encoding='utf-8', xml_declaration=True)
            c,s = math.cos(angle), math.sin(angle)
            queries = []
            for i,(a,h) in enumerate(((.037,.012), (.4/3,.02), (.4/3,.008), (.27,.04), (.071,.004), (.37,.0079))):
                queries.append(dict(id='p'+str(i), x=c*a-s*h, span=.011 if dim==3 else 0., z=s*a+c*h))
            queries.append(dict(id='outside',x=2.,span=0.,z=.03))
            points = directory/'points.csv'; write_csv(points, queries)
            outputs = []
            for ranks in (1,2,4):
                work = directory/f'r{ranks}'; work.mkdir(); shutil.copy2(path, work/'sample.xml')
                out = work/'points.csv'
                self.native('RANSObservationAudit', [work/'sample.xml','--manufactured',points,out], work, 'native', ranks)
                outputs.append(records(out))
            error = inverse = auxiliary_error = 0.
            for r, point in zip(outputs[0], queries):
                if point['id'] == 'outside':
                    require(r['kind']=='outside_resolved', 'Outside queries must not extrapolate'); continue
                x,y,z = point['x'],point['span'],point['z']
                inverse = max(inverse, float(r['inverse_residual']))
                if r['kind'] == 'resolved_native':
                    expected = dict(u=1+.2*x+.3*z+.1*x*z+.02*x*x,
                                    v=.07+.03*y if dim==3 else 0., w=.1+.02*x-.04*z, k=.3, epsilon_model=.2)
                    expected.update({'R_'+name:i+1+.2*x+.3*z+(.4*y if dim==3 else 0.)
                                     for i,name in enumerate(('xx','yy','zz','xy','xz','yz'))})
                    error = max(error, max(abs(float(r[k])-v) for k,v in expected.items()))
                else:
                    require(r['kind']=='auxiliary_velocity' and r['R_xx']==r['k']==r['epsilon_model']=='', 'Auxiliary support declaration')
                    # Independent quadrature of the original AWF momentum ODE.
                    eta = float(r['eta']); H=.008; K=float(r['aux_KP']); nu=1.5e-5; ks=.002
                    rh=ks*math.sqrt(K)/(70*nu)
                    intercept=10.7*(1-rh**max(.5-.4*rh**.7,1-.79*rh**(-.28)))
                    zv=nu*intercept/math.sqrt(K); beta=.09*2.55*math.sqrt(K)
                    D=lambda h:nu+beta*max(0.,h-zv)
                    def integral(h, moment):
                        cuts=sorted(set([0.,h]+[v for v in (ks,zv) if 0<v<h]))
                        return sum(quad(lambda e:(1 if moment==0 else max(e-ks,0.))/D(e),a,b,epsabs=1e-12,epsrel=2e-12)[0] for a,b in zip(cuts,cuts[1:]))
                    foot=float(r['foot_x']); hx=foot-s*H; hz=math.tan(angle)*foot+c*H
                    uh=1+.2*hx+.3*hz+.1*hx*hz+.02*hx*hx
                    wh=.1+.02*hx-.04*hz; th=c*uh+s*wh
                    force=float(r['aux_Ft']); tau=(th-force*integral(H,1))/integral(H,0)
                    tangent=tau*integral(eta,0)+force*integral(eta,1)
                    normal=-eta*(float(r['aux_dsUP'])+float(r['aux_dspanVP']))
                    expected_u=c*tangent-s*normal; expected_w=s*tangent+c*normal
                    auxiliary_error=max(auxiliary_error,abs(float(r['u'])-expected_u),abs(float(r['w'])-expected_w))
                    if dim==3:
                        force=float(r['aux_Fspan']); vh=.07+.03*y
                        tau=(vh-force*integral(H,1))/integral(H,0)
                        span=tau*integral(eta,0)+force*integral(eta,1)
                        auxiliary_error=max(auxiliary_error,abs(float(r['v'])-span),
                            abs(float(r['u'])-(c*tangent-s*normal)),abs(float(r['w'])-(s*tangent+c*normal)))
            require(error < 5e-10 and inverse < 1e-10 and auxiliary_error < 5e-10, 'Native polynomial/frame/quadrature gate '+str((error,inverse,auxiliary_error)))
            rank_error=0.
            for data in outputs[1:]:
                for a,b in zip(outputs[0],data):
                    require(a['id']==b['id'] and a['kind']==b['kind'] and a['element']==b['element'], 'MPI geometry/support mismatch')
                    for key in a:
                        if key not in ('id','kind','origin','owner','element') and a[key] and b[key]:
                            rank_error=max(rank_error,abs(float(a[key])-float(b[key]))/max(1.,abs(float(a[key]))))
            require(rank_error <= 1e-6, 'Native probe rank agreement')
            self.save(directory.name, dict(polynomial_max_error=error, inverse_max=inverse,
                auxiliary_quadrature_max_error=auxiliary_error, rank_scaled_max_error=rank_error, ranks=[1,2,4], queries=queries))
        # Curved native inverse and C1 support use constant physical fields.
        directory=self.output/'curved-probes';directory.mkdir();path=directory/'sample.xml'
        hill(path,model='standard_original',degree=4,layers=4)
        tree=ET.parse(path);root=tree.getroot()
        for e in root.findall("CONDITIONS/FUNCTION[@NAME='InitialConditions']/E"):
            e.attrib['VALUE']={'u':'1','v':'.01','k':'.3','eps':'.2','p':'.2'}[e.attrib['VAR']]
        tree.write(path,encoding='utf-8',xml_declaration=True)
        queries=[]
        for foot,eta in ((-.045,.013),(.025,.003),(.1,.004),(-.1,.004),(0.,.008)):
            h=.02*(1+math.cos(math.pi*foot/.1)) if abs(foot)<.1 else 0.
            hp=-.04*math.pi/.2*math.sin(math.pi*foot/.1) if abs(foot)<.1 else 0.
            n=math.hypot(1,hp)
            queries.append(dict(id=str(len(queries)),x=foot-eta*hp/n,span=0,z=h+eta/n))
        points=directory/'points.csv';write_csv(points,queries);outputs=[]
        for ranks in (1,4):
            work=directory/f'r{ranks}';work.mkdir();shutil.copy2(path,work/'sample.xml');out=work/'points.csv'
            self.native('RANSObservationAudit',[work/'sample.xml','--manufactured',points,out],work,'native',ranks)
            outputs.append(records(out))
        require(outputs[0][0]['kind']=='resolved_native' and abs(float(outputs[0][0]['u'])-1)<5e-10,'Curved constant interpolation')
        require(all(outputs[0][i]['kind']=='unsupported_C1_foot' for i in (2,3)),'C1 feet must be explicitly unsupported below hm')
        error=max(abs(float(a[k])-float(b[k]))/max(1.,abs(float(a[k]))) for a,b in zip(*outputs)
                  for k in ('u','v','w','eta','aux_dsUP') if a[k] and b[k])
        require(error<1e-6,'Curved observation MPI agreement')
        self.save('curved-probes',dict(rank_scaled_error=error, kinds=[r['kind'] for r in outputs[0]]))
        for label,model,forcing,dealias,kind in (
            ('log','standard',False,False,'unsupported_log_layer'),
            ('body','standard_original',True,False,'unsupported_auxiliary_forcing'),
            ('dealias','standard_original',False,True,'unsupported_auxiliary_forcing')):
            directory=self.output/('unsupported-'+label);directory.mkdir();path=directory/'sample.xml'
            plane(path,model=model,degree=3,nx=1,layers=3)
            tree=ET.parse(path);root=tree.getroot()
            if forcing:
                group=element(root,'FORCING');body=element(group,'FORCE',TYPE='Body')
                element(body,'BODYFORCE','ObservationTestBody')
                group=element(root.find('CONDITIONS'),'FUNCTION',NAME='ObservationTestBody')
                for name in ('u','v'):element(group,'E',VAR=name,VALUE='0')
            if dealias:info(root,'SpectralHPDealiasing','True')
            tree.write(path,encoding='utf-8',xml_declaration=True)
            points=directory/'points.csv';write_csv(points,[dict(id='below',x=.13,span=0,z=.004),dict(id='above',x=.13,span=0,z=.02)])
            out=directory/'output.csv'
            self.native('RANSObservationAudit',[path,'--manufactured',points,out],directory,'native')
            data=records(out)
            require(data[0]['kind']==kind and data[0]['u']==data[0]['R_xx']=='' and data[1]['kind']=='resolved_native','Unsupported auxiliary declaration')
            self.save('unsupported-'+label,dict(kinds=[r['kind'] for r in data]))

    def controls_suite(self):
        for model in ('standard_fixed','cls_fixed'):
            high,field=self.solve(model+'-smooth',model=model,smooth=True,steps=3,degree=6,layers=10)
            root=ET.parse(high).getroot();p={e.text.split('=',1)[0].strip():float(e.text.split('=',1)[1]) for e in root.findall('CONDITIONS/PARAMETERS/P')}
            errors={k:0. for k in ('u','v','k','epsilon')}
            for row in records(field):
                eta=float(row['y']);expected=dict(u=p['MatchingVelocity']+p['Ustar']/p['Kappa']*math.log(eta/p['H']),
                    v=0.,k=p['Qstar']*p['Ustar']**2,epsilon=p['Ustar']**3/(p['Kappa']*eta))
                for k,v in expected.items():
                    actual=float(row['eps_tilde' if model.startswith('cls') else 'eps']) if k=='epsilon' else float(row[k])
                    errors[k]=max(errors[k],abs(actual-v)/(p['Ustar'] if k=='v' else abs(v)))
            require(max(errors.values())<1e-3, 'Smooth-anchor neutral preservation')
            wall=records(high.parent/'neutral.terrain-wall.csv')
            require(all(float(r['chi_min'])==float(r['chi_max'])==1 for r in wall),'Smooth correction must remain active')
            self.save(model+'-smooth-errors',dict(normalized_linf=errors, anchor=p['MatchingVelocity'], ks=0., steps=3))
            self.solve(model+'-finite-nu',model=model,diffusion='Molecular',smooth=True,steps=2)
            original,origfield=self.solve(model+'-original',model=model.replace('_fixed','_original'),steps=2)
            zero,zerofield=self.solve(model+'-zero',model=model,steps=2,
                modify=lambda r:parameter(r,'RANSAWFFixedWeight',0.))
            errors=compare(origfield,zerofield);require(max(errors.values())<=1e-6,'Original recovery')
            self.save(model+'-zero-recovery',errors)
        # Both closures keep integrated continuity and expose the C1 max defect.
        for model in ('standard_fixed','cls_fixed'):
            path,field=self.solve(model+'-curved',model=model,steps=1,curved=True,degree=4,layers=5,ranks=4)
            record=self.results[model+'-curved'];require(record['continuity']['maximum_transpiration']>1e-6,'Curved transpiration must be exercised')
            self.save(model+'-C1-compatibility',corner_compatibility(path.parent))

    def generator_suite(self):
        frozen_data(); cases={}; meshes={}
        for surface in ('S','R'):
            for model in MODELS:
                name=surface+'-'+model; output=self.output/name
                result=generate(dict(surface=surface,model=model,ks=0. if surface=='S' else .006,
                    geometry='floor',degree=4,layers=7,steps=2),output)
                initial=prepare(output/'case.xml',self.build)
                cases[name]=dict(reference=result['reference'], preparation=initial)
                tree=ET.parse(output/'case.xml');meshes[name]=ET.tostring(tree.getroot().find('GEOMETRY'))
        require(len(set(meshes.values()))==1,'Surface/model comparisons must use the same mesh')
        for width in (1.,2.,4.,8.):
            result=generate(dict(geometry='hill',ramp_width_h=width,degree=4,layers=7),self.output/('width'+str(width)))
            path=self.output/('width'+str(width))/'case.xml'
            if width==1.:mesh=ET.tostring(ET.parse(path).getroot().find('GEOMETRY'))
            require(ET.tostring(ET.parse(path).getroot().find('GEOMETRY'))==mesh,'Width changed mesh')
        # Full native hex preparation, including more than one spanwise cell.
        result=generate(dict(surface='R',model='cls_original',ks=.006,geometry='hill',
             dimension=3,span_layers=2,degree=3,layers=5),self.output/'R-3d-hill')
        cases['R-3d-hill']=prepare(self.output/'R-3d-hill/case.xml',self.build)
        negatives=[dict(hm=.02),dict(point_p=.02),dict(point_n=.005),dict(ks=.02),
                   dict(ks=.001),dict(ramp_end=-.05),dict(ramp_width_h=30),dict(model='cls_log'),
                   dict(diffusion='unknown'),dict(outer_intensity=0),dict(degree=2)]
        for value in negatives:
            try:config(value)
            except ValueError:pass
            else:raise RuntimeError('Invalid configuration accepted: '+str(value))
        for surface in ('S','R'):
            for scope in ('flat','hill'):
                probes(surface,scope,self.output/(surface+'-'+scope+'-points.csv'))
        self.save('generator',dict(cases=cases, identical_comparison_mesh=True, identical_width_mesh=True,
                                  rejected_configs=negatives, cv1_snapshot_sha256=sha(DATA/'snapshot.json')))

    def workflow_suite(self):
        # Actual finite-depth inlet short solves, observation/statistics/objective,
        # plus experimental-input history rejection. Not a developed calibration.
        for surface,model in (('S','standard_fixed'),('R','cls_fixed')):
            directory=self.output/(surface+'-workflow')
            generate(dict(surface=surface,model=model,ks=0. if surface=='S' else .006,
                          degree=3,layers=6,steps=2),directory)
            prepare(directory/'case.xml',self.build)
            self.native('IncNavierStokesSolver',['case.xml','--use-scotch'],directory,'solver',2,1800)
            stages=stage_metrics(directory/'case.xml')
            require(stages['completed_stages']==2 and stages['minimum_k']>0 and stages['minimum_epsilon']>0
                    and stages['maximum_guard_activations']==0 and max(stages[k] for k in
                    ('maximum_weak_divergence','maximum_projected_normal','maximum_normal_integral_defect','maximum_resolved_mass_defect'))<=1e-7,
                    'Finite-depth stage/normal/mass/positivity gate')
            points=directory/'probes.csv';probes(surface,'flat',points)
            outputs=[]
            for i in (1,2):
                field=directory/f'case_{i}.chk';out=directory/f'points-{i}.csv'
                observe(directory/'case.xml',field,points,out,self.build,1,self.cpus)
                outputs.append(out)
            initial_observation=directory/'points-initial.csv'
            observe(directory/'case.xml',directory/'initial.fld',points,initial_observation,self.build,1,self.cpus)
            try:statistics([initial_observation,outputs[0]],directory/'invalid-initial-statistics.csv')
            except ValueError:pass
            else:raise RuntimeError('Initialization-only field entered statistics')
            statistics(outputs,directory/'statistics.csv')
            objective=score(directory/'statistics.csv',directory/'manifest.json',directory/'candidate.json',True)
            require(len(objective['results'])==3,'Flat objective must use only three calibration series')
            require(any(not r['complete'] for r in objective['results']), 'Below-hm stress support must remain incomplete')
            # Native history: one step plus one step across ranks vs same two steps.
            resumed=directory/'resumed';resumed.mkdir()
            helper=ROOT/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
            execute([sys.executable,helper,directory/'case.xml',directory/'case_1.chk',resumed/'resume.xml','--steps','1'],resumed,'prepare-restart')
            self.native('IncNavierStokesSolver',['resume.xml','--use-scotch'],resumed,'restart',1,1800)
            first=directory/'whole.csv';second=resumed/'resumed.csv'
            self.native('RANSWallFieldAudit',[directory/'case.xml',directory/'case.fld',first],directory,'whole-audit')
            self.native('RANSWallFieldAudit',[resumed/'resume.xml',resumed/'resume.fld',second],resumed,'resumed-audit')
            errors=compare(first,second);require(max(errors.values())<=1e-6,'Experimental history rank continuation')
            resumed_observation=resumed/'points.csv'
            observe(resumed/'resume.xml',resumed/'resume.fld',points,resumed_observation,self.build,1,self.cpus)
            statistics([outputs[0],resumed_observation],resumed/'across-history.csv')
            rejections={}
            def changed_expression(root):
                for pattern in ("CONDITIONS/FUNCTION[@NAME='RANSExperimentalInlet']/E[@VAR='u']",
                                "CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='3']/D[@VAR='u']"):
                    node=root.find(pattern);node.set('VALUE','('+node.attrib['VALUE']+')+0.01*Uinf')
            for label,mutation in (
                ('displacement',lambda r:parameter(r,'RANSExperimentDisplacement',.004)),
                ('profile_parameter',lambda r:parameter(r,'Uinf',1.4)),
                ('profile_expression',changed_expression),
                ('wall_property',lambda r:parameter(r,'RANSRoughnessHeight',.001))):
                work=resumed/label;work.mkdir();tree=ET.parse(resumed/'resume.xml');mutation(tree.getroot())
                tree.write(work/'reject.xml',encoding='utf-8',xml_declaration=True)
                try:self.native('IncNavierStokesSolver',['reject.xml','--use-scotch'],work,'reject',1,120)
                except RuntimeError:
                    text=(work/'reject.log').read_text()
                    require('RANS restart model/configuration mismatch' in text or 'inlet boundary differs' in text,
                            'Wrong history rejection reason')
                    rejections[label]=True
                else:raise RuntimeError('Changed experimental history accepted')
            self.save(surface+'-workflow',dict(stages=stages,history_scaled_errors=errors,rejected_history_changes=rejections,
                       candidate_objective=objective, statistics=json.loads((directory/'statistics.csv.json').read_text())))

    def utilities_suite(self):
        import cao_tamura_experiment as experiment
        from cao_tamura_observables import speedup, reattachment, recovery_zeros
        # Deliberately synthetic native-shaped rows verify only the statistical
        # algebra. Verification output is ineligible for experimental scoring.
        directory=self.output/'statistics-oracle';directory.mkdir();inputs=[]
        for i,(t,u) in enumerate(((0.,1.),(1.,3.),(3.,7.))):
            row=dict(id='oracle',time=t,x=0,span=0,z=.02,origin='checkpoint',kind='resolved_native',
                     u=u,v=0,w=2*u,k=.3,epsilon_model=.2,R_xx=2,R_yy=-.01,R_zz=1,R_xy=0,R_xz=-.5,R_yz=0)
            auxiliary=dict(row,id='auxiliary',kind='auxiliary_velocity',k='',epsilon_model='')
            for key in ('R_xx','R_yy','R_zz','R_xy','R_xz','R_yz'):auxiliary[key]=''
            path=directory/f'input-{i}.csv';write_csv(path,[row,auxiliary]);inputs.append(path)
        receipt=statistics(inputs,directory/'statistics.csv',verification=True)
        row,aux=records(directory/'statistics.csv')
        expected=dict(u=4.,w=8.,R_xx=7.,R_yy=-.01,R_zz=21.,R_xz=9.5)
        error=max(abs(float(row[k])-v) for k,v in expected.items())
        require(error<5e-10 and receipt['negative_diagonal_count']==1 and aux['R_xx']=='', 'Time/covariance/missing-data oracle')
        # Real statistics require unmodified native observer receipts.
        try:statistics(inputs,directory/'invalid.csv')
        except FileNotFoundError:pass
        else:raise RuntimeError('Unbound observations entered statistics')
        frozen_copy=self.output/'frozen-copy';shutil.copytree(DATA,frozen_copy)
        old=experiment.DATA
        try:
            experiment.DATA=frozen_copy
            (frozen_copy/'profiles.csv').write_bytes((frozen_copy/'profiles.csv').read_bytes()+b'\n')
            try:experiment.frozen_data()
            except ValueError:pass
            else:raise RuntimeError('Changed CV1 snapshot accepted')
        finally:experiment.DATA=old
        require(abs(speedup(2.75,2.,3.)['value']-.375)<5e-10 and not speedup(1.,0.,3.)['available'],
                'Speedup and near-zero denominator oracle')
        dividing=[]
        for j,z in enumerate((.002,.004,.006)):
            for i in range(30):
                x=.11+.01*i
                dividing.append(dict(id=f'dividing:{j}:{i}',x=x,span=0,z=z,u=x-(.237-3*z)))
        xr=reattachment(dividing)
        require(xr['available'] and abs(xr['x_R_m']-.237)<5e-10 and xr['tail_method_difference_over_H']<5e-10,
                'Two/three-tail ground extrapolation oracle')
        altered=[dict(r,u=float(r['u'])-.2*float(r['z'])**2) for r in dividing]
        require(reattachment(altered)['tail_method_difference_over_H']>1e-6,'Tail method sensitivity must be visible')
        require(recovery_zeros([(0.,-1.),(1.,None),(2.,1.)])==[], 'Recovery must not bridge missing probes')
        dividing[0]['u']='';require(not reattachment(dividing)['available'],'Incomplete dividing support')
        for surface in ('S','R'):
            for scope in ('hill','reference','local-reference','dividing'):
                probes(surface,scope,self.output/(surface+'-'+scope+'.csv'))
        crest=lambda name:next(r for r in records(self.output/name) if r['id']=='crest_speedup')
        require(float(crest('S-hill.csv')['z'])==.05 and float(crest('S-reference.csv')['z'])==.01,
                'Vertical-height crest reference mapping')
        # A provisional manifest cannot expose holdout scores before CV3 freeze.
        candidate=self.output/'candidate';generate({},candidate)
        try:score(directory/'statistics.csv',candidate/'manifest.json',candidate/'invalid-score.json')
        except ValueError as e:require('Freeze a versioned CV3' in str(e),'Wrong pre-freeze score rejection')
        else:raise RuntimeError('Holdout scores exposed before calibration freeze')
        # End-to-end comparison uses a miniature synthetic dataset and bound
        # receipts. No frozen experimental hill/holdout scores are evaluated.
        fixture=self.output/'comparison-fixture';fixture.mkdir();data=fixture/'data';data.mkdir()
        spec=[('speed','speedup',.375),('rms_ratio','sigma_u_over_flat_same_height',2.),
              ('shear_ratio','uv_over_flat_same_height',2.),('rms','sigma_u_over_Uinf',1.)]
        write_json(data/'series.json',[dict(id=n,case='SS',role='validation',quantity=q) for n,q,_ in spec])
        write_csv(data/'profiles.csv',[dict(series_id=n,point_id='0',value=v,weight=1) for n,_,v in spec])
        write_csv(data/'reported_scalars.csv',[dict(id='crest_S_SS',value=.375),dict(id='reattachment_SS',value=.237/.04)])
        write_json(data/'snapshot.json',dict(version='synthetic-operator-fixture',files={p.name:sha(p) for p in data.iterdir()}))
        hill_cfg=config(dict(geometry='hill'));floor_cfg=config({});ref=dict(Uinf=3.,ustar=.1)
        common=dict(cv1_snapshot_sha256=sha(data/'snapshot.json'),reference=ref,
                    calibration_complete=True,calibration_version='synthetic-only')
        hm={**common,'config':hill_cfg,'physics_fingerprint':'synthetic-hill'}
        fm={**common,'config':floor_cfg,'physics_fingerprint':'synthetic-floor'}
        write_json(fixture/'hill.json',hm);write_json(fixture/'floor.json',fm)
        hrows=[];frows=[]
        for name,_,_ in spec+[('crest_speedup','speedup',.375)]:
            base=dict(id=name if name=='crest_speedup' else name+':0',x=0.,span=0.,kind='resolved_native',
                      origin='checkpoint_statistics',time_begin=1.,time_end=3.)
            height=.01 if name=='crest_speedup' else .02
            hrows.append(dict(base,z=.04+height,u=2.75,R_xx=9. if name=='rms' else 4.,R_xz=-.5))
            frows.append(dict(base,z=height,u=2.,R_xx=1.,R_xz=-.25))
        for name,rows,fingerprint in (('hill',hrows,'synthetic-hill'),('floor',frows,'synthetic-floor'),
                 ('dividing',[dict(r,u=float(r['x'])-(.237-3*float(r['z']))) for r in dividing],'synthetic-hill')):
            path=fixture/(name+'.csv');write_csv(path,rows)
            write_json(str(path)+'.json',dict(output_sha256=sha(path),experimental_fingerprint=fingerprint))
        old=experiment.DATA
        try:
            experiment.DATA=data
            result=score(fixture/'hill.csv',fixture/'hill.json',fixture/'comparison.json',False,
                         fixture/'floor.csv',fixture/'floor.json',fixture/'dividing.csv')
            require(all(r['complete'] and r['complete_range_rmse']<5e-10 for r in result['results'])
                    and result['scalars']['crest_speedup']['absolute_error']<5e-10
                    and result['scalars']['reattachment']['absolute_error_over_H']<5e-10,'End-to-end observable comparison')
            frows[0]['z']+=.001;write_csv(fixture/'floor.csv',frows)
            write_json(fixture/'floor.csv.json',dict(output_sha256=sha(fixture/'floor.csv'),experimental_fingerprint='synthetic-floor'))
            try:score(fixture/'hill.csv',fixture/'hill.json',fixture/'invalid.csv',False,fixture/'floor.csv',fixture/'floor.json')
            except ValueError as e:require('same physical vertical yprime' in str(e),'Wrong reference-coordinate rejection')
            else:raise RuntimeError('Wrong vertical reference accepted')
        finally:experiment.DATA=old
        self.save('utilities',dict(time_covariance_max_error=error,negative_diagonals_retained=True,
                    missing_stress_retained=True,unbound_statistics_rejected=True,corrupt_snapshot_rejected=True,
                    speedup_oracle=.375,reattachment_oracle=xr,missing_dividing_support_rejected=True,
                    holdout_before_freeze_rejected=True,comparison_fixture='synthetic mini-dataset only',
                    comparison_profile_errors=[r['complete_range_rmse'] for r in result['results']],
                    wrong_vertical_reference_rejected=True))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cpus',required=True)
    p.add_argument('--suite',choices=('probes','controls','generator','workflow','utilities'),required=True)
    a=p.parse_args();require(sys.platform.startswith('linux'),'Linux is authoritative')
    os.environ.update(THREADS);os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    v=Validation(a.build,a.output,a.cpus)
    getattr(v,a.suite+'_suite')()
    receipt=dict(status='PASS',suite=a.suite,source_parent=subprocess.check_output(['git','-C',ROOT,'rev-parse','HEAD'],text=True).strip(),
                 source_sha256=v.source_hashes,
                 binaries={n:sha(v.build/n) for n in ('IncNavierStokesSolver','RANSObservationAudit','RANSWallFieldAudit')},
                 geometry_sha256=sha(v.build/'geometry-lib/libSpatialDomains.so'),cpus=a.cpus,maximum_ranks=4,threads=THREADS,results=v.results)
    write_json(v.output/'evidence.json',receipt)


if __name__=='__main__':
    main()
