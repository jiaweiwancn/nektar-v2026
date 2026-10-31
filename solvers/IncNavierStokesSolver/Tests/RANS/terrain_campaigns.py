"""Independent transition, restart and hill-pilot acceptance campaigns.

Called by validate_terrain.py; all solver work remains in its bounded Linux
runner. A short pilot verifies functionality, not an experimental prediction.
"""
import csv
import math
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import xml.etree.ElementTree as ET
from generate_m4 import info, parameter
from generate_m5 import reference


def parameters(root):
    return {p.text.split('=',1)[0].strip():float(p.text.split('=',1)[1])
            for p in root.findall('CONDITIONS/PARAMETERS/P')}


def settings(root):
    return {p.attrib['PROPERTY']:p.attrib['VALUE'] for p in root.findall('CONDITIONS/SOLVERINFO/I')}


def read_points(directory):
    result = []
    for path in sorted(Path(directory).glob('neutral.terrain-points-r*.csv')):
        with path.open() as stream:
            result += [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream) if float(row['phase']) == 1]
    return result


def point_oracle(directory):
    """Reconstruct mask, sensor norm, requested blend and both epsilon values."""
    root = ET.parse(Path(directory)/'neutral.xml').getroot()
    p, s = parameters(root), settings(root)
    original = s.get('RANSWallModel') != 'AWFEquilibrium'
    switched = s.get('RANSAWFWeight') == 'Switched'
    errors = dict(mask=0.,departure=0.,requested_chi=0.,epsilonP=0.,epsilonH=0.,meanDissipation=0.,roughness=0.)
    active = intermediate = inactive = 0
    for row in read_points(directory):
        x = row['x']
        mask = 1.
        if s.get('RANSAWFSpatialWeight') == 'Ramp':
            r = max(0.,min(1.,(x-p['RANSAWFTransitionStart'])/(p['RANSAWFTransitionEnd']-p['RANSAWFTransitionStart'])))
            mask = 1-10*r**3+15*r**4-6*r**5
        errors['mask'] = max(errors['mask'],abs(mask-row['mask']))
        flow = p.get('RANSAWFFixedWeight',1.)
        if switched and row['sensors_evaluated']:
            suffixes = ('U','DeltaU','K','Epsilon','Pressure','Convection','Normal')
            departure = math.sqrt(sum(p.get('RANSAWFWeight'+name,1.)*row['r_'+name]**2 for name in suffixes))
            errors['departure'] = max(errors['departure'],abs(departure-row['departure'])/max(1.,departure))
            a,b = p['RANSAWFOn'],p['RANSAWFOff']
            r = max(0.,min(1.,(departure**2-a*a)/(b*b-a*a)))
            flow = 1-10*r**3+15*r**4-6*r**5
        requested = 0. if original or mask == 0 else mask*flow
        errors['requested_chi'] = max(errors['requested_chi'],abs(requested-row['requested_chi']))
        if row['chi'] > 1e-10: active += 1
        if 1e-10 < mask < 1-1e-10: intermediate += 1
        if mask < 1e-10: inactive += 1
        if s.get('RANSWallModel') in ('AWFOriginal','AWFEquilibrium'):
            k,nu,chi = row['kP'],p['Kinvis'],row['chi']
            cutoff0 = 5.1*nu/math.sqrt(k)
            ref = reference(k,row['H'],row['P'],row['roughness'],nu,p['Qstar'],p['Kappa']) if chi else None
            for key,z in [('epsilonP',row['P']),('epsilonH',row['H'])]:
                value = (1-chi)*k**1.5/(2.55*max(z,cutoff0))
                if chi: value += chi*ref['numerator']/max(z,ref['cutoff'])
                actual = row['epsP' if key == 'epsilonP' else 'epsH']
                errors[key] = max(errors[key],abs(value-actual)/max(1.,abs(value)))
            def metric_integral(cutoff):
                H,c=row['H'],row['curvature']
                return row['V']/cutoff if cutoff >= H else 1+math.log(H/cutoff)-c*(H-cutoff/2)
            mean=(1-chi)*k**1.5/2.55*metric_integral(cutoff0)
            if chi: mean+=chi*ref['numerator']*metric_integral(ref['cutoff'])
            mean/=row['V']
            errors['meanDissipation']=max(errors['meanDissipation'],abs(mean-row['dissipation'])/max(1.,abs(mean)))
        query = row.get('property_x',x)
        expression = root.find("CONDITIONS/FUNCTION[@NAME='RANSWallProperties']/E")
        if expression is not None:
            # These campaign fixtures deliberately use simple, independently
            # evaluated expressions rather than the solver's native interpreter.
            text = expression.attrib['VALUE']
            if 'tanh' in text:
                value = (.0002 if expression.attrib['VAR'] == 'z0' else .002)+(.0001 if expression.attrib['VAR'] == 'z0' else .001)*(1+math.tanh(query/.1))/2
            elif 'sign' in text:
                value = .001+.001*(1+(1 if query > -.1 else -1 if query < -.1 else 0))/2
            else:
                value = p.get('WallKs',float(text)) if text != 'WallKs' else p['WallKs']
            actual = row['z0' if expression.attrib['VAR'] == 'z0' else 'roughness']
            errors['roughness'] = max(errors['roughness'],abs(value-actual))
    if max(errors.values()) > 5e-10:
        raise RuntimeError('Independent wall-point oracle failed: '+str(errors))
    return dict(errors=errors,active_points=active,ramp_points=intermediate,inactive_points=inactive)


def sample_oracle(directory):
    """Interpolate saved primary fields independently at physical normal N.

    These MPI tests use an axis-aligned plane. At each wall quadrature x/span,
    independent one-dimensional Lagrange interpolation along a mesh column
    suffices; it does not call the solver's point search or interpolation API.
    """
    root = ET.parse(Path(directory)/'neutral.xml').getroot()
    dimension = int(root.find('GEOMETRY').attrib['DIM'])
    with (Path(directory)/'field.csv').open() as stream:
        fields = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    points = read_points(directory);time = max(p['time'] for p in points)
    errors = dict(uN=0.,kN=0.,epsilonN=0.)
    for point in (p for p in points if p['time'] == time):
        columns = {}
        for row in fields:
            if abs(row['x']-point['x']) < 1e-11 and (dimension == 2 or abs(row['y']-point['span']) < 1e-11):
                columns.setdefault(int(row['element']),[]).append(row)
        vertical = 'y' if dimension == 2 else 'z'
        candidates = [values for _,values in sorted(columns.items())
            if min(r[vertical] for r in values)-1e-11 <= point['N'] <= max(r[vertical] for r in values)+1e-11]
        if not candidates: raise RuntimeError('No independent N column found')
        column = sorted(candidates[0],key=lambda r:r[vertical])
        z = [r[vertical] for r in column]
        weights = [math.prod((point['N']-v)/(z[i]-v) for j,v in enumerate(z) if j != i) for i in range(len(z))]
        for key,field in [('uN','u'),('kN','k'),('epsilonN','eps' if 'eps' in column[0] else 'eps_tilde')]:
            value = sum(w*r[field] for w,r in zip(weights,column))
            errors[key] = max(errors[key],abs(point[key]-value)/max(1.,abs(value)))
    if max(errors.values()) > 5e-10: raise RuntimeError('Independent physical N-sample check failed: '+str(errors))
    return errors


def parallel(run, jobs, workers):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run,name,**options) for name,options in jobs]
        try:
            for future in as_completed(futures): future.result()
        except Exception:
            for future in futures: future.cancel()
            raise


def transition_campaign(run,cases,comparisons,output,workers,compare):
    def verification(root):
        parameter(root,'TimeStep',2e-5)
        parameter(root,'RANSStageAudit',1)
        parameter(root,'RANSNonlinearTolerance',1e-8)
        parameter(root,'RANSPressureTolerance',1e-8)
    jobs = []
    for closure in ('standard','cls'):
        for mode in ('fixed','switched'):
            jobs.append((closure+'_'+mode+'_floor',dict(kind='hill',model=closure+'_'+mode,
                empty=True,steps=1,ranks=4,modify=verification)))
        def inactive(root):
            verification(root)
            parameter(root,'RANSAWFTransitionStart',-2)
            parameter(root,'RANSAWFTransitionEnd',-1)
            parameter(root,'RANSAWFMinimumRe',1e12)
        jobs.append((closure+'_zero_hill',dict(kind='hill',model=closure+'_fixed',steps=1,ranks=4,modify=inactive)))
        jobs.append((closure+'_original_hill',dict(kind='hill',model=closure+'_original',steps=1,ranks=4,modify=verification)))
        jobs.append((closure+'_smooth_property',dict(kind='hill',model=closure+'_fixed',properties=True,steps=1,ranks=4,modify=verification)))
        def step_property(root):
            verification(root)
            group = ET.SubElement(root.find('CONDITIONS'),'FUNCTION',NAME='RANSWallProperties')
            ET.SubElement(group,'E',VAR='roughness',VALUE='0.001+0.001*(1+sign(x+0.1))/2')
        jobs.append((closure+'_step_property',dict(kind='hill',model=closure+'_original',steps=1,ranks=4,modify=step_property)))
        for ranks in (1,2,4):
            jobs.append((closure+'_mpi'+str(ranks),dict(model=closure+'_fixed',nx=1,layers=8,steps=2,ranks=ranks)))
        def smooth_wall(root):
            group = ET.SubElement(root.find('CONDITIONS'),'FUNCTION',NAME='RANSWallProperties')
            ET.SubElement(group,'E',VAR='roughness',VALUE='0')
        jobs.append((closure+'_smooth_wall',dict(model=closure+'_original',nx=1,layers=6,steps=1,dt=.0002,modify=smooth_wall)))
    def span_perturbation(root):
        parameter(root,'N',parameters(root)['RANSStripPointN'])
        parameter(root,'SpanAmplitude',1e-4)
        parameter(root,'SpanWave',2*math.pi/.032)
        r='((abs(z/N-1)-abs(z/N-2)+1)/2)'
        envelope=f'(1-10*{r}^3+15*{r}^4-6*{r}^5)'
        derivative=f'(-30*{r}^2*(1-{r})^2/N)'
        group=root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
        for node in group:
            if node.attrib['VAR']=='v':node.attrib['VALUE']=f'SpanAmplitude*({envelope}+z*{derivative})*cos(SpanWave*y)'
            if node.attrib['VAR']=='w':node.attrib['VALUE']=f'SpanAmplitude*SpanWave*z*{envelope}*sin(SpanWave*y)'
        exact=root.find("CONDITIONS/FUNCTION[@NAME='ExactSolution']")
        if exact is not None:root.find('CONDITIONS').remove(exact)
    for closure in ('standard','cls'):
        jobs.append((closure+'_span_continuity',dict(model=closure+'_fixed',dimension=3,nx=1,layers=6,
            dt=.0002,steps=2,ranks=4,modify=span_perturbation)))
    jobs.append(('log_smooth_property',dict(kind='hill',model='standard',properties=True,steps=1,ranks=4,modify=verification)))
    def legacy_properties(root):
        group=ET.SubElement(root.find('CONDITIONS'),'FUNCTION',NAME='RANSWallProperties')
        ET.SubElement(group,'E',VAR='roughness',VALUE='0.001+0.001*(1+sign(x-0.2))/2')
    for closure in ('standard','cls'):
        jobs.append((closure+'_legacy_property',dict(model=closure+'_original',legacy=True,dt=.0002,steps=1,modify=legacy_properties,
            expected_rejection='Spatial AWF transitions/properties require RANSWallGeometry=Plane')))
    def legacy_ramp(root):
        info(root,'RANSAWFSpatialWeight','Ramp')
        parameter(root,'RANSAWFTransitionStart',.08)
        parameter(root,'RANSAWFTransitionEnd',.28)
    jobs.append(('legacy_ramp',dict(model='standard_fixed',legacy=True,steps=1,modify=legacy_ramp,
        expected_rejection='Spatial AWF transitions/properties require RANSWallGeometry=Plane')))
    parallel(run,jobs,workers)
    for closure in ('standard','cls'):
        for mode in ('fixed','switched'):
            name = closure+'_'+mode+'_floor'
            receipt = point_oracle(output/name)
            if min(receipt[k] for k in ('active_points','ramp_points','inactive_points')) <= 0:
                raise RuntimeError('Spatial gate was not sampled in every region: '+name)
            comparisons[name+'_oracle'] = receipt
        errors = compare(output/(closure+'_zero_hill/field.csv'),output/(closure+'_original_hill/field.csv'))
        comparisons[closure+'_zero_is_original'] = errors
        comparisons[closure+'_zero_policy_oracle'] = point_oracle(output/(closure+'_zero_hill'))
        if max(errors.values()) > 1e-6: raise RuntimeError('Inactive correction differs from original AWF beyond the declared field tolerance')
        for ranks in (2,4):
            errors = compare(output/(closure+'_mpi1/field.csv'),output/(closure+'_mpi'+str(ranks)+'/field.csv'))
            comparisons[closure+'_mpi'+str(ranks)] = errors
            if max(errors.values()) > 1e-6: raise RuntimeError('Rank-count field agreement failed: '+str(errors))
        if min(cases[closure+'_mpi4']['remote_samples'],cases[closure+'_mpi4'].get('remote_without_wall_samples',0)) <= 0:
            raise RuntimeError('No remote normal samples exercised: '+closure)
        comparisons[closure+'_independent_remote_N'] = sample_oracle(output/(closure+'_mpi4'))
        for kind in ('smooth','step'):
            comparisons[closure+'_'+kind+'_property'] = point_oracle(output/(closure+'_'+kind+'_property'))
        foot = [r for r in read_points(output/(closure+'_step_property')) if abs(r['x']+.1)<1e-10]
        properties = {round(r['roughness'],8) for r in foot}
        if properties != {.001,.002} or any((r['property_x']>-.1) != (r['curvature']>1) for r in foot):
            raise RuntimeError('Terrain property jump did not use the two interior foot traces: '+closure)
        comparisons[closure+'_one_sided_foot_properties'] = dict(shared_node_properties=sorted(properties),
            maximum_value_error=max(abs(r['roughness']-(.002 if r['curvature']>1 else .001)) for r in foot))
    comparisons['log_smooth_property'] = point_oracle(output/'log_smooth_property')

    for closure in ('standard','cls'):
        rows=read_points(output/(closure+'_span_continuity'))
        if max(abs(r['spanP']) for r in rows)<1e-6 or max(abs(r['wP']) for r in rows)<1e-6:
            raise RuntimeError('Nonzero spanwise strip continuity was not exercised: '+closure)

    def restart(root,field):
        group = root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
        group.clear();group.attrib['NAME']='InitialConditions'
        for node in root.findall('CONDITIONS/VARIABLES/V'):
            ET.SubElement(group,'F',VAR=node.text,FILE=str(field))
        info(root,'RANSRestartMode','History')
    def history_settings(root):
        info(root,'RANSAWFSpatialWeight','Ramp')
        parameter(root,'RANSAWFTransitionStart',.08)
        parameter(root,'RANSAWFTransitionEnd',.28)
        parameter(root,'RANSNonlinearTolerance',1e-9)
        parameter(root,'RANSPressureTolerance',1e-9)
        parameter(root,'WallKs',.002)
        group = ET.SubElement(root.find('CONDITIONS'),'FUNCTION',NAME='RANSWallProperties')
        ET.SubElement(group,'E',VAR='roughness',VALUE='WallKs')
    for closure in ('standard','cls'):
        options = dict(model=closure+'_fixed',nx=1,layers=6,dt=.0002,boundary='matched')
        run(closure+'_history_full',steps=4,modify=history_settings,**options)
        run(closure+'_history_prefix',steps=2,modify=history_settings,**options)
        field = output/(closure+'_history_prefix/neutral.fld')
        def resume(root): history_settings(root);restart(root,field)
        run(closure+'_history_resume',steps=2,ranks=4,modify=resume,**options)
        errors = compare(output/(closure+'_history_full/field.csv'),output/(closure+'_history_resume/field.csv'))
        comparisons[closure+'_history_restart'] = errors
        if max(errors.values()) > 1e-6: raise RuntimeError('Terrain full-history restart differs: '+str(errors))
        for label,key,value in [('mask','RANSAWFTransitionEnd',.3),('property','WallKs',.0021),('normal','RANSNormalAugmentation',11)]:
            def mismatch(root,key=key,value=value): resume(root);parameter(root,key,value)
            run(closure+'_history_bad_'+label,steps=1,modify=mismatch,
                expected_rejection='model/configuration',**options)
        def moved_mesh(root):
            resume(root)
            for vertex in root.findall('GEOMETRY/VERTEX/V'):
                xyz = list(map(float,vertex.text.split()));xyz[0] *= 1.0001
                vertex.text = ' '.join(format(v,'.17g') for v in xyz)
        run(closure+'_history_bad_coordinates',steps=1,modify=moved_mesh,
            expected_rejection='element coordinates/order mismatch',**options)
    negative = [
        ('bad_ramp',lambda root:(info(root,'RANSAWFSpatialWeight','Ramp'),parameter(root,'RANSAWFTransitionStart',.2),parameter(root,'RANSAWFTransitionEnd',.1)),'finite begin < end',False),
        ('over_hill',lambda root:parameter(root,'RANSAWFTransitionEnd',-.09),'strictly upstream',True),
        ('thin_layer_limit',lambda root:parameter(root,'RANSMaximumCurvatureRatio',.15),'interface must be a physical normal offset',True),
        ('native_backend',lambda root:info(root,'RANSLinearSolver','Native'),'requires RANSLinearSolver=MatrixFree',False),
        ('active_bad_re',lambda root:parameter(root,'RANSAWFMinimumRe',1e12),'fails cutoff/high-Re preflight',False),
        ('moving_property',lambda root:ET.SubElement(ET.SubElement(root.find('CONDITIONS'),'FUNCTION',NAME='RANSWallProperties'),'E',VAR='roughness',VALUE='0.002+0.0001*t'),'independent of time',False),
    ]
    for name,modify,message,curved in negative:
        run(name,kind='hill' if curved else 'plane',model='standard_fixed',steps=1,
            modify=modify,expected_rejection=message)


def pilot_campaign(run,cases,comparisons,output,workers,compare):
    jobs = [('log_hill',dict(kind='hill',model='standard',dt=2e-5,steps=2,ranks=4))]
    for closure in ('standard','cls'):
        jobs.append((closure+'_original_hill',dict(kind='hill',model=closure+'_original',dt=2e-5,steps=2,ranks=4)))
        for width in (1,2,4,8):
            for empty in (False,True):
                name = closure+'_width'+str(width)+('_floor' if empty else '_hill')
                jobs.append((name,dict(kind='hill',model=closure+'_fixed',dt=2e-5,steps=2,ranks=4,width_ratio=width,empty=empty)))
        jobs.append((closure+'_molecular_hill',dict(kind='hill',model=closure+'_fixed',dt=2e-5,steps=2,ranks=4,diffusion='Molecular')))
        jobs.append((closure+'_3d_hill',dict(kind='hill',model=closure+'_fixed',dt=2e-5,steps=2,ranks=4,dimension=3)))
    jobs.sort(key=lambda job:job[1].get('dimension',2)!=3)
    parallel(run,jobs,workers)
    for closure in ('standard','cls'):
        study = []
        for width in (1,2,4,8):
            pair = {}
            for terrain in ('hill','floor'):
                name = closure+'_width'+str(width)+'_'+terrain
                comparisons[name+'_oracle'] = point_oracle(output/name)
                points = read_points(output/name)
                time = max(p['time'] for p in points)
                points = [p for p in points if p['time'] == time]
                foot = min(abs(p['x']+.1) for p in points)
                incident = [p for p in points if abs(p['x']+.1) <= foot+1e-10]
                wake = [p for p in points if p['x'] >= .1]
                pair[terrain] = dict(incident_uH=sum(p['uH'] for p in incident)/len(incident),
                    incident_kP=sum(p['kP'] for p in incident)/len(incident),
                    incident_epsH=sum(p['epsH'] for p in incident)/len(incident),
                    incident_tau=sum(p['traction_x'] for p in incident)/len(incident),
                    wake_min_traction=min(p['traction_x'] for p in wake),
                    negative_traction_points=sum(p['traction_x'] < 0 for p in wake),
                    chi_max=max(p['chi'] for p in points),
                    maximum_weak_divergence=cases[name]['maximum_weak_divergence'])
                if width != 1:
                    pair[terrain]['complete_field_change_from_width1'] = compare(
                        output/(closure+'_width1_'+terrain)/'field.csv',output/name/'field.csv')
                # Whole-height incident/crest/wake profiles on mesh-aligned
                # physical columns. Collapse element-boundary duplicates.
                with (output/name/'field.csv').open() as stream:
                    fields = [{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
                profiles = []
                for station in (-.1,0.,.2,.6):
                    selected = [r for r in fields if abs(r['x']-station) < 1e-10]
                    groups = {}
                    for r in selected: groups.setdefault(round(r['y'],12),[]).append(r)
                    for z, group in sorted(groups.items()):
                        profiles.append(dict(x=station,z=z,**{k:sum(r[k] for r in group)/len(group)
                            for k in ('u','v','k','eps' if closure == 'standard' else 'eps_tilde')}))
                if not profiles: raise RuntimeError('No physical station profiles were exported')
                with (output/name/'profiles.csv').open('w') as stream:
                    writer = csv.DictWriter(stream,fieldnames=list(profiles[0]))
                    writer.writeheader();writer.writerows(profiles)
            study.append(dict(width_over_hill_height=width,**pair))
        comparisons[closure+'_width_study'] = study
