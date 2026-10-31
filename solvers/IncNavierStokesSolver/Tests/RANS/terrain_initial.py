"""Prepare and independently audit a wall-compatible hill initial condition."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def prepare(session,build,cpu=None):
    session,build = Path(session).resolve(),Path(build).resolve()
    initial = session.parent/'initial.fld'
    command = (["taskset","-c",str(cpu)] if cpu is not None else [])+[str(build/'RANSWallFieldAudit'),str(session),'--terrain-initial',str(initial)]
    with (session.parent/'initial.log').open('w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    tree = ET.parse(session);root = tree.getroot()
    group = root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
    group.clear();group.attrib['NAME'] = 'InitialConditions'
    for variable in root.findall('CONDITIONS/VARIABLES/V'):
        ET.SubElement(group,'F',VAR=variable.text,FILE=initial.name)
    tree.write(session,encoding='utf-8',xml_declaration=True)
    result = audit(str(initial)+'.csv',root)
    return dict(command=command,field_sha256=hashlib.sha256(initial.read_bytes()).hexdigest(),**result)


def audit(path,root):
    p = {e.text.split('=',1)[0].strip():float(e.text.split('=',1)[1]) for e in root.findall('CONDITIONS/PARAMETERS/P')}
    dim = int(root.find('GEOMETRY').attrib['DIM'])
    amplitude = p.get('RANSTerrainHeight',0.) if next((s.attrib['VALUE'] for s in root.findall('CONDITIONS/SOLVERINFO/I') if s.attrib['PROPERTY']=='RANSWallGeometry'),'Plane') == 'CosineHill' else 0.
    H,N,star,kap,matching = (p[k] for k in ('H','RANSStripPointN','Ustar','Kappa','MatchingVelocity'))
    top,half,centre = p['Z'],p.get('RANSTerrainHalfLength',.1),p.get('RANSTerrainCentre',0.)
    def surface(x):
        a = (x-centre)/half
        if abs(a) >= 1: return 0.,0.,0.
        return amplitude/2*(1+math.cos(math.pi*a)),-amplitude*math.pi/(2*half)*math.sin(math.pi*a),-amplitude*(math.pi/half)**2/2*math.cos(math.pi*a)
    def speed(e): return matching+star/kap*math.log(e/H)
    def primitive(e): return matching*e+star/kap*e*(math.log(e/H)-1)
    def base(x,z):
        h,hp,_ = surface(x);gap=top-h
        e=z-h*((top-z)/gap)**2
        return e,speed(e)*(-hp*(top-z)**2*(top+h)/gap**3),speed(e)*(1+2*h*(top-z)/gap**2)
    mismatch = normal = reconstruction = 0.;min_k=min_eps=float('inf')
    for row in csv.DictReader(Path(path).open()):
        name = row.pop('field');row={k:float(v) for k,v in row.items()}
        x,z,foot,eta = row['x'],row['y' if dim==2 else 'z'],row['foot_x'],row['distance']
        h,hp,hpp=surface(foot);length=math.hypot(1.,hp)
        tx,tz,nx,nz=1/length,hp/length,-hp/length,1/length;c=hpp/length**3
        r=max(0.,min(1.,(eta-N)/N));chi=1-10*r**3+15*r**4-6*r**5
        derivative=-30*r*r*(1-r)**2/N
        b=base(x,z)
        a=chi*speed(eta)+derivative*(primitive(eta)-primitive(b[0]))
        gx=(1-chi)*b[1]+a*nx;gz=(1-chi)*b[2]+a*nz
        expected={'u':gz,('v' if dim==2 else 'w'):-gx,'k':p['Qstar']*star**2,'p':0.}
        if dim==3: expected['v']=0.
        value=expected.get(name,star**3/kap*(chi/eta+(1-chi)/b[0]))
        mismatch=max(mismatch,abs(row['input']-value)/max(1.,abs(value)))
        if not math.isfinite(row['input']) or not math.isfinite(row['projected']): raise RuntimeError('Nonfinite terrain initial field')
        if name=='k': min_k=min(min_k,row['projected'])
        if name in ('eps','eps_tilde'): min_eps=min(min_eps,row['projected'])
        if chi>0:
            reconstruction=max(reconstruction,abs(x-foot-eta*nx),abs(z-h-eta*nz))
            if abs(eta-H)<1e-9: normal=max(normal,abs(gz*nx-gx*nz))
    if mismatch>5e-10 or reconstruction>1e-8 or normal>1e-8 or min_k<=0 or min_eps<=0:
        raise RuntimeError(f'Terrain initial audit failed: mismatch={mismatch}, reconstruction={reconstruction}, normal={normal}, min_k={min_k}, min_eps={min_eps}')
    return dict(maximum_relative_mismatch=mismatch,normal_velocity_raw_max=normal,
                normal_coordinate_reconstruction_max=reconstruction,min_projected_k=min_k,min_projected_epsilon=min_eps)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('session',type=Path);parser.add_argument('--build',type=Path,required=True)
    args=parser.parse_args();print(prepare(args.session,args.build))
