#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reproducible H26 verification and interleaved complete-step matrix."""
import argparse,copy,hashlib,json,math,os,statistics,struct,sys
from pathlib import Path
import xml.etree.ElementTree as E
from run import run,THREADS
os.environ.update(THREADS)

def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def write(path,value):
    assert not path.exists(),str(path)
    path.write_text(json.dumps(value,indent=2)+'\n')
def load(path):return json.loads(path.read_text())

class Campaign:
    def __init__(self,root,version):
        self.root=root;self.version=version;self.tools=root/f'source-v{version}/solvers/IncNavierStokesSolver/Tests/RANS/h26';self.runtime=root/f'runtime-v{version}'
        self.matrix=root/f'science/matrix-v{version}'
    def launch(self,label,case,request,stock=False,seconds=180):
        label+=f'-v{self.version}';p=self.root/'inputs'/f'{label}.json';out=self.root/'science'/f'{label}.json';receipt=self.root/'operations'/label/'receipt.json'
        if receipt.exists():
            assert load(receipt)['status']=='PASS',f'failed operation requires a NEW explicit label: {label}'
            assert load(p)==request
            return load(out)
        write(p,request);program='WJStockFlowDriver'if stock else'WJNativeFlowDriver'
        status=run(self.root,label,[str(self.runtime/program),str(case),str(p),str(out)],'science',seconds)
        assert status==0,f'{label} failed; preserve log/result before repair'
        return load(out)
    def analysis(self,label,arguments,seconds=120):
        label+=f'-v{self.version}';receipt=self.root/'operations'/label/'receipt.json'
        if receipt.exists():assert load(receipt)['status']=='PASS';return
        assert run(self.root,label,[sys.executable,'-B',*map(str,arguments)],'analysis',seconds)==0
    def prepare(self):
        assert not self.matrix.exists();self.matrix.mkdir();base=self.root/'science/profile-v3/case.xml';template=E.parse(base).getroot()
        geometry=E.tostring(template.find('GEOMETRY'));assert len(template.find('GEOMETRY/ELEMENT'))==700
        nu=6.423173803526447e-5;mu=nu+.3819178635467331*.051*19.9125/2
        X='((x+10)/25)';Z='((y-0.0875)/19.9125)';u=f'exp(-t)*{X}^2*(1-{X})*(2*{Z}-3*{Z}^2)';v=f'-exp(-t)*(19.9125/25)*(2*{X}-3*{X}^2)*{Z}^2*(1-{Z})';p=f'0.01*exp(-2*t)*({X}^2+{Z}^2)'
        config={}
        for name,model,stock,stress,awf in (('Snu','nu',True,0,0),('Smu','constant',True,0,0),('Nnu','nu',False,0,0),('Nmu','constant',False,0,0),('Nspace','space',False,.001,1),('Ntime','time',False,.001,1)):
            x=copy.deepcopy(template);conditions=x.find('CONDITIONS')
            for item in conditions.find('SOLVERINFO'):
                if item.get('PROPERTY')=='SolverType':item.set('VALUE','VCSWeakPressure')
            for item in conditions.find('PARAMETERS'):
                key=item.text.split('=')[0].strip()
                if key in ('TimeStep','Kinvis','NumSteps','IO_InfoSteps','IO_CheckSteps'):item.text=f'{key} = '+str(dict(TimeStep=.005,Kinvis=mu if model=='constant'else nu,NumSteps=420,IO_InfoSteps=20,IO_CheckSteps=200)[key])
            for item in conditions.findall('BOUNDARYCONDITIONS/REGION/D'):
                item.set('VALUE',u if item.get('VAR')=='u'else v);item.set('USERDEFINEDTYPE','TimeDependent')
            for function in conditions.findall('FUNCTION'):
                if function.get('NAME')=='InitialConditions':
                    for item in function:item.set('VALUE',dict(u=u,v=v,p=p)[item.get('VAR')].replace('exp(-t)','1').replace('exp(-2*t)','1'))
            assert E.tostring(x.find('GEOMETRY'))==geometry
            case=self.matrix/f'{name}.xml';E.indent(x);case.write_text(E.tostring(x,encoding='unicode')+'\n')
            req=dict(mode='advance',model=model,profile='project',steps=420,warmup=20,dt=.005,awf=awf,stress=stress,d0=.0875,mu_ref=mu,checkpoint_stride=200,summary_stride=20,tolerance=1e-11)
            config[name]=dict(case=identity(case),request=req,stock=stock)
        # Independent closed-form identity for forcing time levels, including
        # BE startup; this is not an approximation to the source time level.
        weight_error=0.
        for rate in (1.,2.):
            rho=math.exp(-rate*.005);a=rho*rho/(2*rho-1);w=lambda n:a*rho**n+(rho-a)*.5**n
            assert abs(w(0)-rho)<1e-15
            for n in range(1,421):weight_error=max(weight_error,abs(2*w(n)-w(n-1)-rho**(n+1)))
        assert weight_error<1e-14
        write(self.matrix/'config.json',dict(cases=config,geometry_sha256=hashlib.sha256(geometry).hexdigest(),source_template=identity(base),mu_ref=mu,mu_ref_rule='molecular nu plus retained outer kappa*u_tau*r at the physical mid-height',turbulent_range_envelope=[0.,.3819178635467331*.051*19.9125],time_modulation=[.9,1.1],time_weight_identity_error=weight_error,interval=[0.,2.1],warmed_interval=[.1,2.1],warmup=20,accepted_per_block=420,retained_per_block=400,source=identity(self.root/f'source-v{self.version}.json')))
    def verify(self):
        root=self.root;fixture=root/'science/fixtures-v2'
        for curve,natural in ((0.,0),(0.,1),(.1,0),(.1,1)):
            tag=f'algebra-c{curve}-n{natural}';self.launch(tag,fixture/f'tiny-c{curve}-n{natural}.xml',dict(mode='algebra',curve=curve,awf=1,alpha=300))
            self.analysis(f'independent-c{curve}-n{natural}',[self.tools/'oracle.py','--native',root/'science'/f'{tag}-v{self.version}.json','--curve',curve,'--out',root/'science'/f'independent-c{curve}-n{natural}-v{self.version}.json'])
        for mode in ('adapter','nonlinear','formula'):self.launch(mode,fixture/'tiny-c0.1-n1.xml',dict(mode=mode,curve=.1))
        self.analysis('independent-continuous-formulas',[self.tools/'formula.py','--source',root/'science'/f'formula-v{self.version}.json','--out',root/'science'/f'independent-formulas-v{self.version}.json'])
        for j,dt in enumerate((.04,.02,.01)):
            self.launch(f'temporal-poly-{j}',fixture/'tiny-c0.0-n0.xml',dict(mode='advance',model='time',steps=round(.4/dt),dt=dt,awf=1,checkpoint_stride=0))
        self.launch('tolerance-tight',fixture/'tiny-c0.0-n0.xml',dict(mode='advance',model='time',steps=40,dt=.01,awf=1,tolerance=1e-13,checkpoint_stride=0))
        whole=self.launch('restart-whole',fixture/'tiny-c0.0-n0.xml',dict(mode='advance',steps=80,dt=.005,checkpoint_stride=0));first=self.launch('restart-first',fixture/'tiny-c0.0-n0.xml',dict(mode='advance',steps=37,dt=.005,checkpoint_stride=0));resumed=self.launch('restart-resumed',fixture/'tiny-c0.0-n0.xml',dict(mode='advance',steps=43,dt=.005,checkpoint_stride=0,restart=first['final_checkpoint']))
        assert whole['state_sha256']==resumed['state_sha256']and whole['previous_sha256']==resumed['previous_sha256']
        from fixtures import mesh
        folder=root/f'science/spatial-v{self.version}';folder.mkdir(exist_ok=True)
        for curve in (0.,.1):
            for n in (2,4,8):
                case=folder/f'c{curve}-n{n}.xml'
                if not case.exists():mesh(case,[i/n for i in range(n+1)],[i/n for i in range(n+1)],curve=curve,natural=True,quadrature=12)
                d=self.launch(f'spatial-c{curve}-n{n}',case,dict(mode='advance',steps=100,dt=.0001,curve=curve,solution='trigonometric',model='time',awf=1,checkpoint_stride=0));print('SPATIAL',curve,n,d['L2_errors'],d['H1_velocity_seminorm_error'],d['traction_errors_L2_and_weak_moments'],flush=True)
    def benchmark(self,pilot=False):
        config=load(self.matrix/'config.json')['cases']
        repetitions=1 if pilot else 3
        for repetition in range(repetitions):
            for name in ('Snu','Nnu','Smu','Nmu','Nspace','Ntime'):
                c=config[name];request=dict(c['request'])
                if pilot:request.update(steps=40,warmup=10,checkpoint_stride=20)
                d=self.launch(f'{"pilot"if pilot else"block"}-{repetition}-{name}',Path(c['case']['path']),request,c['stock'],120)
                print(name,repetition,statistics.mean(d['ordinary_samples']),min(d['ordinary_samples']),max(d['ordinary_samples']),flush=True)
        if not pilot:
            c=config['Ntime'];first=load(self.root/'science'/f'block-2-Ntime-v{self.version}.json');request=dict(c['request']);request.update(steps=20,warmup=0,restart=first['final_checkpoint'],checkpoint_stride=0);self.launch('restart-700',Path(c['case']['path']),request,False,120)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);p.add_argument('--mode',choices=('prepare','verify','pilot','benchmark'),required=True);a=p.parse_args();c=Campaign(a.root,a.version)
    if a.mode=='pilot':c.benchmark(True)
    else:getattr(c,a.mode)()
