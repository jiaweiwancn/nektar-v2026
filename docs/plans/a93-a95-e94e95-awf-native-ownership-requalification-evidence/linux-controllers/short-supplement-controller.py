import pathlib,subprocess,os,time
r=pathlib.Path(__file__).parent
while not (r/'ownership-short-qualification/results.json').exists():time.sleep(5)
source=r/'source/solvers/IncNavierStokesSolver/Tests/RANS/run_rans_ownership_short.py'
cmd=['python3','-B',str(source),'--root',str(r),'--geometry-receipt',str(r/'geometric-regressions-v1/results.json'),'--concurrency','2','--output-label','ownership-short-qualification-v2','--audit-prefix','short-audits-v2','--reuse-executions',str(r/'ownership-short-qualification')]
with (r/'short-supplement.log').open('x') as f:c=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
(r/'short-supplement-exit.txt').write_text(str(c.returncode))
