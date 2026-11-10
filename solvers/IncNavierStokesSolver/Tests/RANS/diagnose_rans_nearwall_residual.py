#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One targeted frozen experiment: localize the initial wall dual-residual floor.

No solver advancement, source correction or field resetting is performed.
"""
import argparse,os,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
from rans_preservation_oracle import Reference,weak_balance,representation,physical_grid
from qualify_rans_preservation import load,probe,Q
from run_a93_a95_e94e95_matrix import case_xml

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'nearwall-frozen-diagnosis';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));rows=[]
 for wall in ('EP','Accepted'):
  try:
   ref=Reference('E94E95','Natural',wall);knots=ref.knots();chosen=list(zip(knots[:-1],knots[1:],strict=True))[:8];fine=sorted(set(knots+[(x+y)/2 for x,y in chosen]));folder=out/wall
   case_xml(folder,'E94E95','Natural',wall,steps=0,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=fine,degree=5,quadrature=Q[1],boundary_policy='FrozenReferenceV1')
   base=load(r/'campaigns'/('E-N-'+wall+'-h0p5')/'initial-metrics.json');interface=load(r/'campaigns'/('E-N-'+wall+'-Interface8p5')/'initial-metrics.json');grid,weights=physical_grid(sorted(set(knots+fine)),10,8)
   full=probe(r/'runtime/RANSPreservationProbe',folder,dict(operation='neutral',grid=grid,export_basis=True));o=full['observation'];weak=weak_balance(ref,o);rep=representation(ref,o,weights);audit=o['audit_alpha_zero'];b=base['residual'];h=interface['residual']
   rows.append(dict(wall=wall,status='pass' if weak['status']=='pass' else 'fail',scope='read-only initialization, zero accepted steps',dof=sum(o['sizes']),split_intervals=chosen,analytical_weak_balance=weak,representation=rep,alpha_zero=audit,
    base_vs_interface_vs_nearwall=dict(momentum_x=[b['coupled_norms'][0],h['coupled_norms'][0],audit['coupled_norms'][0]],inner_k=[b['native_k_inner_outer_shared_interface_dual'][0],h['native_k_inner_outer_shared_interface_dual'][0],audit['native_k_inner_outer_shared_interface_dual'][0]],outer_epsilon=[b['coupled_norms'][3],h['coupled_norms'][3],audit['coupled_norms'][3]])))
  except Exception as e:rows.append(dict(wall=wall,status='fail',failure=repr(e)))
  write_json(out/'results.json',dict(status='running',checks=rows))
 status='pass' if len(rows)==2 and all(x['status']=='pass' for x in rows) else 'fail';write_json(out/'results.json',dict(status=status,checks=rows,classification='residual reductions are measured independently of the field-error pilot gates; no new evolution qualification'));print(status);return 0 if status=='pass' else 1
if __name__=='__main__':sys.exit(main())
