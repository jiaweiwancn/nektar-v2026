#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Specialize only CT11's two request literals for an admitted CT14 cohort.

The CT11 file and every scientific check stay unchanged. No history, state,
time or metric is renamed. The AST edit replaces exactly the required-dt list
and required-member set; accepted coverage still uses each actual dt.
"""
import ast,copy,inspect
from functools import lru_cache
import ct_accepted11_target as original

COHORTS=((.001,.0005),(.00025,.000125))

@lru_cache(None)
def classifier(dts):
    if dts not in COHORTS:raise ValueError('unadmitted CT14 timestep cohort')
    tree=ast.parse(inspect.getsource(original.target_pair));counts=[0,0]
    class RequestOnly(ast.NodeTransformer):
        def visit_List(self,node):
            if ast.literal_eval(node)==[.001,.0005]:
                counts[0]+=1;return ast.copy_location(ast.Name(id='CT14_REQUIRED_DT',ctx=ast.Load()),node)
            return self.generic_visit(node)
        def visit_Set(self,node):
            if ast.literal_eval(node)=={'0.001','0.0005'}:
                counts[1]+=1;return ast.copy_location(ast.Name(id='CT14_REQUIRED_NAMES',ctx=ast.Load()),node)
            return self.generic_visit(node)
    # Only literal nodes are candidates; other lists/sets retain their AST.
    old_list,old_set=RequestOnly.visit_List,RequestOnly.visit_Set
    def safe(method):
        def visit(self,node):
            try:ast.literal_eval(node)
            except (ValueError,TypeError):return self.generic_visit(node)
            return method(self,node)
        return visit
    RequestOnly.visit_List=safe(old_list);RequestOnly.visit_Set=safe(old_set)
    tree=RequestOnly().visit(tree)
    if counts!=[1,1]:raise ValueError('CT11 request template changed; no automatic scientific edit')
    namespace=dict(vars(original),CT14_REQUIRED_DT=list(dts),CT14_REQUIRED_NAMES={str(x) for x in dts})
    exec(compile(ast.fix_missing_locations(tree),'<CT14 request-only CT11 specialization>','exec'),namespace)
    return namespace['target_pair']

def cohort_pair(body,dts):return classifier(tuple(dts))(body)

def equivalence_controls(body):
    """Actual first-pair parity plus explicit synthetic request/physical negatives."""
    if original.target_pair(body)!=cohort_pair(body,COHORTS[0]):raise ValueError('actual first-pair classifier changed')
    good=copy.deepcopy(body);good['temporal']['Pi']['difference']['weighted_RMS']=0.
    if original.target_pair(good)['status']!='pass':raise ValueError('positive scientific control incomplete')
    controls={}
    for name,mutate in (
        ('positive',lambda q:None),
        ('absolute_pressure',lambda q:q['temporal']['Pi']['difference'].update(weighted_RMS=1.)),
        ('field_drift',lambda q:q['members']['0.001']['accepted_metrics'][-1]['drift']['U'].update(weighted_RMS=1.)),
        ('mass',lambda q:q['members']['0.001']['accepted_metrics'][-1]['mass'].update(outward_net=1.)),
        ('traction',lambda q:q['members']['0.001']['accepted_metrics'][-1]['traction'].update(normalized_interface_traction_RMS=1.)),
        ('domain',lambda q:q['members']['0.001']['domain'].update(status='fail')),
        ('factor',lambda q:q['members']['0.001']['actual_stage_checks'][0].update(predictor_factor_nonzeros=40000001)),
        ('initial',lambda q:q.update(same_actual_initial=False)),
        ('coverage',lambda q:q['members']['0.001']['accepted_metrics'].pop()),
        ('budget',lambda q:q['members']['0.001']['budgets']['budgets']['native']['equations']['inner_k'].update(storage_backward_difference=1000.)),
        ('extra_solve',lambda q:q['members']['0.001'].update(extra_solves=1))):
        q=copy.deepcopy(good);mutate(q);a=original.target_pair(q);b=cohort_pair(q,COHORTS[0])
        if a!=b or (a['status']=='pass')!=(name=='positive'):raise ValueError('classifier parity/discrimination:'+name)
        controls[name]=a
    synthetic=copy.deepcopy(good);synthetic['synthetic_request_control']=True;synthetic['required_dt']=list(COHORTS[1]);synthetic['members']={}
    for old,dt in zip(good['members'].values(),COHORTS[1],strict=True):
        row=copy.deepcopy(old);row['dt']=dt;row['accepted_metrics']=[dict(copy.deepcopy(old['accepted_metrics'][-1]),time=round(j*dt,12)) for j in range(round(.002/dt)+1)];row['budgets']['budget_binding']['dt']=dt
        for control in row['budgets']['budgets'].values():control['binding']['dt']=dt
        synthetic['members'][str(dt)]=row
    if cohort_pair(synthetic,COHORTS[1])['status']!='pass':raise ValueError('new request positive failed')
    if cohort_pair(body,COHORTS[1])['status']!='fail':raise ValueError('foreign old cohort accepted')
    synthetic['members']['0.000125']['accepted_metrics'].pop()
    if cohort_pair(synthetic,COHORTS[1])['status']!='fail':raise ValueError('new actual-dt coverage negative accepted')
    return dict(status='pass',actual_first_pair=original.target_pair(body),same_cohort_parity=controls,new_request_positive='synthetic request/metrics only,not native-history qualification',foreign_old_cohort='rejected',new_coverage_negative='rejected',AST_edits='exactly required-dt list and required-member set; every scientific check,scale,tolerance and T=.002 unchanged')
