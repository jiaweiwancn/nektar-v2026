#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check the exact bounded producer/observer source differences; no ABI guess."""
import argparse
import hashlib
from pathlib import Path
import subprocess
from wj_immutable import write_once


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--old',required=True);p.add_argument('--new',default='HEAD');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    git=lambda *v:subprocess.check_output(['git','-C',str(a.source),*v],text=True)
    prefix='solvers/IncNavierStokesSolver/EquationSystems/'
    diag=prefix+'RANSWallEvolutionDiagnostics.cpp';field=prefix+'RANSWallGlobalField.cpp';header=prefix+'RANSWallGlobalField.h'
    changed=git('diff','--name-only',a.old,a.new,'--','*.cpp','*.h','*.hpp','*CMakeLists.txt','*.cmake').splitlines()
    if set(changed)!={diag,field,header}:raise ValueError('unexpected native input changes: '+str(changed))
    old={name:git('show',a.old+':'+name) for name in changed};new={name:git('show',a.new+':'+name) for name in changed}
    marker='std::string RANSWallEvolutionPatch::WJMeasurements'
    if old[diag].split(marker)[0]!=new[diag].split(marker)[0]:raise ValueError('governing audit/helper prefix changed')
    start=new[field].index('// This independent native transform deliberately')
    end=new[field].index('RANSWallGlobalField::Vector RANSWallGlobalField::RadialProject',start)
    without=new[field][:start]+new[field][end:]
    if without.strip()!=old[field].strip():raise ValueError('existing global field implementation changed')
    begin=new[header].index('    /** Explicit element ownership')
    finish=new[header].index('    std::array<double,4> NativeEvaluate(',begin)
    if new[header][:begin]+new[header][finish:]!=old[header]:raise ValueError('unexpected public declaration/layout change')
    if old[header].split('private:')[1]!=new[header].split('private:')[1]:raise ValueError('global field data layout declaration changed')
    if 'virtual' in new[header]:raise ValueError('unexpected virtual dispatch')
    h=lambda value:hashlib.sha256(value.encode()).hexdigest()
    result=dict(schema='WJ4BoundedNativeSourceDifferenceV1',producer=git('rev-parse',a.old).strip(),observer=git('rev-parse',a.new).strip(),changed_native_inputs=changed,
        existing_global_field_implementation_unchanged=True,global_field_private_data_unchanged=True,new_method_nonvirtual=True,
        diagnostics_before_readonly_measurements_unchanged=True,
        governing_scope='all other C++/header/CMake inputs identical; differences are an added nonvirtual independent native evaluation method and read-only WJMeasurements diagnostics',
        changed_files={name:dict(old_sha256=h(old[name]),new_sha256=h(new[name])) for name in changed})
    write_once(a.output,result);print(a.output)

if __name__=='__main__':main()
