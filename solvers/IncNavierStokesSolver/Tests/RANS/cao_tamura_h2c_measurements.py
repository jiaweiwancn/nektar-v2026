#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only frozen RR stress/RMS height inventory. No calibration from hill data."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def wall(x):return .04*math.cos(math.pi*x/.2)**2 if abs(x)<.1 else 0.
def distance(x,z):
    # Independent piecewise minimization; include both feet and the flat floor.
    def f(s):return (x-s)**2+(z-wall(s))**2
    candidates=[f(-.1),f(.1)]
    if x<=-.1 or x>=.1:candidates.append(z*z)
    for i in range(100):
        a=-.1+.2*i/100;b=a+.002
        for _ in range(48):
            c=a+(b-a)*.3819660112501051;d=a+(b-a)*.6180339887498949
            if f(c)<f(d):b=d
            else:a=c
        candidates.append(f((a+b)/2))
    return math.sqrt(min(candidates))
def main():
    a=argparse.ArgumentParser();a.add_argument('--data',type=Path,required=True);a.add_argument('--output',type=Path,required=True);q=a.parse_args()
    snapshot=json.loads((q.data/'snapshot.json').read_text())
    # The existing data auditor is authoritative; independently retain the
    # entire snapshot and content hashes used by this inventory.
    series={s['id']:s for s in json.loads((q.data/'series.json').read_text())}
    observations=[]
    for r in csv.DictReader((q.data/'profiles.csv').open()):
        s=series[r['series_id']]
        if s['case'] not in ('R','RR') or not any(k in s['quantity'] for k in ('sigma','uv')):continue
        x=.04*float(r['x_over_H']) if r['x_over_H'] else None
        z=.04*float(r['y_over_H']) if r['y_over_H'] else None
        if s['height_kind']=='probe_fixed_height':
            z=wall(x)+.007 # frozen Fig12 RR vertical-above-floor convention
        if s['panel']=='f11':z=None # downstream delta_7.5 unknown, never use delta_in
        d=(z if s['case']=='R' else distance(x,z)) if x is not None and z is not None else None
        observations.append(dict(series=s['id'],figure=s['panel'],point_id=int(r['point_id']),quantity=s['quantity'],
            geometry='rough-flat' if s['case']=='R' else 'rough-hill',x_m=x,z_global_m=z,nearest_floor_distance_m=d,height_convention=s['height_kind'],
            published_ordinate=float(r['axis_y']),axis_y_extraction_bound=float(r['axis_y_bound']),
            conversion=s['coordinate_conversion'],uncertainty='frozen extraction envelope; nominal scales/displacement; not measurement sigma',
            unresolved='delta_7.5 not supplied' if z is None else None,
            quality=r['quality'],weight=float(r['weight']),role=r['role'],
            displacement_already_applied=True,strip_membership={str(H):None if d is None else ('inside' if d<=H else 'outside') for H in (.0055,.006,.0065)}))
    out=dict(version='CaoTamuraH2cHeightInventoryV1',data_version='cao-tamura-2006-v1.1',
        snapshot=snapshot,files={f.name:digest(f) for f in q.data.iterdir() if f.is_file()},
        known_minimum_RR_fixed_traverse_m=.007,fig11_physical_height='unresolved: local downstream delta_7.5 absent',
        hill_data_use='measurement-clearance inventory only, never fit',observations=observations,
        screens=dict(H_plus_minimum=100,molecular_fraction_maximum=.05,status='provisional project screens, no universal high-Re or roughness-sublayer proof'),
        policy=dict(exact_analytical_provider_relative_error=1e-12,roughness_root_absolute_velocity_error=1e-11,
            exact_flat_preservation_normalized_absolute=1e-4,exact_flat_selected_spatial_block_target=1e-3,
            initial_traction_reference_error_target=.005,initial_gradient_relative_L2_target=.05,
            numerical_comparison_floor=1e-12,scope='representation targets assessed as results, never replace hard acceptance; no response-relative drift criterion for equilibrium',
            hill_dt=2e-6,hill_steps=8,hill_fine_dt=1e-6,hill_fine_steps=16,
            cost_cap_seconds_per_native_leg=1800,no_long_continuation=True))
    q.output.parent.mkdir(parents=True,exist_ok=True);q.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
