#!/usr/bin/env python3
"""Reject the retained native ownership failure, without changing a field."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--packet', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    raw = gzip.decompress(a.packet.read_bytes()) if a.packet.suffix == '.gz' else a.packet.read_bytes()
    packet = json.loads(raw)
    obs = packet['observation']
    assert packet['performed_time_steps'] == 0 and obs['time'] == 0
    assert len(obs['grid']) == len(obs['grid_sides']) == len(obs['grid_geometry_by_field']) == 1
    point = obs['grid'][0]
    containing = [e['geometry'] for e in obs['elements']
                  if all(e['lower'][j] < point[j] < e['upper'][j] for j in (0, 2))]
    assert len(containing) == 1, 'fixture requires a unique strictly containing element'
    selected = obs['grid_geometry_by_field'][0][0]
    ok = selected == containing[0]
    result = dict(status='pass' if ok else 'fail',
                  reason='native_geometry_matches' if ok else 'native_geometry_not_containing_point',
                  point=point, selected_geometry=selected, containing_geometry=containing[0],
                  performed_time_steps=0, packet_sha256=hashlib.sha256(raw).hexdigest(),
                  scope='read-only native ownership check; physical forward validity is independent')
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['reason'])
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
