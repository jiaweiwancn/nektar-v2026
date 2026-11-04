#!/usr/bin/env python3
"""Re-extract the same retained R2 records, preserving native signed zeros."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import struct


H0 = 'b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71'
REFERENCE = '4ad26ed0167bb1765bd0be2c150928280b69dbd82ace8a27ea5f9c8c7f4d0409'


def signed_int(token):
    return -0.0 if token == '-0' else int(token)


def load(text):
    return json.loads(text, parse_int=signed_int,
        parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


def pin(path):
    path = Path(path).resolve(); digest = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            digest.update(part)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


@lru_cache(maxsize=2)
def record(path, index, expected):
    if pin(path)['sha256'] != expected:
        raise ValueError('Historical source record changed')
    with Path(path).open() as stream:
        for i, line in enumerate(stream):
            if i == index:
                return load(line)
    raise ValueError('Missing complete native record')


def packed(values):
    if len(values) != 2165:
        raise ValueError('Unexpected outer6 vector length')
    return struct.pack('<2165d', *values)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior-ledger', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); old = load(a.prior_ledger.read_text())
    a.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema=2, scope='Same archived R2 records; binary64 signed zero preserved',
        prior_ledger=pin(a.prior_ledger), extractor=pin(__file__),
        historical_native_source=old['historical_native_source'], native_solver=old['native_solver'],
        inputs=old['inputs'], states=[], native_calls=0, scientific_reconstructions=0)
    references, starts = [], []
    for row in old['states'][:2]:
        ref = row['preflight_report']
        if pin(ref['path']) != ref:
            raise ValueError('Historical native initialization report changed')
        report = load(Path(ref['path']).read_text())
        references.append(report['initial']); starts.append(report['state'])
    initial = packed(starts[1])
    if (hashlib.sha256(initial).hexdigest() != H0 or packed(starts[0]) != initial or
        starts[1] != references[1] or
        any(hashlib.sha256(packed(ref)).hexdigest() != REFERENCE for ref in references)):
        raise ValueError('Original outer6 native initialization binding differs')
    names = ('linear-auxiliary', 'saturatedr2-auxiliary', 'boundary-path-2', 'boundary-path-3', 'frozen-terminal')
    for name in names:
        previous = a.prior_ledger.parent/(name+'.json')
        packet = load(previous.read_text()); origin = packet['source_record']; trace = origin['trace']
        native = record(trace['path'], origin['record_index'], trace['sha256'])
        if native['alpha'] != packet['alpha'] or native['time'] != packet['time']:
            raise ValueError('Historical stage scope differs')
        internal = native['iterations'][origin['iteration_index']] if 'iteration_index' in origin else native
        current = dict(state=internal['state'], history=native['history'], residual=internal['residual'])
        changes = {}
        for field, values in current.items():
            if values != packet[field]:
                raise ValueError('A saved coefficient changed numerically: '+name+':'+field)
            before, after = packed(packet[field]), packed(values)
            changes[field] = sum(x != y for x, y in zip(struct.iter_unpack('<Q', before), struct.iter_unpack('<Q', after)))
            packet[field] = values
            vector = a.output/(name+'-'+field+'.bin'); vector.write_bytes(after)
            packet[field+'_binary'] = pin(vector)
        if packet['history'] != starts[1]:
            raise ValueError('An archived stage history differs in value from h0')
        # The resolution probe rebuilds the frozen reference, before the
        # zero-amplitude initialization expression changes zero signs.
        packet['require_history_h0'] = packet['history_binary']['sha256'] == REFERENCE
        packet['history_equals_original_h0_coefficient_values'] = True
        path = a.output/(name+'.json'); path.write_text(json.dumps(packet, separators=(',', ':'), allow_nan=False)+'\n')
        result['states'].append(dict(name=name, packet=pin(path), previous_packet=pin(previous),
            sign_bit_changes_only=changes, source_record=origin,
            vectors={field:packet[field+'_binary'] for field in current}))
    h0 = a.output/'h0.bin'; h0.write_bytes(initial)
    result['initial_h0'] = pin(h0)
    result['frozen_reference_sha256'] = REFERENCE
    result['initial_outer_k_direction'] = old['initial_outer_k_direction']
    result['success'] = True
    (a.output/'input-ledger.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(states=len(result['states']), initial_h0=result['initial_h0'])))


if __name__ == '__main__':
    main()
