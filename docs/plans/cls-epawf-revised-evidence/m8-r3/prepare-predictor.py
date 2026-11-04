#!/usr/bin/env python3
"""Prepare a hash-bound R3 numerical predictor; metadata/binary checks only."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import xml.etree.ElementTree as ET

H0 = '4ad26ed0167bb1765bd0be2c150928280b69dbd82ace8a27ea5f9c8c7f4d0409'
NUMERICAL_SELECTORS = {'RANSWALLSTAGETEST', 'RANSWALLNONLINEARPOLICY', 'RANSWALLACTIONPOLICY',
    'RANSWALLR3PROGRESSFILE', 'RANSWALLR3PREDICTORFILE', 'RANSWALLFEASIBLEFILE',
    'RANSWALLSOLVEDIAGNOSTICS', 'RANSWALLSOLVEDIAGNOSTICREPLAY', 'RANSWALLSOLVEDIAGNOSTICFILE'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def identity(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def vector(values):
    require(isinstance(values, list) and all(type(x) in (int, float) and math.isfinite(x) for x in values),
        'Nonfinite or malformed coefficient vector')
    return hashlib.sha256(struct.pack('<'+str(len(values))+'d', *values)).hexdigest()


def load(path):
    # Native iostream output represents negative zero as the integer token -0.
    return json.loads(Path(path).read_text(), parse_int=lambda x: -0.0 if x == '-0' else int(x),
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def protected_xml(path):
    root = ET.parse(path).getroot()
    for parent in root.iter():
        for child in list(parent):
            if child.tag == 'I' and child.get('PROPERTY', '').upper() in NUMERICAL_SELECTORS:
                parent.remove(child)
            elif child.tag == 'P' and (child.text or '').split('=', 1)[0].strip().upper() == 'RANSWALLMAXITERATIONS':
                parent.remove(child)
    value = ET.canonicalize(ET.tostring(root, encoding='unicode'), strip_text=True)
    return hashlib.sha256(value.encode()).hexdigest()


def prepare(args):
    packet_pin, ledger_pin = identity(args.packet), identity(args.input_ledger)
    packet, current, runtime = (load(path) for path in (args.packet, args.current_preflight, args.runtime_pin))
    require(packet.get('schema') == 1 and packet.get('inner_stress_model') == 'SaturatedR2' and
        packet.get('time') == .0025 and packet.get('alpha') == 3200 and
        packet.get('source_record', {}).get('converged') is False, 'Predictor is not the saved failed Sat auxiliary')
    require(len(packet['state']) == len(packet['history']) == len(packet['residual']) == 2165,
        'Predictor changed the fixed outer6 coefficient shape')
    pins = {}
    for field in ('state', 'history', 'residual'):
        pin = packet[field+'_binary']
        require(identity(pin['path']) == pin and pin['bytes'] == 2165*8 and vector(packet[field]) == pin['sha256'],
            'Historical packet/binary binding changed: '+field)
        pins[field] = pin
    require(vector(packet['history']) == H0 and vector(current['initial']) == H0 and
        current.get('time') == current.get('steps') == current.get('eta') == 0 and current.get('dt') == .0025 and
        current.get('inner_stress_model') == 'SaturatedR2' and current.get('nonlinear_policy') == 'ScalarFeasibleR3',
        'Current initialization changed original h0, target scope or selected policy')
    origin_trace = packet['source_record']['trace']
    require(identity(origin_trace['path']) == origin_trace, 'Historical native trace changed')
    # The root's immutable extraction ledger binds the origin packet. The
    # large historical JSONL is hashed, never reparsed as a new trajectory.
    ledger = load(args.input_ledger)
    def contains_pin(value):
        if isinstance(value, dict):
            if value.get('path') == packet_pin['path'] and value.get('sha256') == packet_pin['sha256']:
                return True
            return any(contains_pin(item) for item in value.values())
        return isinstance(value, list) and any(contains_pin(item) for item in value)
    require(contains_pin(ledger), 'Immutable input ledger does not bind this exact historical packet')
    protected = [protected_xml(path) for path in (args.reference_session, args.current_preflight_session, args.target_session)]
    require(len(set(protected)) == 1, 'R3 target changed protected physics, geometry, expansions or initial functions')
    target = ET.parse(args.target_session).getroot()
    selectors = target.findall('.//SOLVERINFO/I')
    info = {row.get('PROPERTY', '').upper():row.get('VALUE') for row in selectors}
    require(len(info) == len(selectors), 'Repeated target solver-info property')
    def parameters(path):
        values = {}
        for row in ET.parse(path).getroot().findall('.//PARAMETERS/P'):
            name, value = (part.strip() for part in (row.text or '').split('=', 1))
            require(name.upper() not in values, 'Repeated numerical parameter')
            values[name.upper()] = value
        return values
    current_parameters, target_parameters = (parameters(path) for path in
        (args.current_preflight_session, args.target_session))
    require(current_parameters.get('RANSWALLMAXITERATIONS') == target_parameters.get('RANSWALLMAXITERATIONS'),
        'Current native signature was initialized with a different Newton budget')
    require(info.get('RANSWALLSTAGETEST') in ('ExactTargetR3', 'JacobianBenchmarkR3') and
        info.get('RANSWALLNONLINEARPOLICY') == 'ScalarFeasibleR3' and
        info.get('RANSWALLACTIONPOLICY') == 'VerifiedR3Action' and info.get('RANSWALLINNERSTRESS') == 'SaturatedR2',
        'Target is not an explicitly selected R3 Sat stage/benchmark')
    require(not any(row.get('NAME', '').upper() in ('BODYFORCEINNER', 'BODYFORCEOUTER', 'SCALARSOURCEINNER', 'SCALARSOURCEOUTER')
        for row in target.findall('.//FUNCTION')), 'R3 target introduced prescribed forcing')
    require(runtime.get('success') is True and runtime.get('clean') is True and
        runtime.get('source_commit') == args.expected_source, 'Unsuccessful or changed current runtime source pin')
    key = 'build/IncNavierStokesSolver'
    executable = identity(runtime['runtime_resolved_paths'][key])
    require(executable['sha256'] == runtime['runtime_artifacts'][key], 'Current executable changed')
    return dict(schema=1, role='nonconverged_auxiliary_predictor_only',
        encoding='little-endian IEEE754 binary64', state=pins['state'], history=pins['history'],
        expected_h0_sha256=H0, current_executable=executable, current_target_xml=identity(args.target_session),
        current_native_signature_sha256=hashlib.sha256(current['signature'].encode()).hexdigest(),
        current_preflight=identity(args.current_preflight), current_preflight_session=identity(args.current_preflight_session),
        current_runtime=identity(args.runtime_pin), producing_source_commit=args.expected_source,
        protected_physical_xml_sha256=protected[0], historical_input_ledger=ledger_pin, historical_packet=packet_pin,
        origin=dict(alpha=3200., time=.0025, converged=False, inner_stress_model='SaturatedR2',
            source_record=packet['source_record'], residual=pins['residual'],
            scaled_block_norms=packet['scaled_block_norms'], historical_failure_is_not_reclassified=True),
        target=dict(alpha=400., time=.0025, eta=0., tolerance=2e-12, history_advanced=False),
        predictor_is_not_a_physical_time_sample=True, coefficients_are_not_lifted_or_clipped=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('packet', 'input-ledger', 'reference-session', 'current-preflight',
                 'current-preflight-session', 'target-session', 'runtime-pin', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--expected-source', required=True)
    args = p.parse_args()
    require(not args.output.exists(), 'Use a fresh predictor manifest path')
    result = prepare(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(manifest=identity(args.output), native_calls=0, original_h0=H0)))


if __name__ == '__main__':
    main()
