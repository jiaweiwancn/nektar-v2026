#!/usr/bin/env python3
"""Prepare a neutral-ABL example using the verified RANS study generator."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'Tests' / 'RANS'))
from generate_m6 import generate_abl
from generate_m4 import parameter


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('Use a finite positive value')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('session', type=Path)
    parser.add_argument('--model', choices=('standard', 'standard_original', 'standard_fixed', 'standard_switched',
                                          'cls_original', 'cls_fixed', 'cls_switched'), default='cls_switched')
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=2)
    parser.add_argument('--degree', type=int, default=4)
    parser.add_argument('--layers', type=int, default=12)
    parser.add_argument('--dt', type=positive, default=.5)
    duration = parser.add_mutually_exclusive_group()
    duration.add_argument('--steps', type=int)
    duration.add_argument('--flow-times', type=positive)
    parser.add_argument('--checkpoint', type=int, default=0, help='Checkpoint interval in steps; zero disables intermediate checkpoints')
    parser.add_argument('--backend', choices=('matrixfree', 'native'), default='matrixfree')
    parser.add_argument('--boundary', choices=('outflow', 'periodic'), default='outflow')
    parser.add_argument('--height', type=positive, default=2.)
    parser.add_argument('--roughness', type=float, default=.1,
                        help='AWF equivalent sand-grain height; standard log-law choice uses z0=0.02 m')
    parser.add_argument('--perturbation', type=float, default=0.)
    parser.add_argument('--linear-tolerance', type=positive, default=1e-13)
    args = parser.parse_args()
    if args.steps is not None and args.steps <= 0:
        parser.error('--steps must be positive')
    if args.checkpoint < 0 or not math.isfinite(args.roughness) or args.roughness < 0:
        parser.error('Checkpoint interval and roughness must be nonnegative')
    if not math.isfinite(args.perturbation) or args.perturbation < 0:
        parser.error('--perturbation must be finite and nonnegative')
    session = args.session.resolve()
    if session.suffix.lower() != '.xml':
        parser.error('The generated session must use an .xml output path')
    session.parent.mkdir(parents=True, exist_ok=True)
    try:
        reference = generate_abl(session, model=args.model, dimension=args.dimension,
            degree=args.degree, layers=args.layers, nx=1, ny=1, dt=args.dt,
            steps=args.steps if args.steps is not None else 4, flow_times=args.flow_times,
            checkpoint=args.checkpoint, height=args.height, roughness=args.roughness,
            perturbation=args.perturbation, boundary=args.boundary,
            backend='matrixfree' if args.backend == 'matrixfree' else 'static',
            tolerance=1e-8, monitor=True)
        tree = ET.parse(session)
        parameter(tree.getroot(), 'IterativeSolverTolerance', args.linear_tolerance)
        tree.write(session, encoding='utf-8', xml_declaration=True)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(dict(session=str(session), reference=reference,
        fixture_sha256=hashlib.sha256(session.read_bytes()).hexdigest()), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
