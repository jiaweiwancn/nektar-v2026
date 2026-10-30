#!/usr/bin/env python3
"""Prepare a continuation XML without changing the model, mesh or timestep."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True


def setting(root, key, value):
    infos = root.find('CONDITIONS/SOLVERINFO')
    node = infos.find("I[@PROPERTY='" + key + "']")
    if node is None:
        node = ET.SubElement(infos, 'I', PROPERTY=key)
    node.set('VALUE', str(value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('session', type=Path, help='Original XML used to produce the field')
    parser.add_argument('field', type=Path, help='Complete .fld or checkpoint file/directory')
    parser.add_argument('output', type=Path, help='New continuation XML')
    parser.add_argument('--steps', type=int, required=True, help='Additional steps, not the original total')
    parser.add_argument('--mode', choices=('History', 'Auto', 'Cold'), default='History')
    parser.add_argument('--history-path', type=Path, help='Relocated complete history directory, when needed')
    args = parser.parse_args()
    if args.steps <= 0:
        parser.error('--steps must be positive')
    session, field, output = (p.resolve() for p in (args.session, args.field, args.output))
    if not session.is_file() or not field.exists():
        parser.error('The original XML and complete field must exist')
    if output in (session, field) or output.suffix.lower() != '.xml':
        parser.error('Use a new .xml output path to retain the original case and field')
    tree = ET.parse(session)
    root = tree.getroot()
    variables = root.findall('CONDITIONS/VARIABLES/V')
    initial = root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
    steps = next((p for p in root.findall('CONDITIONS/PARAMETERS/P')
                  if p.text.split('=', 1)[0].strip() == 'NumSteps'), None)
    if initial is None or not variables or steps is None:
        parser.error('Expected a RANS session with variables, InitialConditions and NumSteps')
    initial.clear()
    initial.set('NAME', 'InitialConditions')
    for variable in variables:
        ET.SubElement(initial, 'F', VAR=variable.text.strip(), FILE=str(field))
    steps.text = 'NumSteps = ' + str(args.steps)
    # The continuation uses a step count even if the original XML included FinalTime.
    for node in list(root.findall('CONDITIONS/PARAMETERS/P')):
        if node.text.split('=', 1)[0].strip() == 'FinalTime':
            root.find('CONDITIONS/PARAMETERS').remove(node)
    setting(root, 'RANSRestartMode', args.mode)
    infos = root.find('CONDITIONS/SOLVERINFO')
    previous = infos.find("I[@PROPERTY='RANSHistoryPath']")
    if previous is not None:
        infos.remove(previous)
    if args.history_path is not None:
        history = args.history_path.resolve()
        if not history.is_dir():
            parser.error('--history-path must be a complete existing directory')
        setting(root, 'RANSHistoryPath', history)
    output.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root, space='  ')
    tree.write(output, encoding='utf-8', xml_declaration=True)
    print(json.dumps(dict(session=str(output), field=str(field), mode=args.mode,
        additional_steps=args.steps, fixture_sha256=hashlib.sha256(output.read_bytes()).hexdigest()), indent=2))


if __name__ == '__main__':
    main()
