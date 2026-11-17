#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Code-generation controls; not native pressure-space qualification."""
import argparse, importlib.util, shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from cls_epawf_m8_r1_cases import parameters, solver_info
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha
from wp2_neutral_cases import case


def main():
    p = argparse.ArgumentParser()
    for name in ('baseline', 'h1', 'preflight', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(exist_ok=False)
    spec = importlib.util.spec_from_file_location('retained_wp2_factory', a.baseline)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    pre = load(a.preflight)
    records = {}
    for model in ('CLS', 'standard'):
        folder = a.output/model
        folder.mkdir()
        common = dict(level=4, order=16, dt=.01, projection='RadialEnergyV1')
        xml = old.case(folder/'baseline', a.h1, pre, model, **common)
        old_xml, old_definition = xml.read_bytes(), (xml.parent/'wp2-definition.json').read_bytes()
        shutil.move(xml.parent, folder/'retained-default')
        default = case(folder/'baseline', a.h1, pre, model, **common)
        require(default.read_bytes() == old_xml and (default.parent/'wp2-definition.json').read_bytes() == old_definition,
                'default WP2 XML/definition changed')
        enriched = case(folder/'inner-p2', a.h1, pre, model, pressure_space='inner-p2-v1', **common)
        left, right = ET.parse(default).getroot(), ET.parse(enriched).getroot()
        for tag in ('VERTEX', 'EDGE', 'ELEMENT'):
            require(ET.tostring(left.find('./GEOMETRY/'+tag)) == ET.tostring(right.find('./GEOMETRY/'+tag)),
                    'pressure control changed physical mesh '+tag)
        require(parameters(left) == parameters(right), 'pressure control changed a physical/temporal parameter')
        li, ri = solver_info(left), solver_info(right)
        for key in ('RANSWALLCHECKPOINT', 'RANSWALLRESTART', 'RANSWALLMODELMANIFESTFILE'):
            li.pop(key, None);ri.pop(key, None)
        require(li == ri, 'pressure control changed another solver policy')
        require((default.parent/'model-manifest.json').read_bytes() == (enriched.parent/'model-manifest.json').read_bytes(),
                'pressure control changed model calibration or physical boundary identity')
        pressure = [(e.get('COMPOSITE'), e.get('NUMMODES')) for e in right.findall('./EXPANSIONS/E[@FIELDS="p"]')]
        require(pressure == [('C[1]', '3,3'), ('C[2]', '2,2')] and right.find('./GEOMETRY/DOMAIN').text == 'C[1,2]',
                'invalid regional expansion/domain')
        definition = load(enriched.parent/'wp2-definition.json')
        require(definition['estimated_slots'] == 9093 and definition['native_dimensions_required'],
                'wrong inner-P2 size estimate or missing native confirmation requirement')
        fine = case(folder/'level5-inner-p2', a.h1, pre, model, level=5, order=16, dt=.01,
                    projection='RadialEnergyV1', pressure_space='inner-p2-v1')
        require(load(fine.parent/'wp2-definition.json')['estimated_slots'] == 11609, 'wrong confirmation estimate')
        records[model] = dict(default_XML_and_definition_bitwise=True, physical_mesh_and_parameters_unchanged=True,
            model_manifest_unchanged=True, pressure_expansions=pressure,
            estimated_level4_slots=9093, estimated_level5_slots=11609, native_execution='UNEXECUTED')
    for kwargs in (dict(level=3, pressure_space='inner-p2-v1'), dict(level=4, pressure_space='all-p2-v1')):
        try:
            case(a.output/'invalid', a.h1, pre, 'standard', **kwargs)
        except ValueError as error:
            require('pressure representation' in str(error), 'unexpected scope rejection')
        else:
            raise AssertionError('undeclared pressure space accepted')
    atomic(a.output/'result.json', dict(version='WP2PressureGenerationChecksV1', passed=True, models=records,
        baseline_sha256=sha(a.baseline), source_sha256=sha(Path(__file__)),
        case_factory_sha256=sha(Path(__file__).with_name('wp2_neutral_cases.py')), governing_solves=0,
        scope='Native compilation, actual dimensions, rank/pressure action and evolution remain separate gates.'))


if __name__ == '__main__':
    main()
