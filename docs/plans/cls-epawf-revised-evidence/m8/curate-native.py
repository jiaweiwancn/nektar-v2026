"""Archive completed M8 native receipts without native fields or histories.

Run on Linux. Full receipts are copied verbatim. Runtime manifests are reduced
only as declared below; the immutable originals and hashes remain authoritative.
An unfinished campaign is refused. Use a fresh output directory each time.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--campaign', action='store_true')
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    index = dict(schema=1, kind='M8 native integration archive', root=str(root),
                 generated_utc=datetime.now(timezone.utc).isoformat(),
                 generator_sha256=digest(Path(__file__).read_bytes()), files=[])

    def archive(name, relative, runtime=False):
        source = root/relative
        raw = source.read_bytes()
        data, transform = raw, 'byte-for-byte copy'
        if runtime:
            full = json.loads(raw)
            omit = ('source_inputs', 'documentation_tools', 'ldd')
            value = {k: v for k, v in full.items() if k not in omit}
            value['curation'] = dict(
                original_path=str(source), original_sha256=digest(raw),
                original_bytes=len(raw), omitted_fields=list(omit),
                source_input_count=len(full.get('source_inputs', {})),
                source_inputs_canonical_sha256=digest(json.dumps(
                    full.get('source_inputs', {}), sort_keys=True,
                    separators=(',', ':')).encode()),
                documentation_tool_count=len(full.get('documentation_tools', {})),
                policy='All runtime images, resolved paths, actual loader bindings, '
                       'tool images, packages, mandatory originals and comparison '
                       'results retained; full input manifests and ldd text on Linux')
            data = (json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False)+'\n').encode()
            transform = 'Omit source_inputs, documentation_tools and raw ldd; add explicit curation metadata'
        (output/name).write_bytes(data)
        index['files'].append(dict(name=name, original_path=str(source),
            original_sha256=digest(raw), original_bytes=len(raw),
            curated_sha256=digest(data), curated_bytes=len(data), transform=transform))

    for model in ('scaled', 'regularized'):
        directory = 'native-completion-001/'+model+'-native/'
        for filename in ('validation.json', 'preregistration.json'):
            archive('native-'+model+'-'+filename, directory+filename)
    archive('native-completion-runner.json', 'native-completion-001/receipt.json')
    for period in ('before', 'after'):
        archive('native-completion-runtime-'+period+'.json',
                'native-completion-001/runtime-'+period+'.json', True)
    for filename in ('validation.json', 'preregistration.json'):
        archive('mms-completion-'+filename, 'mms-completion-001/'+filename)
        archive('flow-'+filename, 'flow-001/qualification/'+filename)
    archive('flow-runner.json', 'flow-001/receipt.json')
    for period in ('before', 'after'):
        archive('flow-runtime-'+period+'.json', 'flow-001/runtime-'+period+'.json', True)
    for number in ('001', '002'):
        archive('initialization-build-'+number+'.json',
                'initialization/build-'+number+'/receipt.json')
        archive('resource-adjustment-'+number+'.json', 'resource-adjustment-'+number+'.json')
    archive('initialization-build-runtime.json',
            'initialization/build-002/runtime-fingerprint.json', True)
    for number in ('001', '002', '003'):
        archive('pilot-'+number+'-runner.json', 'pilot-'+number+'/receipt.json')
        if (root/('pilot-'+number+'/qualification/validation.json')).is_file():
            archive('pilot-'+number+'-validation.json', 'pilot-'+number+'/qualification/validation.json')
    archive('mms-completion-runner.json', 'mms-completion-001-runner/receipt.json')
    archive('original-controls-attempt-001.json', 'initialization/original-controls-001/failure-receipt.json')
    archive('run-original-controls-attempt-001.py', 'run-original-controls-001.py')
    archive('original-controls-attempt-002.json', 'initialization/original-controls-002/receipt.json')
    archive('original-controls-validation.json', 'initialization/original-controls-003/receipt.json')
    archive('original-controls-ctest.xml', 'initialization/original-controls-003/ctest.xml')
    for period in ('before', 'after'):
        archive('original-controls-runtime-'+period+'.json',
                'initialization/original-controls-003/runtime-'+period+'.json', True)
    if args.campaign:
        receipt = json.loads((root/'campaign-001/receipt.json').read_text())
        if receipt['status'] in ('running', 'pending'):
            raise RuntimeError('The campaign has not completed')
        archive('campaign-runner.json', 'campaign-001/receipt.json')
        for period in ('before', 'after'):
            archive('campaign-runtime-'+period+'.json', 'campaign-001/runtime-'+period+'.json', True)
        archive('legacy-ctest.xml', 'campaign-001/ctest.xml')
        archive('legacy-regression.log', 'campaign-001/legacy-regression.log')
    (output/'native-curation-index.json').write_text(
        json.dumps(index, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    print(json.dumps(dict(files=len(index['files']), bytes=sum(
        f['curated_bytes'] for f in index['files']), output=str(output))))


if __name__ == '__main__':
    main()
