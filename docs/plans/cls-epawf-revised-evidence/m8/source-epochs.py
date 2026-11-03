"""Bind stored M8 runs to immutable Git source epochs (no numerical work)."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo = args.repository.resolve()
    def git(*argv):
        return subprocess.check_output(['git', '-C', str(repo), *argv])
    def sha(raw):
        return hashlib.sha256(raw).hexdigest()
    labels = ('22d5c930', '0736319', 'aa854f4', '11aa94c', '2a93b04', '6284787', 'ed34052')
    commits = {label: git('rev-parse', label).decode().strip() for label in labels}
    trees = {}
    for label, commit in commits.items():
        entries = {}
        for record in git('ls-tree', '-r', '-z', commit).split(b'\0'):
            if not record:
                continue
            metadata, name = record.split(b'\t', 1)
            mode, kind, object_id = metadata.decode().split()
            if kind == 'blob':
                entries[name.decode()] = object_id
        trees[label] = entries
    suffixes = {'.cpp', '.h', '.hpp', '.c', '.cc', '.cxx', '.hxx', '.inl'}
    def cpp_delta(before, after):
        paths = sorted(name for name in set(trees[before]) | set(trees[after])
                       if name.startswith(('library/', 'solvers/')) and
                       Path(name).suffix in suffixes)
        return [name for name in paths if trees[before].get(name) != trees[after].get(name)]
    expected = {
        ('aa854f4', '6284787'): ['library/UnitTests/SolverUtils/RANSCurvedWallProbe.cpp'],
        ('6284787', 'ed34052'): ['solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp'],
    }
    deltas = []
    for (before, after), paths in expected.items():
        actual = cpp_delta(before, after)
        assert actual == paths, (before, after, actual)
        deltas.append(dict(before=commits[before], after=commits[after],
                           changed_cpp_or_header_paths=actual))
    kernel = ['library/SolverUtils/Turbulence/RANSRobustVariant.cpp',
              'library/SolverUtils/Turbulence/RANSRobustVariant.h',
              'library/UnitTests/SolverUtils/RANSRobustVariantProbe.cpp',
              'library/UnitTests/SolverUtils/RANSTransverseEvolutionProbe.cpp']
    kernel_records = {}
    for path in kernel:
        assert all(trees[label][path] == trees['0736319'][path] for label in labels[1:]), path
        raw = git('show', commits['ed34052']+':'+path)
        kernel_records[path] = dict(sha256=sha(raw), git_blob_id=trees['ed34052'][path],
                                    unchanged_since=commits['0736319'])
    mandatory = [
        'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex',
        'docs/plans/cls-epawf-revised-contract-v1.json',
        'library/SolverUtils/Turbulence/RANSModels.cpp',
        'library/SolverUtils/Turbulence/RANSModels.h',
        'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp',
        'library/SolverUtils/Turbulence/AnalyticalWallFunction.h',
        'library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp',
        'library/SolverUtils/Turbulence/EquilibriumAWFProfiles.h',
        'library/SolverUtils/Turbulence/EquilibriumAWFTransition.cpp',
        'library/SolverUtils/Turbulence/EquilibriumAWFTransition.h',
        'library/LibUtilities/Interpreter/Interpreter.cpp',
    ]
    originals = {}
    for path in mandatory:
        assert all(trees[label][path] == trees[labels[0]][path] for label in labels), path
        originals[path] = dict(sha256=sha(git('show', commits['ed34052']+':'+path)),
                               git_blob_id=trees['ed34052'][path])
    result = dict(schema=1, kind='Static source-epoch binding; no numerical execution',
                  generated_utc=datetime.now(timezone.utc).isoformat(),
                  generator_sha256=sha(Path(__file__).read_bytes()),
                  git_object_hash_algorithm='SHA-1 blob IDs; explicit content SHA256 also supplied',
                  commits=commits, cpp_deltas=deltas, kernel_inputs=kernel_records,
                  frozen_inputs_unchanged_from_M7=originals,
                  limitations='This verifies source identity, not compiled execution. '
                  'Use the producing runtime fingerprints for binary/loader/environment identity. '
                  'Python proof/audit changes keep their separate producing commits.')
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(dict(success=True, epochs=len(commits), frozen_files=len(originals),
                          kernel_files=len(kernel_records), cpp_deltas=len(deltas))))


if __name__ == '__main__':
    main()
