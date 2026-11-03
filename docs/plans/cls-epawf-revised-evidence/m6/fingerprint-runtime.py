#!/usr/bin/env python3
"""Pin M6 numerical inputs and every resolved Linux build runtime image.

ROOT contains the isolated source/ and build/ directories. A later comparison
admits a documentation-only closure commit only when the complete selected
source sets, runtime image sets, dependency bindings and numerical environment
remain identical. There is no exemption for changed probes or libraries.
Documentation postprocessors are recorded separately from numerical inputs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import sys
import uuid

sys.dont_write_bytecode = True

SUFFIXES = {'.cpp', '.h', '.hpp', '.c', '.cc', '.cxx', '.hxx', '.inl', '.py', '.xml', '.tst', '.cmake'}
REQUIRED = (
    'docs/plans/cls-epawf-revised-contract-v1.json',
    'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex',
    'library/SolverUtils/Turbulence/RANSModels.cpp',
    'library/SolverUtils/Turbulence/RANSModels.h',
    'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp',
    'library/SolverUtils/Turbulence/AnalyticalWallFunction.h',
    'library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp',
    'library/SolverUtils/Turbulence/EquilibriumAWFProfiles.h',
    'library/SolverUtils/Turbulence/EquilibriumAWFTransition.cpp',
    'library/SolverUtils/Turbulence/EquilibriumAWFTransition.h',
    'library/LibUtilities/Interpreter/Interpreter.cpp',
)
ENVIRONMENT = (
    'LD_LIBRARY_PATH', 'LD_PRELOAD', 'LD_AUDIT',
    'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_PROC_BIND',
    'OMP_PLACES', 'OMP_SCHEDULE', 'OMP_MAX_ACTIVE_LEVELS', 'MKL_DYNAMIC',
    'MKL_CBWR', 'OPENBLAS_CORETYPE', 'PYTHONPATH', 'PYTHONHOME', 'PYTHONHASHSEED',
    'PYTHONDONTWRITEBYTECODE', 'LANG', 'LC_ALL',
)
NUMERICAL_PACKAGES = ('numpy', 'scipy', 'mpmath', 'sympy')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def stamp(path):
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def digest(path):
    path = Path(path)
    before = stamp(path)
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    require(stamp(path) == before, 'File changed while being hashed: '+str(path))
    return result.hexdigest()


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *arguments],
                                   text=True)


def elf_type(path):
    """Return ELF e_type, avoiding filename/execute-bit assumptions."""
    try:
        with path.open('rb') as stream:
            header = stream.read(18)
    except (OSError, PermissionError):
        return None
    if len(header) < 18 or header[:4] != b'\x7fELF' or header[5] not in (1, 2):
        return None
    return struct.unpack(('<' if header[5] == 1 else '>')+'H', header[16:18])[0]


def build_images(build):
    images = {}
    # Include all executable and shared-object ELF images, even when a shared
    # object lacks its execute bit. Object files and archives are link inputs,
    # not runtime images; their final linked images are hashed below.
    for directory, _, files in os.walk(build, followlinks=False):
        for name in files:
            path = Path(directory)/name
            if path.is_file() and elf_type(path) in (2, 3):
                images['build/'+path.relative_to(build).as_posix()] = path
    return dict(sorted(images.items()))


def dynamic_dependencies(path, environment):
    result = subprocess.run(['ldd', str(path)], env=environment, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=30, check=False)
    text = result.stdout
    require('not found' not in text, 'Unresolved dependency for '+str(path)+': '+text)
    static = 'statically linked' in text or 'not a dynamic executable' in text
    require(result.returncode == 0 or static, 'ldd failed for '+str(path)+': '+text)
    bindings, virtual = {}, []
    for line in text.splitlines():
        line = line.strip()
        if '=>' in line:
            soname, destination = (part.strip() for part in line.split('=>', 1))
            destination = re.sub(r'\s+\(0x[0-9a-fA-F]+\)\s*$', '', destination)
            if destination.startswith('/'):
                dependency = Path(destination)
                require(dependency.is_file(), 'Missing resolved library: '+destination)
                bindings[soname] = str(dependency)
        elif line.startswith('/'):
            destination = re.sub(r'\s+\(0x[0-9a-fA-F]+\)\s*$', '', line)
            dependency = Path(destination)
            require(dependency.is_file(), 'Missing ELF interpreter: '+destination)
            bindings['ELF interpreter'] = str(dependency)
        elif line.startswith('linux-vdso') or line.startswith('linux-gate'):
            virtual.append(line.split()[0])
    return dict(bindings=dict(sorted(bindings.items())), virtual=sorted(virtual),
                static=static, exit_code=result.returncode, raw=text)


def package_identity(name):
    distribution = importlib.metadata.distribution(name)
    entries, total_bytes = [], 0
    files=distribution.files
    if files is None:
        # Debian/Ubuntu .egg-info packages need not ship a RECORD file.
        # Hash their actual installed package trees instead of treating an
        # absent packaging manifest as absent numerical code.
        module=importlib.import_module(name)
        roots=[Path(value) for value in getattr(module,'__path__',[])]
        if not roots and getattr(module,'__file__',None):roots=[Path(module.__file__).parent]
        located={str(path):path for root in roots for path in root.rglob('*') if path.is_file()}
        files=sorted(located)
    for item in sorted(files or [], key=str):
        name_on_disk = str(item)
        if not (name_on_disk.endswith('.py') or '.so' in Path(name_on_disk).name):
            continue
        path = Path(distribution.locate_file(item))
        if path.is_file():
            total_bytes += path.stat().st_size
            entries.append([name_on_disk, digest(path)])
    require(bool(entries), 'No installed numerical package files found: '+name)
    encoded = json.dumps(entries, separators=(',', ':'), ensure_ascii=True).encode()
    return dict(version=distribution.version, files=len(entries), bytes=total_bytes,
                installed_code_sha256=hashlib.sha256(encoded).hexdigest())


def compare(previous, current):
    kinds = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
             'runtime_environment', 'python_packages', 'tool_images')
    differences = {}
    for kind in kinds:
        before, after = previous.get(kind), current.get(kind)
        require(isinstance(before, dict) and isinstance(after, dict),
                'Previous fingerprint lacks comparison field '+kind)
        added, removed = sorted(set(after)-set(before)), sorted(set(before)-set(after))
        changed = sorted(key for key in set(before)&set(after) if before[key] != after[key])
        if added or removed or changed:
            differences[kind] = dict(added=added, removed=removed, changed=changed)
    return differences


def capture(root):
    require(platform.system() == 'Linux', 'ELF runtime qualification belongs on Linux')
    source, build = root/'source', root/'build'
    require(source.is_dir() and build.is_dir(), 'ROOT must contain source/ and build/')
    status = git(source, 'status', '--porcelain=v1', '--untracked-files=all')
    require(not status.strip(), 'Source checkout is not clean: '+status)
    commit = git(source, 'rev-parse', 'HEAD').strip()
    tracked = git(source, 'ls-files', '-z').split('\0')
    selected = [name for name in tracked if name and not name.startswith('docs/') and
                (Path(name).suffix.lower() in SUFFIXES or Path(name).name == 'CMakeLists.txt')]
    for name in REQUIRED:
        require(name in tracked and (source/name).is_file(), 'Missing tracked frozen/numerical input '+name)
    source_inputs = {name: digest(source/name) for name in sorted(set(selected)|set(REQUIRED))}
    documentation_tools = {name: digest(source/name) for name in sorted(tracked)
                           if name.startswith('docs/') and Path(name).suffix.lower() == '.py'}
    contract_path = 'docs/plans/cls-epawf-revised-contract-v1.json'
    contract = json.loads((source/contract_path).read_text())
    require(source_inputs[contract['note']['path']] == contract['note']['sha256'],
            'Revised note does not match the frozen model contract')
    images = build_images(build)
    require(any(path.name == 'IncNavierStokesSolver' for path in images.values()), 'Missing built evolution solver')
    require(any(path.name == 'RANSEquilibriumTransitionProbe' for path in images.values()), 'Missing built transition probe')
    environment = dict(os.environ)
    # Match the actual native qualification runner rather than relying on an
    # interactive shell's BLAS/OpenMP defaults.
    environment.update({name: '1' for name in
                        ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')})
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    overlay = str(build/'geometry-lib')
    ld_paths = [value for value in environment.get('LD_LIBRARY_PATH', '').split(':') if value]
    environment['LD_LIBRARY_PATH'] = ':'.join([overlay]+[value for value in ld_paths if value != overlay])
    runtime_environment = {key: environment.get(key) for key in ENVIRONMENT}
    runtime_environment.update({key: value for key, value in environment.items()
                                if key.startswith(('OMPI_MCA_', 'PMIX_MCA_', 'UCX_', 'I_MPI_', 'MPICH_'))})
    artifacts, bindings, ldd_output, snapshots, resolved = {}, {}, {}, {}, {}
    for key, path in images.items():
        artifacts[key], snapshots[path] = digest(path), stamp(path)
        resolved[key] = str(path.resolve())
        dependencies = dynamic_dependencies(path, environment)
        bindings[key] = {name: dependencies[name] for name in ('bindings', 'virtual', 'static')}
        ldd_output[key] = dependencies['raw']
        for destination in dependencies['bindings'].values():
            library = Path(destination)
            library_key = 'library/'+destination
            if library_key not in artifacts:
                artifacts[library_key], snapshots[library] = digest(library), stamp(library)
                resolved[library_key] = str(library.resolve())
    tools = {}
    for tool in ('python3', 'mpiexec', 'ctest', 'cmake'):
        path = shutil.which(tool)
        require(path is not None, 'Missing qualification tool '+tool)
        target = Path(path)
        tools[tool] = dict(path=path, resolved=str(target.resolve()), sha256=digest(target))
        snapshots[target] = stamp(target)
        if elf_type(target) in (2, 3):
            dependencies = dynamic_dependencies(target, environment)
            bindings['tool/'+tool] = {name: dependencies[name] for name in ('bindings', 'virtual', 'static')}
            ldd_output['tool/'+tool] = dependencies['raw']
            for destination in dependencies['bindings'].values():
                library, key = Path(destination), 'library/'+destination
                if key not in artifacts:
                    artifacts[key], snapshots[library] = digest(library), stamp(library)
                    resolved[key] = str(library.resolve())
    actual_python = Path(sys.executable)
    tools['fingerprint_interpreter'] = dict(path=str(actual_python), resolved=str(actual_python.resolve()),
                                           sha256=digest(actual_python))
    snapshots[actual_python] = stamp(actual_python)
    packages = {name: package_identity(name) for name in NUMERICAL_PACKAGES}
    require(set(build_images(build)) == set(images), 'Build runtime image set changed during fingerprint')
    for path, before in snapshots.items():
        require(stamp(path) == before, 'Runtime artifact changed during fingerprint: '+str(path))
    require(git(source, 'rev-parse', 'HEAD').strip() == commit, 'Source revision changed during fingerprint')
    require(not git(source, 'status', '--porcelain=v1', '--untracked-files=all').strip(),
            'Unexpected source work appeared during fingerprint')
    return dict(schema=2, milestone='M6', source_root=str(source), build_root=str(build),
                source_commit=commit, clean=True, source_inputs=source_inputs,
                source_selection='all tracked non-documentation C/C++/Python/XML/test/CMake inputs plus frozen contract/note',
                mandatory_originals={name: source_inputs[name] for name in REQUIRED},
                documentation_tools=documentation_tools, runtime_artifacts=dict(sorted(artifacts.items())),
                runtime_resolved_paths=dict(sorted(resolved.items())), runtime_bindings=dict(sorted(bindings.items())),
                runtime_environment=dict(sorted(runtime_environment.items())), python_packages=packages, tool_images=tools,
                effective_environment_policy='same fixed one-thread numerical/Python settings and geometry-lib overlay as the native qualification runner',
                ldd=ldd_output, host=dict(platform=platform.platform(), python=sys.version, executable=sys.executable,
                                        path=os.environ.get('PATH'), temporary_directory=os.environ.get('TMPDIR')),
                comparison_policy='all numerical source/runtime sets, hashes, dependency bindings, numerical env and package/tool code must match; documentation-only commits allowed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compareprevious', '--compare', dest='previous', type=Path)
    args = parser.parse_args()
    result = dict(schema=2, milestone='M6', success=False, clean=False,
                  captured_utc=datetime.now(timezone.utc).isoformat(), failures=[])
    try:
        result.update(capture(args.root.resolve()))
        if args.previous:
            previous = json.loads(args.previous.read_text())
            require(previous.get('success') and previous.get('clean'), 'Previous fingerprint did not succeed cleanly')
            differences = compare(previous, result)
            result.update(compared_fingerprint=str(args.previous.resolve()),
                          compared_fingerprint_sha256=digest(args.previous),
                          compared_source_commit=previous['source_commit'], differences=differences,
                          identical_qualification_inputs_and_runtime=not differences)
            require(not differences, 'Qualification input/runtime mismatch: '+json.dumps(differences, sort_keys=True))
        result['success'] = True
    except Exception as error:
        result['failures'].append(str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name+'.tmp-'+uuid.uuid4().hex)
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    temporary.replace(args.output)
    print(json.dumps(dict(success=result['success'], source_commit=result.get('source_commit'),
                          inputs=len(result.get('source_inputs', {})), runtime_files=len(result.get('runtime_artifacts', {})),
                          identical=result.get('identical_qualification_inputs_and_runtime'), failures=result['failures'])))
    return 0 if result['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
