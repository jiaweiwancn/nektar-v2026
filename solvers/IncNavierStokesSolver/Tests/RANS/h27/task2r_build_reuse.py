#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reuse verified object files that have no dependency on the R3 source changes."""

import argparse
from functools import lru_cache
import hashlib
from pathlib import Path
import re
import shlex
import shutil
import time

from resources import atomic, identity
from task2r_advancing import BINDING, read, require, verified


TARGETS = ("WJNativeMixedFlow", "IncNavierStokesSolver")


@lru_cache(maxsize=None)
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def uses_fingerprint(path):
    return b"WJ_H27_SOURCE_ID" in path.read_bytes()


def flags(path, source, build):
    flag_file = path / "flags.make"
    code = flag_file.read_text().replace(str(source), "<source>").replace(str(build), "<build>")
    code = re.sub(r'(WJ_H27_SOURCE_ID=\\?")[0-9a-f]{64}',r'\g<1><source-id>',code)
    return tuple(row for row in code.splitlines() if row.startswith(("CXX_DEFINES","CXX_INCLUDES","CXX_FLAGS")))


def compiler(build):
    cache = (build / "CMakeCache.txt").read_text()
    row, = (line for line in cache.splitlines() if line.startswith("CMAKE_CXX_COMPILER:"))
    return Path(row.split("=",1)[1]).resolve(strict=True)


def reuse(root, out):
    began = time.monotonic()
    basis = read(root/"evidence"/BINDING)
    verified(basis["x15_parity"])
    old_build, new_build = root/"build-outlet-campaign-extent", root/"build-r3"
    old_source, new_source = root/"source-outlet-campaign-extent", root/"source-r3"
    driver = compiler(new_build)
    require(driver == compiler(old_build), "configured Release compiler differs")
    seal = read(root/"evidence/outlet-campaign-extent-runtime.json")
    executable = verified(seal["application"])
    linked = old_build/"IncNavierStokesSolver"
    require(digest(linked) == digest(executable),
            "build archive does not contain the verified x=15 replay application")
    linked_at = linked.stat().st_mtime_ns
    records = {}
    for target in TARGETS:
        sub = Path("CMakeFiles")/f"{target}.dir"
        old, new = old_build/sub, new_build/sub
        require(flags(old,old_source,old_build) == flags(new,new_source,new_build),
                f"{target} compiler flags differ beyond the declared source fingerprint")
        makefile = (new/"build.make").read_text()
        records[target] = entries = []
        for dependencies in sorted(old.rglob("*.o.d")):
            obj = dependencies.with_suffix("")
            rel = obj.relative_to(old_build)
            mapped = Path(str(rel).replace(str(old_source).lstrip("/"),str(new_source).lstrip("/")))
            dest = new_build / mapped
            if not obj.is_file() or dest.exists() or str(mapped)+":" not in makefile:
                continue
            text = dependencies.read_text()
            _, _, names = text.replace("\\\n"," ").partition(":")
            dep_files = [Path(name) for name in shlex.split(names)]
            if not dep_files or not any(p.is_relative_to(old_source) for p in dep_files):
                continue
            before = obj.stat().st_mtime_ns
            with obj.open("rb") as stream:
                elf = stream.read(4) == b"\x7fELF"
            if not elf or not driver.stat().st_mtime_ns <= before <= linked_at:
                continue
            used = []
            reusable = True
            for entry in dep_files:
                if (not entry.is_absolute() or not entry.is_file() or entry.stat().st_mtime_ns > before
                        or uses_fingerprint(entry)):
                    reusable = False
                    break
                if entry.is_relative_to(old_source) or entry.is_relative_to(old_build):
                    parent = old_source if entry.is_relative_to(old_source) else old_build
                    replacement = (new_source if parent == old_source else new_build)/entry.relative_to(parent)
                    if not replacement.is_file() or digest(entry) != digest(replacement):
                        reusable = False
                        break
                    used.append(dict(path=str(replacement),sha256=digest(replacement)))
            if not reusable:
                continue
            dest.parent.mkdir(parents=True,exist_ok=True)
            obj_tmp = dest.with_suffix(".o.r3-stage")
            dep_tmp = dest.with_suffix(".o.d.r3-stage")
            shutil.copyfile(obj,obj_tmp)
            dep_tmp.write_text(text.replace(str(old_source),str(new_source)).replace(str(old_build),str(new_build)))
            require(digest(obj) == digest(obj_tmp),"reused compile object changed during copy")
            dep_tmp.replace(dest.with_suffix(".o.d"))
            obj_tmp.replace(dest)
            entries.append(dict(source=identity(obj),reused=identity(dest),local_dependencies=used))
    require(any(records.values()),"no verified compile objects fit the new target inventory")
    atomic(out,dict(schema="WJH27R3VerifiedBuildObjectReuseV1",application=seal["application"],
        prior_replay=identity(verified(basis["x15_parity"])),compiler=identity(driver),records=records,
        elapsed_seconds=time.monotonic()-began,diagnostic_path_note="Reused objects retain archived diagnostic source locations. "
        "Changed sources, changed headers and units consuming the scientific fingerprint compile inside the new build."))
    print({target:len(rows) for target,rows in records.items()})


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("out",type=Path)
    reuse(**vars(parser.parse_args()))
