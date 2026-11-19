#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Losslessly compact unbound failure reports from settled outlet attempts."""

import argparse
import gzip
import hashlib
import os
from pathlib import Path
import shutil

from resources import atomic, identity
from task2_outlet_baseline import read
from task2r_run import GIB, require, storage


def decoded_identity(path):
    digest = hashlib.sha256()
    count = 0
    with gzip.open(path, "rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
            count += len(block)
    return count,digest.hexdigest()


def referenced_paths(root, candidates):
    paths = set()
    for path in (root/"operations").glob("*/receipt.json"):
        receipt = path.read_text()
        paths.update(candidate for candidate in candidates if str(candidate) in receipt)
    return paths


def compact(root, labels, output):
    before = storage(root)
    sources = []
    require(len(labels)==len(set(labels)) and not output.exists(),"duplicate case or compaction result")
    for label in labels:
        require(Path(label).name==label and label.startswith("outlet-qual-"),"invalid named outlet attempt")
        receipt = root/"operations"/f"{label}-steady"/"receipt.json"
        failed = read(receipt)["result"]
        original = root/"science"/label/"state.failure.json"
        require(failed["status"]=="FAIL" and original.is_file(),
                "only a raw failure report from a settled attempt may be compacted")
        sources.append((receipt,original))
    references = referenced_paths(root,[path for _,path in sources])
    require(not references,f"registered inputs or outputs still refer to {sorted(references)}")
    require(before["science_bytes"]+sum(path.stat().st_size for _,path in sources)+1024*1024<=8*GIB,
            "temporary compressed copies may exceed the science storage limit")
    records = []
    for receipt,original in sources:
        archived = original.with_suffix(".json.gz")
        staged = archived.with_name(f"{archived.name}.pending-{os.getpid()}")
        observed = identity(original)
        if not archived.exists():
            try:
                with original.open("rb") as stream, staged.open("xb") as compressed:
                    with gzip.GzipFile(filename="",mode="wb",fileobj=compressed,
                                       mtime=0,compresslevel=1) as destination:
                        shutil.copyfileobj(stream,destination,1024*1024)
                    compressed.flush()
                    os.fsync(compressed.fileno())
                require(decoded_identity(staged)==(observed["bytes"],observed["sha256"])
                        and staged.stat().st_size<observed["bytes"],
                        "archive did not preserve the original failure bytes")
                staged.rename(archived)
                fd=os.open(archived.parent,os.O_RDONLY)
                try:
                    os.fsync(fd)
                finally:
                    os.close(fd)
            finally:
                staged.unlink(missing_ok=True)
        else:
            require(decoded_identity(archived)==(observed["bytes"],observed["sha256"]),
                    "existing compressed failure report differs from the original")
        records.append(dict(result=identity(receipt),original=observed,compressed=identity(archived)))
    atomic(output,dict(schema="WJH27FlatCampaignFailureCompactionV1",status="VERIFIED",
        records=records,prior_storage=before,
        restoration="Decode each .json.gz into a temporary file, verify SHA-256, then rename to the original path."))
    for row in records:
        Path(row["original"]["path"]).unlink()
    print("Compacted",len(records),"raw failure reports",
          before["science_bytes"]-storage(root)["science_bytes"],"bytes recovered",flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("labels",nargs="+")
    args = parser.parse_args()
    compact(args.root,args.labels,args.output)
