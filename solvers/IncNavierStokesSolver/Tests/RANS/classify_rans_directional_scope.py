#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import argparse,json
from pathlib import Path
from rans_directional_scope import classify
from rans_qualification import write_json


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('required-scope output already exists')
    try:value=classify(json.loads(a.input.read_text()))
    except Exception as error:value=dict(status='failed_unresolved_or_incomplete',failure=repr(error),response_extension_eligible=False)
    write_json(a.output,value);print(value['status']);return 0 if value['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
